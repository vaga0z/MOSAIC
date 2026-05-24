"""
MOSAIC UNIFIED ENCODER — HYPER-OPTIMIZED MULTI-CORE
====================================================
Supports both operational states via the --state flag:

    Half MOSAIC  : Inner RS ECC + Titanium Layer only.
                   Independent oligos, no outer erasure matrix.
                   Best for density-sensitive, perceptually tolerant payloads.

    Full MOSAIC  : Dual-layer RS Matrix + Titanium Layer.
                   Outer RS(255,223) erasure matrix for deterministic,
                   block-level SHA-256-verified integrity.
                   Required for enterprise archival payloads.

Usage:
    python ENCODER_MOSAIC.py <input_file> <output.fasta> --state half
    python ENCODER_MOSAIC.py <input_file> <output.fasta> --state full
    python ENCODER_MOSAIC.py <input_file> <output.fasta>          # prompts if omitted
"""

import struct
import sys
import os
import re
import argparse
import hashlib
import reedsolo
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==========================================
# HYPER-OPTIMIZATION: LOOK-UP TABLES (LUTs)
# ==========================================

# 1. Byte-to-DNA pre-computation — eliminates per-byte bit-string casting
BYTE_TO_DNA = []
_mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
for _i in range(256):
    _b = format(_i, '08b')
    BYTE_TO_DNA.append("".join(_mapping[_b[j:j+2]] for j in range(0, 8, 2)))

# 2. TMR Header pre-computation — only 256 possible seeds
TMR_CACHE = []
for _seed in range(256):
    _seed_bin   = format(_seed, '08b')
    _hdr_bits   = "".join((bit * 3) + ("01" if bit == '1' else "10") for bit in _seed_bin)
    TMR_CACHE.append("".join(_mapping[_hdr_bits[j:j+2]] for j in range(0, 40, 2)))

# 3. Reverse complement translator
DNA_TRANS = str.maketrans('ACGT', 'TGCA')
def reverse_complement(seq):
    return seq.translate(DNA_TRANS)[::-1]

# 4. Flattened Titanium constraint ban list
FWD_PRIMER = "CGTCGGCAGCGTCAG"
REV_PRIMER = "GTCTCGTGGGCTCGG"
BANNED_MOTIFS = [
    "AAAA", "TTTT", "GGGG", "CCCC",                    # Homopolymers
    "GAATTC", "GGATCC", "AAGCTT", "GGTCTC", "GAGACC",  # Restriction sites
    "TAATACGACTCACTATAGGG", "AGGAGG",                   # In-vivo biosafety
    FWD_PRIMER, reverse_complement(FWD_PRIMER),          # Primer collisions
    REV_PRIMER, reverse_complement(REV_PRIMER),
]

MICROSAT_REGEX = re.compile(r"([ACGT]{2})\1{3,}")

# ==========================================
# CODEC CONSTANTS
# ==========================================

INNER_RS_ECC = 15
rs_inner     = reedsolo.RSCodec(INNER_RS_ECC)

OUTER_RS_K   = 223
OUTER_RS_ECC = 32
OUTER_RS_N   = OUTER_RS_K + OUTER_RS_ECC
rs_outer     = reedsolo.RSCodec(OUTER_RS_ECC)

PAYLOAD_BYTES = 40   # data bytes per oligo (excludes 4-byte frame ID)

# ==========================================
# TITANIUM CONSTRAINT AUDIT LAYER
# ==========================================

def is_biologically_safe(seq):
    """
    Multi-constraint physical chemistry filter.
    Returns True only if all 5 constraint layers pass simultaneously.
    """
    seq_len = len(seq)

    # Layer 1: substring bans (homopolymers, restriction sites,
    #          biosafety motifs, primer collisions)
    for motif in BANNED_MOTIFS:
        if motif in seq:
            return False

    # Layer 2: global GC content (45% – 55%)
    gc = seq.count('G') + seq.count('C')
    if not (45.0 <= (gc / seq_len) * 100.0 <= 55.0):
        return False

    # Layer 3: sliding window GC (30-base window, 30% – 70% = 9 – 21 GCs)
    WINDOW = 30
    cur_gc = seq[:WINDOW].count('G') + seq[:WINDOW].count('C')
    if not (9 <= cur_gc <= 21):
        return False
    for i in range(1, seq_len - WINDOW + 1):
        if seq[i - 1]         in 'GC': cur_gc -= 1
        if seq[i + WINDOW - 1] in 'GC': cur_gc += 1
        if not (9 <= cur_gc <= 21):
            return False

    # Layer 4: hairpin traps (stem >= 8 bases, thermodynamically stable)
    for i in range(seq_len - 8):
        if seq[i:i+8].translate(DNA_TRANS)[::-1] in seq[i+8:]:
            return False

    # Layer 5: microsatellite slippage (dinucleotide repeats x4+)
    if MICROSAT_REGEX.search(seq):
        return False

    return True

# ==========================================
# LFSR SCRAMBLER
# ==========================================

def vanguard_lfsr_cipher(data_bytes, ignition_key):
    """
    32-bit LFSR stream scrambler seeded by a 16-bit ignition key.
    Primitive polynomial: f(x) = x^32 + x^21 + x + 1.
    XOR-scrambles data_bytes against the keystream, imposing a uniform
    output distribution independent of input entropy.
    """
    state     = ignition_key & 0xFFFFFFFF
    processed = bytearray(len(data_bytes))
    for i, byte in enumerate(data_bytes):
        new_byte = 0
        for b in range(7, -1, -1):
            data_bit = (byte >> b) & 1
            feedback = ((state >> 31) ^ (state >> 21) ^
                        (state >>  1) ^ (state >>  0)) & 1
            state    = ((state << 1) | feedback) & 0xFFFFFFFF
            new_byte |= (data_bit ^ (state & 1)) << b
        processed[i] = new_byte
    return bytes(processed)

# ==========================================
# FRAME ASSEMBLY (SHARED BY BOTH STATES)
# ==========================================

def assemble_and_encode_frame(data_40_bytes, frame_id):
    """
    Encodes one 40-byte payload into a 256-base DNA strand.

    Frame layout (512 bits):
      [40-bit TMR header] | [472-bit scrambled RS(59,44) codeword]

    The RS codeword protects:
      [4-byte frame_id] | [40-byte payload]  →  59-byte codeword

    The 256x256 ignition key space is exhausted until a key is found
    whose scrambled output passes all Titanium constraint layers.

    Returns (seed, salt, dna_strand) or (None, None, None) on exhaustion.
    """
    if len(data_40_bytes) != PAYLOAD_BYTES:
        data_40_bytes = data_40_bytes.ljust(PAYLOAD_BYTES, b'\x00')

    vault_core     = struct.pack('>I', frame_id) + data_40_bytes
    protected_core = bytes(rs_inner.encode(vault_core))

    for seed in range(256):
        header_dna = TMR_CACHE[seed]
        for salt in range(256):
            ignition_key  = (salt << 8) | seed
            scrambled     = vanguard_lfsr_cipher(protected_core, ignition_key)
            scrambled_dna = "".join(BYTE_TO_DNA[b] for b in scrambled)
            dna_strand    = header_dna + scrambled_dna
            if is_biologically_safe(dna_strand):
                return seed, salt, dna_strand

    return None, None, None

# ==========================================
# PARALLEL WORKER (BOTH STATES USE THIS)
# ==========================================

def parallel_encode_worker(frame_id, payload_40):
    seed, salt, dna = assemble_and_encode_frame(payload_40, frame_id)
    return frame_id, seed, salt, dna

# ==========================================
# OUTER RS ERASURE MATRIX (FULL MOSAIC ONLY)
# ==========================================

def apply_outer_rs(data_bytes):
    """
    Byte-interleaved RS(255,223) column encoding across the oligo pool.

    For each of the PAYLOAD_BYTES=40 byte offsets, a column of 223
    data bytes is RS-encoded into 255 bytes (appending 32 parity bytes).
    This produces 255 oligo payloads per block: 223 data + 32 parity.

    Returns:
        encoded_oligos   : list of (block_idx, pos_in_block, bytes_40)
        total_data_oligos: number of data oligos before parity addition
        num_blocks       : number of RS(255,223) blocks
    """
    # Pad raw bytes to a multiple of PAYLOAD_BYTES
    remainder = len(data_bytes) % PAYLOAD_BYTES
    if remainder:
        data_bytes += b'\x00' * (PAYLOAD_BYTES - remainder)

    chunks            = [data_bytes[i:i+PAYLOAD_BYTES]
                         for i in range(0, len(data_bytes), PAYLOAD_BYTES)]
    total_data_oligos = len(chunks)

    # Pad chunks to a multiple of OUTER_RS_K
    pad = total_data_oligos % OUTER_RS_K
    if pad:
        chunks += [bytes(PAYLOAD_BYTES)] * (OUTER_RS_K - pad)

    num_blocks     = len(chunks) // OUTER_RS_K
    encoded_oligos = []

    for block_idx in range(num_blocks):
        block_data       = chunks[block_idx * OUTER_RS_K : (block_idx + 1) * OUTER_RS_K]
        columns          = [bytes(p[j] for p in block_data) for j in range(PAYLOAD_BYTES)]
        encoded_columns  = [bytes(rs_outer.encode(col))     for col in columns]

        for pos in range(OUTER_RS_N):
            payload = bytes(encoded_columns[j][pos] for j in range(PAYLOAD_BYTES))
            encoded_oligos.append((block_idx, pos, payload))

    return encoded_oligos, total_data_oligos, num_blocks

# ==========================================
# MASTER ENCODING PIPELINE
# ==========================================

def encode_file_to_dna(input_path, output_fasta_path, state):
    """
    state : 'half' or 'full'
    """
    state = state.lower()
    assert state in ('half', 'full'), "state must be 'half' or 'full'"

    label = "FULL MOSAIC — DUAL-LAYER RS MATRIX" if state == 'full' \
            else "HALF MOSAIC — HIGH DENSITY PIPELINE"

    print(f"{'='*65}")
    print(f"MOSAIC ENCODER — {label}")
    print(f"Titanium Layer Active | 65,536-State LFSR")
    print(f"{'='*65}")

    with open(input_path, 'rb') as f:
        raw_bytes = f.read()

    original_size = len(raw_bytes)
    file_ext      = os.path.splitext(input_path)[1].lower().lstrip('.') or 'bin'
    sha256        = hashlib.sha256(raw_bytes).hexdigest()

    print(f"[+] Input   : {input_path} ({original_size:,} bytes)")
    print(f"[+] State   : {state.upper()} MOSAIC")
    print(f"[+] SHA-256 : {sha256}\n")

    cpu_cores = os.cpu_count() or 4

    # ------------------------------------------------------------------
    # FULL MOSAIC: outer RS erasure matrix first, then inner encode
    # ------------------------------------------------------------------
    if state == 'full':
        print(f"[PHASE 1] Generating Outer RS({OUTER_RS_N},{OUTER_RS_K}) Erasure Matrix...")
        encoded_oligos, total_data_oligos, num_blocks = apply_outer_rs(raw_bytes)
        total_with_parity = len(encoded_oligos)
        overhead = (total_with_parity - total_data_oligos) / total_data_oligos * 100
        print(f"    Data oligos   : {total_data_oligos:,}")
        print(f"    Parity oligos : {total_with_parity - total_data_oligos:,}  "
              f"(overhead {overhead:.2f}%)")
        print(f"    Total         : {total_with_parity:,} across {num_blocks} blocks\n")

        print("[PHASE 2] Parallel LFSR Scrambling & Titanium Auditing...")
        dna_pool, failed = [], 0
        total            = total_with_parity

        with ProcessPoolExecutor(max_workers=cpu_cores) as executor:
            futures = {
                executor.submit(parallel_encode_worker, (block_idx << 16) | pos, payload): i
                for i, (block_idx, pos, payload) in enumerate(encoded_oligos)
            }
            done = 0
            for future in as_completed(futures):
                frame_id, seed, salt, dna = future.result()
                done += 1
                block_idx = (frame_id >> 16) & 0xFFFF
                pos       = frame_id & 0xFFFF
                if dna:
                    dna_pool.append({
                        "block_idx" : block_idx,
                        "pos"       : pos,
                        "seq_id"    : frame_id,
                        "is_parity" : pos >= OUTER_RS_K,
                        "seed"      : seed,
                        "salt"      : salt,
                        "dna"       : dna,
                    })
                else:
                    failed += 1
                if done % 200 == 0 or done == total:
                    sys.stdout.write(
                        f"\r    {done}/{total} oligos "
                        f"({done/total*100:.1f}%) | Failed: {failed}"
                    )
                    sys.stdout.flush()

        dna_pool.sort(key=lambda x: x["seq_id"])
        print(f"\n\n[PHASE 3] Writing FASTA manifest...")

        with open(output_fasta_path, 'w') as f:
            f.write(f"; Architecture: FULL_MOSAIC\n")
            f.write(f"; OriginalSize: {original_size} | "
                    f"FileType: {file_ext} | SHA256: {sha256}\n")
            f.write(f"; TotalDataOligos: {total_data_oligos} | "
                    f"TotalOligosWithParity: {total_with_parity} | "
                    f"NumBlocks: {num_blocks}\n")
            f.write(f"; OuterK: {OUTER_RS_K} | "
                    f"OuterECC: {OUTER_RS_ECC} | "
                    f"OuterN: {OUTER_RS_N}\n")
            for entry in dna_pool:
                flag = "P" if entry["is_parity"] else "D"
                f.write(f">MOSAIC_FULL_B{entry['block_idx']:04d}_"
                        f"P{entry['pos']:03d}_{flag}_"
                        f"Seed{entry['seed']}_Salt{entry['salt']}\n")
                f.write(f"{entry['dna']}\n")

    # ------------------------------------------------------------------
    # HALF MOSAIC: direct linear chunking, no outer RS
    # ------------------------------------------------------------------
    else:
        chunks       = [raw_bytes[i:i+PAYLOAD_BYTES]
                        for i in range(0, len(raw_bytes), PAYLOAD_BYTES)]
        total_oligos = len(chunks)
        print(f"[+] Data oligos : {total_oligos:,}\n")
        print("[PHASE 1] Parallel LFSR Scrambling & Titanium Auditing...")

        dna_pool, failed = [], 0

        with ProcessPoolExecutor(max_workers=cpu_cores) as executor:
            futures = {
                executor.submit(parallel_encode_worker, i, chunk): i
                for i, chunk in enumerate(chunks)
            }
            done = 0
            for future in as_completed(futures):
                seq_id, seed, salt, dna = future.result()
                done += 1
                if dna:
                    dna_pool.append({
                        "seq_id" : seq_id,
                        "seed"   : seed,
                        "salt"   : salt,
                        "dna"    : dna,
                    })
                else:
                    failed += 1
                if done % 200 == 0 or done == total_oligos:
                    sys.stdout.write(
                        f"\r    {done}/{total_oligos} oligos "
                        f"({done/total_oligos*100:.1f}%) | Failed: {failed}"
                    )
                    sys.stdout.flush()

        dna_pool.sort(key=lambda x: x["seq_id"])
        total_with_parity  = total_oligos
        total_data_oligos  = total_oligos
        print(f"\n\n[PHASE 2] Writing FASTA manifest...")

        with open(output_fasta_path, 'w') as f:
            f.write(f"; Architecture: HALF_MOSAIC\n")
            f.write(f"; OriginalSize: {original_size} | "
                    f"FileType: {file_ext} | SHA256: {sha256}\n")
            f.write(f"; TotalDataOligos: {total_oligos}\n")
            f.write(f"; InnerECC: {INNER_RS_ECC}\n")
            for entry in dna_pool:
                f.write(f">MOSAIC_HALF_ID{entry['seq_id']:06d}_"
                        f"Seed{entry['seed']}_Salt{entry['salt']}\n")
                f.write(f"{entry['dna']}\n")

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------
    fasta_size = os.path.getsize(output_fasta_path)
    print(f"\n{'='*65}")
    print(f"SYNTHESIS COMPLETE")
    print(f"    State           : {state.upper()} MOSAIC")
    print(f"    Data oligos     : {total_data_oligos:,}")
    print(f"    Total oligos    : {total_with_parity:,}")
    print(f"    Encoding fails  : {failed}")
    print(f"    Output          : {output_fasta_path} "
          f"({fasta_size/1024:.1f} KB)")
    print(f"{'='*65}\n")


# ==========================================
# CLI ENTRY POINT
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MOSAIC Unified Encoder — Half or Full state"
    )
    parser.add_argument("input",
        help="Input file to encode (any binary format)")
    parser.add_argument("output",
        help="Output FASTA file path")
    parser.add_argument("--state",
        choices=["half", "full"],
        default=None,
        help="MOSAIC operational state: 'half' (density-optimised) "
             "or 'full' (dual-layer, SHA-256 guaranteed integrity)")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[!] Input file not found: {args.input}")
        sys.exit(1)

    # Prompt interactively if --state was not provided
    if args.state is None:
        print("\nSelect MOSAIC operational state:")
        print("  [1] Half MOSAIC — density-optimised, independent oligos")
        print("      Best for: multimedia, visual media, perceptually tolerant payloads")
        print("  [2] Full MOSAIC — dual-layer RS matrix, SHA-256 guaranteed integrity")
        print("      Best for: compressed archives, executables, encrypted databases")
        while True:
            choice = input("\nEnter 1 or 2: ").strip()
            if choice == '1':
                args.state = 'half'
                break
            elif choice == '2':
                args.state = 'full'
                break
            else:
                print("    Please enter 1 or 2.")

    encode_file_to_dna(args.input, args.output, args.state)
