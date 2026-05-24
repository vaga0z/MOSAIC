"""
VANGUARD — UNIFIED MOSAIC DECODER
===================================
Auto-detects architecture from FASTA metadata and routes to the correct
decoding pipeline.  A single script replaces both DECODER_HALF_MOSAIC.py
and DECODER_FULL_MOSAIC.py.

Architecture detection
----------------------
The decoder reads the '; Architecture: ...' comment line at the top of the
FASTA file before touching any oligo data.

    ; Architecture: HALF_MOSAIC   →  V3 single-layer inner-RS pipeline
    ; Architecture: FULL_MOSAIC   →  V4 two-layer inner + outer-RS pipeline

HALF_MOSAIC FASTA header example
---------------------------------
; Architecture: HALF_MOSAIC
; OriginalSize: 1000000 | FileType: png | SHA256: a3f9...
; TotalDataOligos: 25000
; InnerECC: 15

>MOSAIC_HALF_ID000000_Seed12_Salt201
ACGT...

FULL_MOSAIC FASTA header example
---------------------------------
; Architecture: FULL_MOSAIC
; OriginalSize: 1000000 | FileType: png | SHA256: a3f9...
; TotalDataOligos: 25000 | TotalOligosWithParity: 28560 | NumBlocks: 112
; OuterK: 223 | OuterECC: 32 | OuterN: 255

>MOSAIC_FULL_B0000_P000_D_Seed12_Salt201
ACGT...

Usage
-----
    python DECODER_MOSAIC_UNIFIED.py <input.fasta> [options]
    python DECODER_MOSAIC_UNIFIED.py *.fasta --batch --coverage-range 3 10

Options
-------
    --dropout FLOAT       Fraction of oligos to simulate as physically lost (0.0-1.0)
    --coverage INT        Sequencing coverage depth per oligo (default: 5)
    --coverage-range MIN MAX   Sweep coverage range (for --batch mode)
    --model {R9.4,R10.4}  Nanopore error profile (default: R10.4)
    --output PATH         Output filename (auto-derived from metadata if omitted)
    --stress              Run dropout stress sweep (0%–30%)
    --batch               Run multi-FASTA coverage sweep
"""

import struct
import sys
import os
import random
import time
import hashlib
import argparse
import glob
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed

import reedsolo


# ==========================================
# SHARED CODEC CONSTANTS — MATCH ENCODERS
# ==========================================

INNER_RS_ECC  = 15
PAYLOAD_BYTES = 40          # 40 payload bytes per oligo (both architectures)

rs_inner = reedsolo.RSCodec(INNER_RS_ECC)

# FULL_MOSAIC outer RS defaults (overridden dynamically from metadata)
OUTER_RS_K   = 223
OUTER_RS_ECC = 32
OUTER_RS_N   = OUTER_RS_K + OUTER_RS_ECC
rs_outer     = reedsolo.RSCodec(OUTER_RS_ECC)

# Nanopore error profiles
ERROR_PROFILES = {
    "R9.4":  {"sub": 0.030, "ins": 0.010, "del": 0.020},
    "R10.4": {"sub": 0.005, "ins": 0.002, "del": 0.004},
}


# ==========================================
# ARCHITECTURE DETECTION
# ==========================================

def detect_architecture(fasta_path):
    """
    Fast pre-scan: reads only the ';' comment lines at the top of the FASTA
    file and returns 'HALF_MOSAIC', 'FULL_MOSAIC', or raises ValueError if
    the Architecture field is missing or unrecognised.
    """
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>') or (not line.startswith(';')):
                # Reached oligo data without finding Architecture
                break
            if 'Architecture' in line:
                # e.g. '; Architecture: FULL_MOSAIC'
                _, _, value = line.partition('Architecture:')
                arch = value.strip().split()[0].upper()
                if arch in ('HALF_MOSAIC', 'FULL_MOSAIC'):
                    return arch
                raise ValueError(
                    f"Unknown Architecture value '{arch}' in {fasta_path}"
                )

    raise ValueError(
        f"No 'Architecture:' field found in FASTA metadata of {fasta_path}. "
        "Expected '; Architecture: HALF_MOSAIC' or '; Architecture: FULL_MOSAIC'."
    )


# ==========================================
# SHARED: BIOLOGICAL SIMULATION & CONSENSUS
# ==========================================

def nanopore_markov_simulator(clean_dna, sub_rate, ins_rate, del_rate,
                               homopolymer_multiplier=2.0):
    """Simulate Nanopore sequencing errors on a single clean DNA strand."""
    noisy, bases      = "", ['A', 'C', 'G', 'T']
    consecutive, prev = 1, ""
    for base in clean_dna:
        if base == prev: consecutive += 1
        else:            consecutive  = 1
        prev    = base
        eff_del = del_rate * homopolymer_multiplier if consecutive >= 3 else del_rate
        r       = random.random()
        if   r < eff_del:                         continue
        elif r < eff_del + sub_rate:              noisy += random.choice([b for b in bases if b != base])
        elif r < eff_del + sub_rate + ins_rate:   noisy += base + random.choice(bases)
        else:                                     noisy += base
    return noisy


BAND    = 40
NEG_INF = float('-inf')


def needleman_wunsch_align(seq1, seq2, band=BAND):
    """
    Banded Needleman-Wunsch alignment.
    Falls back to full O(n²) NW when the length delta exceeds the band.
    """
    n, m = len(seq1), len(seq2)
    if abs(n - m) > band:
        return _needleman_wunsch_full(seq1, seq2)

    score = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    for i in range(min(n + 1, band + 1)): score[i][0] = -i
    for j in range(min(m + 1, band + 1)): score[0][j] = -j

    for i in range(1, n + 1):
        j_lo = max(1, i - band)
        j_hi = min(m, i + band)
        for j in range(j_lo, j_hi + 1):
            diag = score[i-1][j-1] if score[i-1][j-1] != NEG_INF else NEG_INF
            up   = score[i-1][j]   if score[i-1][j]   != NEG_INF else NEG_INF
            left = score[i][j-1]   if score[i][j-1]   != NEG_INF else NEG_INF
            match_score = diag + (1 if seq1[i-1] == seq2[j-1] else -1)
            score[i][j] = max(match_score, up - 1, left - 1)

    a1, a2, i, j = "", "", n, m
    while i > 0 and j > 0:
        cur      = score[i][j]
        diag_val = score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1)
        if cur == diag_val:
            a1 += seq1[i-1]; a2 += seq2[j-1]; i -= 1; j -= 1
        elif score[i-1][j] != NEG_INF and cur == score[i-1][j] - 1:
            a1 += seq1[i-1]; a2 += '-'; i -= 1
        else:
            a1 += '-'; a2 += seq2[j-1]; j -= 1
    while i > 0: a1 += seq1[i-1]; a2 += '-'; i -= 1
    while j > 0: a1 += '-'; a2 += seq2[j-1]; j -= 1
    return a1[::-1], a2[::-1]


def _needleman_wunsch_full(seq1, seq2):
    """Full O(n²) NW — fallback only."""
    n, m  = len(seq1), len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1): score[i][0] = -i
    for j in range(m + 1): score[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match       = score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1)
            score[i][j] = max(match, score[i-1][j] - 1, score[i][j-1] - 1)
    a1, a2, i, j = "", "", n, m
    while i > 0 and j > 0:
        cur = score[i][j]
        if cur == score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1):
            a1 += seq1[i-1]; a2 += seq2[j-1]; i -= 1; j -= 1
        elif cur == score[i-1][j] - 1:
            a1 += seq1[i-1]; a2 += '-'; i -= 1
        else:
            a1 += '-'; a2 += seq2[j-1]; j -= 1
    while i > 0: a1 += seq1[i-1]; a2 += '-'; i -= 1
    while j > 0: a1 += '-'; a2 += seq2[j-1]; j -= 1
    return a1[::-1], a2[::-1]


def build_consensus(reads, iterations=1):
    """Multi-read consensus via iterative Needleman-Wunsch voting."""
    if not reads: return ""
    reads.sort(key=len)
    master = reads[len(reads) // 2]
    for _ in range(iterations):
        votes      = [{'A':0,'C':0,'G':0,'T':0,'-':0} for _ in range(len(master))]
        insertions = [{} for _ in range(len(master) + 1)]
        for r in reads:
            a_al, r_al  = needleman_wunsch_align(master, r)
            m_idx, cur_ins = 0, ""
            for a_c, r_c in zip(a_al, r_al):
                if a_c == '-':
                    if r_c != '-': cur_ins += r_c
                else:
                    if cur_ins:
                        insertions[m_idx][cur_ins] = insertions[m_idx].get(cur_ins, 0) + 1
                        cur_ins = ""
                    votes[m_idx][r_c] += 1
                    m_idx += 1
            if cur_ins:
                insertions[m_idx][cur_ins] = insertions[m_idx].get(cur_ins, 0) + 1
        new_master = ""
        for i in range(len(master)):
            if insertions[i]:
                best = max(insertions[i], key=insertions[i].get)
                if insertions[i][best] >= len(reads) * 0.4: new_master += best
            best_base = max(votes[i], key=votes[i].get)
            if best_base != '-': new_master += best_base
        if insertions[-1]:
            best = max(insertions[-1], key=insertions[-1].get)
            if insertions[-1][best] >= len(reads) * 0.4: new_master += best
        master = new_master
    return master


# ==========================================
# SHARED: BIT / BYTE PRIMITIVES
# ==========================================

def dna_to_binary(dna):
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    try:    return "".join(mapping[b] for b in dna)
    except: return None


def bits_to_bytes(bits):
    return bytearray(int(bits[i:i+8], 2) for i in range(0, len(bits), 8))


def vanguard_lfsr_cipher(data_bytes, ignition_key):
    """32-bit LFSR stream cipher — identical in both architectures."""
    state     = ignition_key & 0xFFFFFFFF
    processed = bytearray()
    for byte in data_bytes:
        new_byte = 0
        for i in range(7, -1, -1):
            data_bit   = (byte >> i) & 1
            feedback   = ((state >> 31) ^ (state >> 21) ^ (state >> 1) ^ (state >> 0)) & 1
            state      = ((state << 1) | feedback) & 0xFFFFFFFF
            cipher_bit = state & 1
            new_byte  |= (data_bit ^ cipher_bit) << i
        processed.append(new_byte)
    return bytes(processed)


# ==========================================
# SHARED: FASTA METADATA PARSER
# ==========================================

def _parse_metadata_lines(fasta_path):
    """
    Returns a dict of all key:value pairs from the ';' comment lines.
    Stops parsing at the first non-comment, non-blank line.
    """
    metadata = {}
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not line.startswith(';'):
                break
            parts = line.lstrip(';').strip().split('|')
            for part in parts:
                part = part.strip()
                if ':' in part:
                    key, val = part.split(':', 1)
                    metadata[key.strip()] = val.strip()
    return metadata


# ==========================================
# HALF_MOSAIC: FASTA PARSER
# ==========================================

def parse_fasta_half(fasta_path):
    """
    Returns:
        metadata   : dict
        oligo_list : list of (frame_id: int, dna: str)

    Header format: >MOSAIC_HALF_ID000123_Seed55_Salt201
    """
    metadata   = _parse_metadata_lines(fasta_path)
    oligo_list = []

    with open(fasta_path, 'r') as f:
        current_header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(';'):
                continue
            if line.startswith('>'):
                current_header = line[1:]
            else:
                if current_header:
                    try:
                        parts    = current_header.split('_')
                        id_part  = next(p for p in parts if p.startswith('ID') and p[2:].isdigit())
                        frame_id = int(id_part[2:])
                        oligo_list.append((frame_id, line))
                    except (IndexError, ValueError, StopIteration):
                        oligo_list.append((0, line))
                    current_header = None

    return metadata, oligo_list


# ==========================================
# FULL_MOSAIC: FASTA PARSER
# ==========================================

def parse_fasta_full(fasta_path):
    """
    Returns:
        metadata   : dict
        oligo_list : list of (block_idx: int, pos: int, dna: str)

    Header format: >MOSAIC_FULL_B0000_P000_D_Seed12_Salt201
                   >MOSAIC_FULL_B0000_P223_P_Seed88_Salt004
    """
    metadata   = _parse_metadata_lines(fasta_path)
    oligo_list = []

    with open(fasta_path, 'r') as f:
        current_header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(';'):
                continue
            if line.startswith('>'):
                current_header = line[1:]
            else:
                if current_header:
                    try:
                        parts      = current_header.split('_')
                        block_part = next(p for p in parts if p.startswith('B') and p[1:].isdigit())
                        pos_part   = next(p for p in parts if p.startswith('P') and p[1:].isdigit())
                        block_idx  = int(block_part[1:])
                        pos        = int(pos_part[1:])
                        oligo_list.append((block_idx, pos, line))
                    except (IndexError, ValueError, StopIteration):
                        oligo_list.append((0, 0, line))
                    current_header = None

    return metadata, oligo_list


# ==========================================
# HALF_MOSAIC: INNER FRAME DECODER  (V3)
# ==========================================

def decode_frame_half(dna_molecule, expected_frame_id):
    """
    Decodes a HALF_MOSAIC oligo.

    Returns: (frame_id, payload_bytes, success)
    On failure uses the 'Analog Degradation' fallback — corrupted payload
    bytes are passed through to preserve the V3 ablation characteristic.
    """
    binary = dna_to_binary(dna_molecule)
    if not binary or len(binary) != 512:
        return expected_frame_id, os.urandom(PAYLOAD_BYTES), False

    # Recover 8-bit seed from 40-bit TMR interleaved header
    header_bits = binary[0:40]
    seed_bits   = ""
    for i in range(0, 40, 5):
        tmr_block  = header_bits[i:i+3]
        seed_bits += '1' if tmr_block.count('1') > 1 else '0'
    seed = int(seed_bits, 2)

    vault_bytes = bits_to_bytes(binary[40:512])   # 59 bytes

    for salt in range(256):
        ignition_key = (salt << 8) | seed
        decrypted    = vanguard_lfsr_cipher(vault_bytes, ignition_key)
        try:
            decoded_core = bytes(rs_inner.decode(decrypted)[0])
            frame_id     = struct.unpack('>I', decoded_core[0:4])[0]
            payload      = decoded_core[4:4+PAYLOAD_BYTES]
            return frame_id, payload, True
        except reedsolo.ReedSolomonError:
            continue

    # Analog Degradation fallback — biological noise passes through
    decrypted_garbage = vanguard_lfsr_cipher(vault_bytes, 0)
    noisy_payload     = decrypted_garbage[4:4+PAYLOAD_BYTES]
    noisy_payload     = noisy_payload.ljust(PAYLOAD_BYTES, b'\xFF')[:PAYLOAD_BYTES]
    return expected_frame_id, noisy_payload, False


# ==========================================
# FULL_MOSAIC: INNER FRAME DECODER  (V4)
# ==========================================

def decode_frame_full(dna_molecule):
    """
    Decodes a FULL_MOSAIC oligo.

    Returns: (payload_bytes, block_idx, pos, success, mutations_corrected)
    On failure returns (None, None, None, False, 0) — treated as erasure
    for the outer RS layer.
    """
    binary = dna_to_binary(dna_molecule)
    if not binary or len(binary) != 512:
        return None, None, None, False, 0

    # Recover 8-bit seed from 40-bit TMR interleaved header
    header_bits = binary[0:40]
    seed_bits   = ""
    for i in range(0, 40, 5):
        tmr_block  = header_bits[i:i+3]
        seed_bits += '1' if tmr_block.count('1') > 1 else '0'
    seed = int(seed_bits, 2)

    vault_bytes = bits_to_bytes(binary[40:512])   # 59 bytes

    for salt in range(256):
        ignition_key = (salt << 8) | seed
        decrypted    = vanguard_lfsr_cipher(vault_bytes, ignition_key)
        try:
            decoded_core = bytes(rs_inner.decode(decrypted)[0])   # 44 bytes

            # seq_id encodes (block_idx << 16) | pos
            seq_id    = struct.unpack('>I', decoded_core[0:4])[0]
            payload   = decoded_core[4:4+PAYLOAD_BYTES]
            block_idx = (seq_id >> 16) & 0xFFFF
            pos       = seq_id & 0xFFFF

            repaired  = bytearray(rs_inner.encode(decoded_core))
            mutations = sum(1 for a, b in zip(decrypted, repaired) if a != b)

            return payload, block_idx, pos, True, mutations
        except reedsolo.ReedSolomonError:
            continue

    return None, None, None, False, 0


# ==========================================
# PARALLEL WORKERS  (must be module-level)
# ==========================================

def _worker_half(clean_dna, config, expected_frame_id):
    """ProcessPoolExecutor worker for HALF_MOSAIC."""
    noisy_pool = [
        nanopore_markov_simulator(clean_dna, config['sub'], config['ins'], config['del'])
        for _ in range(config['coverage'])
    ]
    consensus             = build_consensus(noisy_pool, iterations=1)
    frame_id, payload, ok = decode_frame_half(consensus, expected_frame_id)
    return frame_id, payload, ok


def _worker_full(clean_dna, config, original_block, original_pos):
    """ProcessPoolExecutor worker for FULL_MOSAIC."""
    noisy_pool = [
        nanopore_markov_simulator(clean_dna, config['sub'], config['ins'], config['del'])
        for _ in range(config['coverage'])
    ]
    consensus                         = build_consensus(noisy_pool, iterations=1)
    payload, block_idx, pos, ok, muts = decode_frame_full(consensus)
    if ok:
        return block_idx, pos, payload, True, muts
    else:
        return original_block, original_pos, None, False, 0


# ==========================================
# HALF_MOSAIC: OUTER RS RECOVERY  (none)
# ==========================================

def apply_outer_rs_decode(block_payloads, num_blocks, outer_k, outer_ecc, outer_n):
    """
    Column-wise outer Reed-Solomon erasure recovery for FULL_MOSAIC.
    Missing positions in a block are treated as erasures.
    Returns flat bytearray of recovered data bytes.
    """
    recovered_data = bytearray()

    for block_idx in range(num_blocks):
        block             = block_payloads.get(block_idx, {})
        present_positions = set(block.keys())
        missing_positions = sorted(set(range(outer_n)) - present_positions)

        # Build one column per payload byte offset
        received_columns = []
        for j in range(PAYLOAD_BYTES):
            col = bytearray(outer_n)
            for pos in range(outer_n):
                if pos in block:
                    col[pos] = block[pos][j]
            received_columns.append(bytes(col))

        recovered_columns = []
        decode_ok         = True

        for col in received_columns:
            try:
                if missing_positions:
                    decoded = bytes(rs_outer.decode(col, erase_pos=missing_positions)[0])
                else:
                    decoded = bytes(rs_outer.decode(col)[0])
                recovered_columns.append(decoded)
            except reedsolo.ReedSolomonError:
                recovered_columns.append(bytes(outer_k))
                decode_ok = False

        if not decode_ok:
            print(f"\n[!] Outer RS decode failed for block {block_idx} "
                  f"({len(missing_positions)} erasures > {outer_ecc} budget)")

        for pos in range(outer_k):
            payload = bytes(recovered_columns[j][pos] for j in range(PAYLOAD_BYTES))
            recovered_data.extend(payload)

    return recovered_data


# ==========================================
# PIPELINE: HALF_MOSAIC
# ==========================================

def run_decoder_half(fasta_path, coverage, model, dropout_rate, output_path=None):
    """
    Single-layer (inner RS only) decoding pipeline.
    Uncorrectable oligos trigger the Analog Degradation fallback —
    corrupted bytes are written directly into the output stream.
    """
    print(f"{'='*65}")
    print(f"MOSAIC UNIFIED — HALF_MOSAIC (V3 SINGLE-LAYER) PIPELINE")
    print(f"{'='*65}")

    metadata, oligo_list = parse_fasta_half(fasta_path)

    if not oligo_list:
        print("[!] No oligos found in FASTA.")
        return False, 0.0, 0

    original_size   = int(metadata.get('OriginalSize', 0))
    total_frames    = int(metadata.get('TotalDataOligos', len(oligo_list)))
    file_ext        = metadata.get('FileType', 'bin')
    original_sha256 = metadata.get('SHA256', '')

    if output_path is None:
        base        = os.path.splitext(fasta_path)[0]
        output_path = f"{base}_recovered.{file_ext}"

    config = {**ERROR_PROFILES[model], 'coverage': coverage}

    print(f"[+] Architecture        : HALF_MOSAIC")
    print(f"[+] Total Frames        : {total_frames:,}")
    print(f"[+] Original Size       : {original_size:,} bytes (.{file_ext})")
    print(f"[+] Sequencing model    : {model} | Coverage: {coverage}x")
    print(f"[+] Dropout simulation  : {dropout_rate*100:.1f}%")
    print(f"[+] Biological Fallback : ACTIVE (uncorrected mutations pass to file)\n")

    # Physical oligo dropout simulation
    if dropout_rate > 0:
        n_drop     = int(len(oligo_list) * dropout_rate)
        drop_idx   = set(random.sample(range(len(oligo_list)), n_drop))
        oligo_list = [o for i, o in enumerate(oligo_list) if i not in drop_idx]
        print(f"[SIM] Dropped {n_drop} oligos ({dropout_rate*100:.1f}%). "
              f"{len(oligo_list)} remaining.\n")

    # Linear file buffer — frame_id directly indexes payload position
    file_buffer     = bytearray(total_frames * PAYLOAD_BYTES)
    success_count   = 0
    corrupted_count = 0
    cpu_cores       = os.cpu_count() or 4
    start           = time.time()

    with ProcessPoolExecutor(max_workers=cpu_cores) as executor:
        futures = {
            executor.submit(_worker_half, dna, config, frame_id): frame_id
            for frame_id, dna in oligo_list
        }

        done        = 0
        total_tasks = len(oligo_list)

        for future in as_completed(futures):
            frame_id, payload, success = future.result()
            done += 1

            start_idx = frame_id * PAYLOAD_BYTES
            file_buffer[start_idx : start_idx + PAYLOAD_BYTES] = payload

            if success: success_count   += 1
            else:       corrupted_count += 1

            if done % 200 == 0 or done == total_tasks:
                sys.stdout.write(
                    f"\r    Processed: {done}/{total_tasks} | "
                    f"Clean: {success_count} | "
                    f"Corrupted (analog noise): {corrupted_count}"
                )
                sys.stdout.flush()

    elapsed  = time.time() - start
    survival = (success_count / total_frames * 100) if total_frames > 0 else 0.0

    final_data      = bytes(file_buffer[:original_size])
    recovered_sha256 = hashlib.sha256(final_data).hexdigest()
    match           = recovered_sha256 == original_sha256

    print(f"\n\n[DONE] Decoded in {elapsed:.1f}s")
    print(f"    Inner survival rate : {survival:.2f}%")
    print(f"    Corrupted frames    : {corrupted_count} (analog noise injected)")
    print(f"    SHA-256 original    : {original_sha256[:32]}...")
    print(f"    SHA-256 recovered   : {recovered_sha256[:32]}...")
    print(f"    Integrity check     : {'PERFECT MATCH' if match else 'MISMATCH — data loss detected'}")
    print(f"    Output file         : {output_path}")
    print(f"{'='*65}\n")

    with open(output_path, 'wb') as f:
        f.write(final_data)

    return match, survival, corrupted_count


# ==========================================
# PIPELINE: FULL_MOSAIC
# ==========================================

def run_decoder_full(fasta_path, coverage, model, dropout_rate, output_path=None):
    """
    Two-layer (inner RS + outer RS erasure) decoding pipeline.

    Phase 1 — parallel inner RS decode across all oligos.
    Phase 2 — outer RS erasure recovery, reconstructs lost oligos per block.
    """
    print(f"{'='*65}")
    print(f"MOSAIC UNIFIED — FULL_MOSAIC (V4 TWO-LAYER RS) PIPELINE")
    print(f"{'='*65}")

    metadata, oligo_list = parse_fasta_full(fasta_path)

    if not oligo_list:
        print("[!] No oligos found in FASTA.")
        return False, 0.0, 0

    original_size     = int(metadata.get('OriginalSize', 0))
    total_data_oligos = int(metadata.get('TotalDataOligos', 0))
    outer_k           = int(metadata.get('OuterK',   OUTER_RS_K))
    outer_ecc         = int(metadata.get('OuterECC', OUTER_RS_ECC))
    outer_n           = int(metadata.get('OuterN',   OUTER_RS_N))
    file_ext          = metadata.get('FileType', 'bin')
    original_sha256   = metadata.get('SHA256', '')

    # Detect number of RS blocks from metadata or from the oligo count
    detected_blocks = len(oligo_list) // outer_n
    num_blocks      = detected_blocks if detected_blocks > 0 else int(metadata.get('NumBlocks', 1))

    if original_size == 0:
        original_size = num_blocks * outer_k * PAYLOAD_BYTES
        print(f"[+] Warning: OriginalSize missing — calculated {original_size} bytes.")

    # Dynamically update the outer RS codec if metadata differs from defaults
    global rs_outer
    if outer_ecc != OUTER_RS_ECC:
        rs_outer = reedsolo.RSCodec(outer_ecc)

    if output_path is None:
        base        = os.path.splitext(fasta_path)[0]
        output_path = f"{base}_recovered.{file_ext}"

    config = {**ERROR_PROFILES[model], 'coverage': coverage}

    print(f"[+] Architecture        : FULL_MOSAIC")
    print(f"[+] Oligos in pool      : {len(oligo_list):,}")
    print(f"[+] Original file size  : {original_size:,} bytes (.{file_ext})")
    print(f"[+] Outer RS matrix     : K={outer_k}, ECC={outer_ecc}, N={outer_n}")
    print(f"[+] RS blocks           : {num_blocks}")
    print(f"[+] Erasure budget      : {outer_ecc} lost oligos per block")
    print(f"[+] Sequencing model    : {model} | Coverage: {coverage}x")
    print(f"[+] Dropout simulation  : {dropout_rate*100:.1f}%\n")

    # Physical oligo dropout simulation
    if dropout_rate > 0:
        n_drop     = int(len(oligo_list) * dropout_rate)
        drop_idx   = set(random.sample(range(len(oligo_list)), n_drop))
        oligo_list = [o for i, o in enumerate(oligo_list) if i not in drop_idx]
        print(f"[SIM] Dropped {n_drop} oligos ({dropout_rate*100:.1f}%). "
              f"{len(oligo_list)} remaining.\n")

    # ------------------------------------------------------------------
    # PHASE 1 — Parallel inner RS decode
    # ------------------------------------------------------------------
    print("[PHASE 1] Inner RS decode — parallel across all CPU cores...")
    start     = time.time()
    cpu_cores = os.cpu_count() or 4
    total     = len(oligo_list)

    block_payloads  = {}    # {block_idx: {pos: bytes_40}}
    inner_success   = 0
    inner_fail      = 0
    total_mutations = 0

    with ProcessPoolExecutor(max_workers=cpu_cores) as executor:
        futures = {
            executor.submit(_worker_full, dna, config, block_idx, pos): (block_idx, pos)
            for block_idx, pos, dna in oligo_list
        }

        done = 0
        for future in as_completed(futures):
            block_idx, pos, payload, success, mutations = future.result()
            done += 1

            if success and payload is not None:
                if block_idx not in block_payloads:
                    block_payloads[block_idx] = {}
                block_payloads[block_idx][pos] = payload
                inner_success   += 1
                total_mutations += mutations
            else:
                inner_fail += 1

            if done % 200 == 0 or done == total:
                sys.stdout.write(
                    f"\r    {done}/{total} oligos | "
                    f"Decoded: {inner_success} | "
                    f"Erased: {inner_fail} | "
                    f"Mutations corrected: {total_mutations}"
                )
                sys.stdout.flush()

    elapsed_p1     = time.time() - start
    inner_survival = inner_success / total * 100 if total > 0 else 0.0

    print(f"\n\n[PHASE 1] Done in {elapsed_p1:.1f}s")
    print(f"    Inner decode survival   : {inner_survival:.2f}%")
    print(f"    Erased oligos           : {inner_fail}")
    print(f"    Total mutations fixed   : {total_mutations}\n")

    # ------------------------------------------------------------------
    # PHASE 2 — Outer RS erasure recovery
    # ------------------------------------------------------------------
    print("[PHASE 2] Outer RS erasure recovery across pool...")
    start_p2 = time.time()

    for block_idx in range(num_blocks):
        present = len(block_payloads.get(block_idx, {}))
        missing = outer_n - present
        budget  = outer_ecc
        if missing > 0:
            status = "OK" if missing <= budget else "OVER BUDGET"
            print(f"    Block {block_idx:04d}: {missing:3d} erasures / {budget} budget — {status}")

    recovered_bytes = apply_outer_rs_decode(
        block_payloads, num_blocks, outer_k, outer_ecc, outer_n
    )

    recovered_bytes  = bytes(recovered_bytes[:original_size])
    elapsed_p2       = time.time() - start_p2
    recovered_sha256 = hashlib.sha256(recovered_bytes).hexdigest()
    match            = recovered_sha256 == original_sha256

    print(f"\n[PHASE 2] Done in {elapsed_p2:.2f}s")
    print(f"\n{'='*65}")
    print(f"RECOVERY COMPLETE")
    print(f"    Recovered bytes     : {len(recovered_bytes):,} / {original_size:,}")
    print(f"    SHA-256 original    : {original_sha256[:32]}...")
    print(f"    SHA-256 recovered   : {recovered_sha256[:32]}...")
    print(f"    Integrity check     : {'PERFECT MATCH' if match else 'MISMATCH — data loss detected'}")
    print(f"    Inner survival rate : {inner_survival:.2f}%")
    print(f"    Output file         : {output_path}")
    print(f"{'='*65}\n")

    with open(output_path, 'wb') as f:
        f.write(recovered_bytes)

    return match, inner_survival, inner_fail


# ==========================================
# UNIFIED ENTRY POINT
# ==========================================

def run_decoder(fasta_path, coverage, model, dropout_rate, output_path=None):
    """
    Auto-detects architecture from FASTA metadata and dispatches to the
    correct decoding pipeline.  Returns (match, survival_pct, erased_count).
    """
    arch = detect_architecture(fasta_path)
    print(f"[AUTO-DETECT] Architecture detected: {arch}\n")

    if arch == 'HALF_MOSAIC':
        return run_decoder_half(fasta_path, coverage, model, dropout_rate, output_path)
    else:
        return run_decoder_full(fasta_path, coverage, model, dropout_rate, output_path)


# ==========================================
# DROPOUT STRESS SWEEP
# ==========================================

def run_dropout_stress_test(fasta_path, coverage, model):
    """
    Sweeps oligo dropout rates 0%–30% and records recovery outcome.
    Produces a CSV suitable for plotting the paper's dropout figure.
    """
    print(f"\n{'='*65}")
    print(f"MOSAIC UNIFIED — DROPOUT STRESS TEST")
    print(f"Model: {model} | Coverage: {coverage}x")
    print(f"{'='*65}\n")

    dropout_rates = [0.00, 0.02, 0.04, 0.06, 0.08,
                     0.10, 0.12, 0.14, 0.16, 0.20, 0.25, 0.30]

    results  = []
    csv_path = f"dropout_stress_{model}_{coverage}x.csv"

    with open(csv_path, 'w') as csvf:
        csvf.write("dropout_pct,inner_survival_pct,erased_oligos,integrity_ok\n")

    for rate in dropout_rates:
        print(f"\n--- Dropout: {rate*100:.0f}% ---")
        ok, survival, erased = run_decoder(
            fasta_path, coverage, model, rate,
            output_path=f"_stress_tmp_{rate:.2f}.bin"
        )
        results.append((rate, survival, erased, ok))

        with open(csv_path, 'a') as csvf:
            csvf.write(f"{rate*100:.0f},{survival:.2f},{erased},{'1' if ok else '0'}\n")

        try:    os.remove(f"_stress_tmp_{rate:.2f}.bin")
        except: pass

    print(f"\n[>>>] Stress test complete. Results written to {csv_path}")
    print(f"\n{'Dropout':>10} {'Inner surv':>12} {'Erased':>8} {'Integrity':>10}")
    print("-" * 45)
    for rate, surv, erased, ok in results:
        print(f"{rate*100:>9.0f}% {surv:>11.2f}% {erased:>8} {'PASS' if ok else 'FAIL':>10}")

    return results


# ==========================================
# MULTI-FASTA COVERAGE SWEEP
# ==========================================

def run_coverage_sweep_batch(fasta_files, coverage_min, coverage_max,
                              model, dropout_rate=0.0,
                              output_csv="coverage_sweep_results.csv"):
    """
    Iterates over multiple FASTA files across a coverage range.
    Each file is auto-detected for architecture.
    Results are written incrementally to a CSV.
    """
    print(f"\n{'='*70}")
    print("MOSAIC UNIFIED — MULTI-FASTA COVERAGE SWEEP")
    print(f"{'='*70}")

    with open(output_csv, 'w', newline='') as csvfile:
        csv.writer(csvfile).writerow([
            "fasta_file", "architecture", "payload_bytes",
            "coverage", "model", "dropout_rate",
            "sha_verified", "inner_survival_pct", "erased_oligos"
        ])

    results = []

    for fasta_path in fasta_files:
        print(f"\n{'-'*70}")
        print(f"DATASET: {os.path.basename(fasta_path)}")

        try:
            arch = detect_architecture(fasta_path)
        except ValueError as e:
            print(f"[!] Skipping: {e}")
            continue

        meta         = _parse_metadata_lines(fasta_path)
        payload_size = int(meta.get('OriginalSize', 0))

        print(f"Architecture: {arch}")
        print(f"{'-'*70}")

        for coverage in range(coverage_min, coverage_max + 1):
            print(f"\n[TEST] Coverage = {coverage}x")
            try:
                ok, survival, erased = run_decoder(
                    fasta_path=fasta_path,
                    coverage=coverage,
                    model=model,
                    dropout_rate=dropout_rate,
                    output_path=f"_tmp_recovered_{coverage}.bin"
                )

                row = [
                    os.path.basename(fasta_path), arch, payload_size,
                    coverage, model, dropout_rate,
                    1 if ok else 0, round(survival, 4), erased
                ]
                results.append(row)

                with open(output_csv, 'a', newline='') as csvfile:
                    csv.writer(csvfile).writerow(row)

                try:    os.remove(f"_tmp_recovered_{coverage}.bin")
                except: pass

            except Exception as e:
                print(f"[!] FAILURE @ {coverage}x : {e}")

    print(f"\n{'='*70}")
    print(f"COVERAGE SWEEP COMPLETE — Results: {output_csv}")
    print(f"{'='*70}\n")

    return results


# ==========================================
# CLI
# ==========================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "MOSAIC Unified Decoder — auto-detects HALF_MOSAIC or FULL_MOSAIC "
            "architecture from FASTA metadata and routes to the correct pipeline."
        )
    )

    parser.add_argument(
        "fasta", nargs='+',
        help="Input FASTA file(s) or wildcard glob"
    )
    parser.add_argument(
        "--dropout", type=float, default=0.0,
        help="Fraction of oligos to simulate as physically lost (0.0–1.0)"
    )
    parser.add_argument(
        "--coverage", type=int, default=5,
        help="Coverage depth per oligo"
    )
    parser.add_argument(
        "--coverage-range", nargs=2, type=int, metavar=('MIN', 'MAX'),
        help="Coverage sweep range (requires --batch)"
    )
    parser.add_argument(
        "--model", default="R10.4", choices=["R9.4", "R10.4"],
        help="Nanopore error profile"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output filename (auto-derived from metadata if omitted)"
    )
    parser.add_argument(
        "--stress", action="store_true",
        help="Run dropout stress sweep (0%%–30%%) on the first FASTA"
    )
    parser.add_argument(
        "--batch", action="store_true",
        help="Run multi-FASTA coverage sweep"
    )

    args = parser.parse_args()

    # Expand wildcards
    fasta_files = []
    for pattern in args.fasta:
        expanded = glob.glob(pattern)
        fasta_files.extend(expanded if expanded else [pattern])
    fasta_files = list(dict.fromkeys(fasta_files))   # deduplicate, preserve order

    # Validate existence
    missing = [f for f in fasta_files if not os.path.exists(f)]
    if missing:
        print("\n[!] Missing FASTA files:")
        for m in missing: print(f"    {m}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # BATCH MODE
    # ------------------------------------------------------------------
    if args.batch:
        if not args.coverage_range:
            print("[!] --coverage-range MIN MAX is required for --batch mode.")
            sys.exit(1)
        run_coverage_sweep_batch(
            fasta_files=fasta_files,
            coverage_min=args.coverage_range[0],
            coverage_max=args.coverage_range[1],
            model=args.model,
            dropout_rate=args.dropout
        )

    # ------------------------------------------------------------------
    # STRESS TEST MODE
    # ------------------------------------------------------------------
    elif args.stress:
        run_dropout_stress_test(fasta_files[0], args.coverage, args.model)

    # ------------------------------------------------------------------
    # NORMAL SINGLE DECODE
    # ------------------------------------------------------------------
    else:
        run_decoder(
            fasta_path=fasta_files[0],
            coverage=args.coverage,
            model=args.model,
            dropout_rate=args.dropout,
            output_path=args.output
        )
