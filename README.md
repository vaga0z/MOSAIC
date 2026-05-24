# MOSAIC

MOSAIC (Molecular Object Storage Architecture via Independent Coding) is a dual-layer fault-tolerant molecular storage architecture designed for deterministic recovery under insertion/deletion/substitution (IDS)-dominated nanopore sequencing channels.

The system combines:

* Inner Reed–Solomon correction
* Outer Reed–Solomon erasure recovery
* LFSR-based data whitening
* Constraint-aware synthesis validation
* Coverage threshold characterization
* Parallelized nanopore consensus simulation

to achieve SHA-256-verified file recovery under stochastic nanopore degradation.

---

# Architecture Overview

MOSAIC supports two operational recovery modes:

## FULL MOSAIC

Deterministic archival recovery architecture with:

* outer erasure isolation,
* block-level recovery guarantees,
* probabilistic convergence modeling,
* and SHA-256-verified reconstruction.

Best suited for:

* archival storage,
* scientific datasets,
* cryptographic payloads,
* and integrity-critical recovery scenarios.

---

## HALF MOSAIC

High-density continuous-stream recovery mode without outer erasure isolation.

This configuration permits graceful analog degradation and perceptual recovery under severe nanopore corruption.

Best suited for:

* perceptually tolerant payloads,
* multimedia recovery experiments,
* and degradation analysis.

---

# Features

* Dual-layer fault-tolerant recovery
* Reed–Solomon inner and outer coding
* Nanopore IDS channel simulation
* Coverage sweep benchmarking
* Entropy invariance analysis
* SHA-256 verification
* Multi-FASTA batch processing
* Parallelized encoding and decoding
* Automatic architecture detection from FASTA metadata

---

# FASTA Metadata Standard

MOSAIC embeds architecture and recovery metadata directly into FASTA manifests.

Example:

```text
; Architecture: FULL_MOSAIC
; OriginalSize: 1000000
; FileType: png
; SHA256: ...
; TotalDataOligos: ...
; OuterECC: ...
```

The decoder automatically detects:

* FULL_MOSAIC
* HALF_MOSAIC

and routes recovery accordingly.

---

# Installation

```bash
pip install -r requirements.txt
```

---

# Encoding

## FULL MOSAIC

```bash
python MOSAIC_ENCODER.py input.bin output.fasta full
```

## HALF MOSAIC

```bash
python MOSAIC_ENCODER.py input.bin output.fasta half
```

---

# Decoding

The decoder automatically detects architecture state from FASTA metadata.

```bash
python MOSAIC_DECODER.py file.fasta --coverage 10 --model R9.4
```

---

# Coverage Sweep Benchmarking

```bash
python MOSAIC_DECODER.py datasets/*.fasta \
    --coverage-range 7 13 \
    --model R9.4 \
    --batch
```

Outputs:

* SHA-256 verification
* Inner survival rate
* Coverage convergence behavior
* Oligo erasure statistics
* CSV benchmark logs

---


# Research Context

This repository accompanies the manuscript:

“MOSAIC: A Dual-Layer Fault-Tolerant Molecular Storage Architecture with Empirical Coverage Threshold Characterization for Nanopore Sequencing Channels”

submitted to:

IEEE Transactions on Molecular, Biological, and Multi-Scale Communications (TMBMC)

---

# Disclaimer

This repository reflects the reference implementation associated with the submitted manuscript and remains under active refinement.
