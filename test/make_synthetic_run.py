#!/usr/bin/env python3
"""Simulate a Q16Seq Nanopore run: 32 samples x 3 lambda standards = 96 inner-barcoded wells.

Writes a MinKNOW-like folder (fastq_pass/barcode01/*.fastq.gz) that RT-RubyRed / Q16Seq can
watch, plus the ground truth:

  <out>/fastq_pass/barcode01/SYNTH0001_pass_barcode01_q16seqsim_N.fastq.gz
  <out>/truth_copies_per_reaction.tsv   true 16S copies per reaction, taxon x sample
  <out>/truth_reads_per_well.tsv        simulated reads per taxon (+ lambda) per BRK well
  <out>/lambda_map.tsv                  BRK well -> lambda standard (L1/L2/L3) and sample (group)
  <out>/samples.tsv                     sample design (total copies, taxa present)

Design
  * Sample i (S01..S32) has total bacterial 16S copies per reaction log-spaced 1e3..1e7 and a
    random community of 3-10 of the reference taxa (Dirichlet weights).
  * Each sample is split into three wells spiked with lambda1 = L1 copies, lambda2 = L1/10 and
    lambda3 = L1/100 (default L1 = 7.7e6), placed in BRK(3i-2), BRK(3i-1), BRK(3i).
  * Reads per well ~ lognormal around --reads-per-well; reads are drawn multinomially in
    proportion to copies (no PCR bias), so Q16Seq should recover the truth up to sampling noise.
  * Read = spacer + inner barcode + linker + primer tail + amplicon + reverse tail, with ~2 %
    Nanopore-like errors (substitutions/insertions/deletions); half the reads reverse-complemented.
    Bacterial amplicons span 27F..1391R; the lambda-IS construct carries 27F/1492R tails.
  * A few percent short fragments and reads without barcode are added, as in real runs.

Usage
  python make_synthetic_run.py -o synthetic_run [--reads-per-well 10000] [--lambda1 7.7e6] [--seed 16]
"""
import argparse
import gzip
import os
import re
import uuid

import numpy as np

HERE = os.path.dirname(os.path.realpath(__file__))
REFS = os.path.join(HERE, "synthetic_refs")
COMP = str.maketrans("ACGTN", "TGCAN")
F27 = "AGAGTTTGATCATGGCTCAG"


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                out[name] = []
            elif name:
                out[name].append(line.upper())
    iupac = str.maketrans("RYSWKMBDHVN", "ACGAGACAAAA")  # ambiguity codes -> one valid base
    return {k: "".join(v).translate(iupac) for k, v in out.items()}


def barcodes(path):
    """torchlex barcode FASTA: entries of 4 lines (spacer, barcode, linker, primer tail)."""
    out, name, lines = {}, None, []
    with open(path) as fh:
        for line in list(fh) + [">end"]:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out[name] = lines
                name, lines = line[1:], []
            elif line:
                lines.append(re.sub("[^ACGT]", "", line.upper()))  # the file has a few stray quotes
    return out


def amplicon(seq):
    """Bacterial 16S trimmed to the 27F..1391R amplicon (27F primer sequence prepended)."""
    m = re.search("AGAGTTTGAT[CT][AC]TGGCTCAG", seq)
    start = m.end() if m else 0
    e = re.search("TG[CT]ACACACCGCCCGT", seq)
    end = e.end() if e else len(seq)
    return F27 + seq[start:end]


def taxon_name(header):
    w = header.split()
    return f"{w[1]} {w[2]}"


def mutate(rng, seq, err):
    """Nanopore-like errors: 60 % substitutions, 20 % insertions, 20 % deletions; returns seq, qual."""
    acgt = np.frombuffer(b"ACGT", dtype=np.uint8)
    s = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    n = len(s)
    q = np.clip(rng.normal(22, 4, n), 10, 40).astype(np.uint8)
    ev = rng.random(n) < err
    kind = rng.random(n)
    sub = ev & (kind < 0.6)
    ins = ev & (kind >= 0.6) & (kind < 0.8)
    dele = ev & (kind >= 0.8)
    s[sub] = acgt[rng.integers(0, 4, sub.sum())]
    q[sub] = rng.integers(4, 12, sub.sum())
    pos = np.where(ins)[0] + 1
    s = np.insert(s, pos, acgt[rng.integers(0, 4, len(pos))])
    q = np.insert(q, pos, rng.integers(4, 12, len(pos)).astype(np.uint8))
    keep = np.insert(~dele, pos, True)
    s, q = s[keep], q[keep]
    return s.tobytes().decode(), (q + 33).tobytes().decode()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--samples", type=int, default=32)
    ap.add_argument("--reads-per-well", type=int, default=10000)
    ap.add_argument("--lambda1", type=float, default=7.7e6, help="lambda1 copies per reaction")
    ap.add_argument("--min-copies", type=float, default=1e3)
    ap.add_argument("--max-copies", type=float, default=1e7)
    ap.add_argument("--error-rate", type=float, default=0.02)
    ap.add_argument("--reads-per-file", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=16)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    refs = read_fasta(os.path.join(REFS, "16S_refs.fasta"))
    taxa = [taxon_name(h) for h in refs]
    amps = [amplicon(s) for s in refs.values()]
    lam = list(read_fasta(os.path.join(REFS, "lambda_IS.fasta")).values())[0]
    bcs = barcodes(os.path.join(REFS, "torchlex_barcodes.fasta"))
    n_wells = 3 * a.samples
    if n_wells > 192:
        raise SystemExit("at most 64 samples (192 inner barcodes)")

    # ---- design
    totals = np.logspace(np.log10(a.min_copies), np.log10(a.max_copies), a.samples)
    comp = np.zeros((a.samples, len(taxa)))
    for i in range(a.samples):
        k = rng.integers(3, min(10, len(taxa)) + 1)
        idx = rng.choice(len(taxa), k, replace=False)
        comp[i, idx] = rng.dirichlet(np.full(k, 0.8))
    copies = comp * totals[:, None]                                  # sample x taxon
    lam_levels = [a.lambda1, a.lambda1 / 10, a.lambda1 / 100]

    os.makedirs(a.out, exist_ok=True)
    snames = [f"S{i + 1:02d}" for i in range(a.samples)]
    with open(os.path.join(a.out, "truth_copies_per_reaction.tsv"), "w") as fh:
        fh.write("taxon\t" + "\t".join(snames) + "\n")
        for t, name in enumerate(taxa):
            fh.write(name + "\t" + "\t".join(f"{copies[i, t]:.6g}" for i in range(a.samples)) + "\n")
        fh.write("Total bacterial 16S copies\t" + "\t".join(f"{x:.6g}" for x in totals) + "\n")
    with open(os.path.join(a.out, "samples.tsv"), "w") as fh:
        fh.write("sample\ttotal_16S_copies_per_reaction\tn_taxa\ttaxa_present\twells\n")
        for i, s in enumerate(snames):
            present = [taxa[t] for t in np.argsort(-copies[i]) if copies[i, t] > 0]
            wells = ",".join(f"BRK{3 * i + j + 1:02d}" for j in range(3))
            fh.write(f"{s}\t{totals[i]:.6g}\t{len(present)}\t{'; '.join(present)}\t{wells}\n")
    with open(os.path.join(a.out, "lambda_map.tsv"), "w") as fh:
        fh.write("sample\tlambda\tgroup\n")
        for i, s in enumerate(snames):
            for j in range(3):
                fh.write(f"BRK{3 * i + j + 1:02d}\tL{j + 1}\t{s}\n")

    # ---- reads per well
    reads, truth_rows = [], []
    for i in range(a.samples):
        for j in range(3):
            well = 3 * i + j + 1
            brk = f"BRK{well:02d}"
            w = np.append(copies[i], lam_levels[j])
            n = max(200, int(rng.lognormal(np.log(a.reads_per_well), 0.3)))
            counts = rng.multinomial(n, w / w.sum())
            truth_rows.append((brk, snames[i], f"L{j + 1}", lam_levels[j], counts))
            s1, s2 = bcs[f"{brk}_1st"], bcs[f"{brk}_2nd"]
            left = s1[0] + s1[1] + s1[2] + s1[3]
            right = rc(s2[3]) + rc(s2[2]) + rc(s2[1]) + rc(s2[0])
            for t, c in enumerate(counts):
                src = amps[t] if t < len(taxa) else lam
                reads += [(left + src + right, brk)] * int(c)
    # short fragments and reads without inner barcode (~2 % each)
    junk = []
    for _ in range(int(0.02 * len(reads))):
        s = amps[rng.integers(len(amps))]
        st = rng.integers(0, len(s) - 700)
        junk.append((s[st:st + rng.integers(250, 700)], "junk"))
    for _ in range(int(0.02 * len(reads))):
        junk.append((amps[rng.integers(len(amps))], "no_barcode"))
    reads += junk
    order = rng.permutation(len(reads))

    with open(os.path.join(a.out, "truth_reads_per_well.tsv"), "w") as fh:
        fh.write("well\tsample\tlambda\tlambda_copies\t" + "\t".join(taxa) + "\tlambda_IS\ttotal\n")
        for brk, s, lv, lc, c in truth_rows:
            fh.write(f"{brk}\t{s}\t{lv}\t{lc:.4g}\t" + "\t".join(map(str, c)) + f"\t{c.sum()}\n")

    # ---- MinKNOW-like FASTQ chunks
    d = os.path.join(a.out, "fastq_pass", "barcode01")
    os.makedirs(d, exist_ok=True)
    fh, nfile = None, 0
    for k, r in enumerate(order):
        if k % a.reads_per_file == 0:
            if fh:
                fh.close()
            fh = gzip.open(os.path.join(d, f"SYNTH0001_pass_barcode01_q16seqsim_{nfile}.fastq.gz"),
                           "wt", compresslevel=4)
            nfile += 1
        seq, _ = reads[r]
        if rng.random() < 0.5:
            seq = rc(seq)
        seq, qual = mutate(rng, seq, a.error_rate)
        fh.write(f"@{uuid.UUID(int=int(rng.integers(2**63)))} runid=q16seqsim sample_id=Q16Seq_synthetic "
                 f"barcode=barcode01\n{seq}\n+\n{qual}\n")
    fh.close()
    print(f"{len(reads):,} reads ({len(junk):,} junk) in {nfile} files -> {d}")
    print(f"{a.samples} samples x 3 lambda standards = {n_wells} wells (BRK01-BRK{n_wells:02d}); "
          f"lambda1 = {a.lambda1:.3g} copies; truth tables in {a.out}")


if __name__ == "__main__":
    main()
