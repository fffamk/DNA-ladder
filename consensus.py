#!/usr/bin/env python3
"""Reconstruct DNA ladder bands from basecalled Nanopore reads."""

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
import hashlib
import math
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory

import edlib
import pyabpoa

BANDS = [100, 200, 300, 400, 500, 650, 850, 1000, 1500, 2000,
         3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000]
ADAPTER = "ACTTCGTTCAGTTACGTATTGCT"


@dataclass
class Read:
    name: str
    seq: str
    score: float
    partition: str
    adapter: bool


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def fastq(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while header := handle.readline():
            seq = handle.readline().strip().upper()
            plus = handle.readline()
            qual = handle.readline().rstrip("\r\n")
            if (not header.startswith("@") or not header[1:].strip()
                    or not plus.startswith("+") or not seq
                    or len(seq) != len(qual) or set(seq) - set("ACGTN")
                    or any(not 33 <= ord(c) <= 126 for c in qual)):
                raise ValueError(f"Malformed FASTQ record in {path}: {header.strip()}")
            yield header[1:].split()[0], seq, qual


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def trim(seq, qual):
    left, right = 0, len(seq)
    match = edlib.align(ADAPTER, seq[:75], mode="HW", task="locations", k=5)
    if match["editDistance"] >= 0:
        start, end = min(match["locations"], key=lambda loc: loc[0])
        if start <= 35:
            left = end + 1
    adapter = revcomp(ADAPTER)
    for length in range(len(adapter), 7, -1):
        allowance = 0 if length < 12 else int(length * .16)
        candidates = []
        for extra in range(3):
            for observed in range(max(8, length - allowance), length + allowance + 1):
                start = len(seq) - extra - observed
                if start <= left:
                    continue
                target = seq[start:len(seq) - extra if extra else len(seq)]
                match = edlib.align(adapter[:length], target, mode="NW",
                                    task="distance", k=allowance)
                if match["editDistance"] >= 0:
                    candidates.append((length, -extra, -abs(length - observed), start))
        if candidates:
            right = max(candidates)[3]
            break
    return seq[left:right], qual[left:right], left > 0


def load_reads(paths):
    reads, seen = [], set()
    for path in sorted(paths):
        for name, seq, qual in fastq(path):
            if name in seen:
                raise ValueError(f"Duplicate read ID: {name}")
            seen.add(name)
            seq, qual, adapter = trim(seq, qual)
            if len(seq) < 50:
                continue
            score = -10 * math.log10(sum(10 ** (-(ord(c) - 33) / 10)
                                         for c in qual) / len(qual))
            if score >= 12:
                part = "AB"[int(hashlib.sha256(name.encode()).hexdigest()[:8], 16) % 2]
                reads.append(Read(name, seq, score, part, adapter))
        print(f"Read {path.name}", flush=True)
    return reads


def align(reference, queries, path, band, threads):
    short = band <= 200
    command = ["minimap2", "-c", "-x", "map-ont", "-t", str(threads),
               "-k", "7" if short else "13", "-w", "3" if short else "5",
               "-m", "15" if short else "40", "-s", "15" if short else "40",
               "-n", "2", "--secondary=yes", "-N", "100", "-p", "0.5",
               str(reference), str(queries)]
    with open(path, "w") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr)


def paf(path):
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            yield (fields[0], int(fields[1]), int(fields[2]), int(fields[3]),
                   fields[4], fields[5], int(fields[6]), int(fields[7]),
                   int(fields[8]), int(fields[9]), int(fields[10]))


def full_length(row):
    _, qlen, qstart, qend, _, _, tlen, tstart, tend, matches, block = row
    return ((qend - qstart) / qlen >= .88 and (tend - tstart) / tlen >= .88
            and matches / max(1, block) >= .80)


def select_reads(pool, band, work, threads):
    if len(pool) < 5:
        raise ValueError(f"Band {band}: fewer than five qualifying reads in a partition")
    pool = sorted(pool, key=lambda r: (r.score, -abs(len(r.seq) - band)), reverse=True)
    probes = pool[:24] + sorted(pool[24:], key=lambda r: hashlib.sha256(
        r.name.encode()).hexdigest())[:16]
    ref, query, mapping = work / "seeds.fa", work / "reads.fa", work / "reads.paf"
    write_fasta(ref, ((r.name, r.seq) for r in probes))
    write_fasta(query, ((r.name, r.seq) for r in pool[:1200]))
    align(ref, query, mapping, band, threads)
    supporters = defaultdict(set)
    for row in paf(mapping):
        if full_length(row):
            supporters[row[5]].add(row[0])
    seed = max(probes, key=lambda r: (len(supporters[r.name]), r.score,
                                     -abs(len(r.seq) - band)))
    write_fasta(ref, [(f"band_{band}", seed.seq)])
    write_fasta(query, ((r.name, r.seq) for r in pool))
    align(ref, query, mapping, band, threads)
    support = {row[0] for row in paf(mapping) if full_length(row)}
    selected = [r for r in pool if r.name in support][:400]
    if len(selected) < 5:
        raise ValueError(f"Band {band}: fewer than five full-length supporting reads")
    return seed.seq, selected


def orient(seq, reference):
    reverse = revcomp(seq)
    if edlib.align(reference, reverse)["editDistance"] < edlib.align(
            reference, seq)["editDistance"]:
        return reverse
    return seq


def consensus(reads, seed, band, work, threads):
    ref, query, mapping = work / "seed.fa", work / "selected.fa", work / "oriented.paf"
    write_fasta(ref, [("seed", seed)])
    write_fasta(query, ((r.name, r.seq) for r in reads))
    align(ref, query, mapping, band, threads)
    best = {}
    for row in paf(mapping):
        if row[0] not in best or row[10] > best[row[0]][1]:
            best[row[0]] = (row[4], row[10])
    reads = [r for r in reads if r.name in best]
    if len(reads) < 5:
        raise ValueError(f"Band {band}: fewer than five oriented reads")
    sequences = [r.seq if best[r.name][0] == "+" else revcomp(r.seq) for r in reads]
    alignment = pyabpoa.msa_aligner().msa(sequences, out_cons=True, out_msa=True)
    rows, draft = alignment.msa_seq[:-1], alignment.msa_seq[-1]
    left, right = [], []
    for read, row in zip(reads, rows):
        if read.adapter:
            if best[read.name][0] == "+":
                left.append(len(row) - len(row.lstrip("-")))
            else:
                right.append(len(row.rstrip("-")))
    bounds, support = [], []
    for positions, fallback in [(left, 0), (right, len(draft))]:
        mode, count = Counter(positions).most_common(1)[0] if positions else (fallback, 0)
        bounds.append(mode if count >= 3 else fallback)
        support.append(count)
    seq = draft[bounds[0]:bounds[1]].replace("-", "")
    if not seq:
        raise ValueError(f"Band {band}: empty consensus")
    return seq, len(reads), support


def uncertain_intervals(reference, comparator):
    cigar = edlib.align(reference, comparator, mode="NW", task="path")["cigar"]
    pos = 0
    for count, op in re.findall(r"(\d+)([=XID])", cigar):
        count = int(count)
        if op != "=":
            # I consumes reference; D consumes comparator. Include indel anchors.
            start = max(0, pos - 1) if op in "ID" else pos
            end = min(len(reference), pos + count + 1 if op in "IX" else pos + 1)
            yield start, end
        if op in "=XI":
            pos += count


def build_band(band, reads, work, threads):
    pool = [r for r in reads if abs(len(r.seq) - band) <= max(10, band * .06)]
    seeds, selected = {}, {}
    for part in ["all", "A", "B"]:
        subset = [r for r in pool if part == "all" or r.partition == part]
        seeds[part], selected[part] = select_reads(subset, band, work, threads)
    union = {r.name: r for part in ["all", "A", "B"] for r in selected[part]}
    # Match the four-decimal quality scores used in the reconstruction metadata.
    selected["all"] = sorted(union.values(), key=lambda r: (-float(f"{r.score:.4f}"), r.name))
    sequences = {}
    summary = {"nominal_band_size_bp": band, "candidate_reads": len(pool)}
    for part in ["all", "A", "B"]:
        seed = orient(seeds[part], seeds["all"])
        seq, count, ends = consensus(selected[part], seed, band, work, threads)
        sequences[part] = seq if part == "all" else orient(seq, sequences["all"])
        summary[f"{part}_reads_used"] = count
        summary[f"{part}_left_endpoint_support"] = ends[0]
        summary[f"{part}_right_endpoint_support"] = ends[1]
    summary["actual_size_bp"] = len(sequences["all"])
    mask = list(sequences["all"])
    intervals = []
    for part in ["A", "B"]:
        summary[f"all_vs_{part}_edits"] = edlib.align(
            sequences["all"], sequences[part])["editDistance"]
        for start, end in uncertain_intervals(sequences["all"], sequences[part]):
            intervals.append((start, end, f"all_vs_{part}"))
            mask[start:end] = "N" * (end - start)
    summary["masked_bases"] = mask.count("N")
    summary["status"] = ("concordant" if sequences["all"] == sequences["A"] == sequences["B"]
                         else "local_uncertainty")
    sequences["masked"] = "".join(mask)
    return sequences, summary, intervals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", nargs="+", type=Path, required=True, help="SUP pass FASTQ files")
    parser.add_argument("--out", type=Path, required=True, help="New or empty output directory")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--bands", nargs="+", type=int, default=BANDS, help="Nominal band sizes in bp")
    args = parser.parse_args()
    if args.threads < 1 or any(b <= 0 for b in args.bands) or len(set(args.bands)) != len(args.bands):
        parser.error("Threads and band sizes must be positive; band sizes must be unique")
    if not shutil.which("minimap2"):
        parser.error("minimap2 is not on PATH")
    if any(not p.is_file() for p in args.reads):
        parser.error("One or more FASTQ files do not exist")
    if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        parser.error("Output directory must be new or empty")
    args.out.mkdir(parents=True, exist_ok=True)
    reads = load_reads(args.reads)
    print(f"Retained {len(reads)} reads with length >=50 bp and Q>=12", flush=True)
    results = []
    with TemporaryDirectory(prefix="work_", dir=args.out) as temp:
        for band in sorted(args.bands):
            sequences, summary, intervals = build_band(band, reads, Path(temp), args.threads)
            results.append((band, sequences, summary, intervals))
            print(f"{band} bp: {summary['actual_size_bp']} bp, {summary['status']}", flush=True)
    for part, filename in [("all", "ladder_sequences.fasta"), ("A", "split_A.fasta"),
                           ("B", "split_B.fasta"), ("masked", "ladder_sequences_masked.fasta")]:
        write_fasta(args.out / filename, (
            (f"ladder_{band}bp nominal_band_size_bp={band} actual_size_bp={len(seqs[part])}", seqs[part])
            for band, seqs, _, _ in results))
    with open(args.out / "band_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0][2]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary for _, _, summary, _ in results)
    with open(args.out / "uncertainty.bed", "w") as handle:
        for band, _, _, intervals in results:
            for start, end, label in sorted(intervals):
                handle.write(f"ladder_{band}bp\t{start}\t{end}\t{label}\n")


if __name__ == "__main__":
    main()
