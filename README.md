# DNA ladder consensus

Reconstruct the 18 bands of the Invitrogen 1 Kb Plus DNA Ladder from Nanopore reads.

## Install

```bash
conda env create -f environment.yml
conda activate ladder-consensus
```

## Run

Use the SUP pass FASTQ files from `ladder_sup_basecalls.tar.gz`:

```bash
python consensus.py --reads /path/to/basecalled_SUP/pass/*.fastq.gz --out consensus
```

The script accepts compressed or uncompressed FASTQ files. Use a new or empty
output directory. To reconstruct selected bands, add, for example, `--bands 100 1000`.
`--threads 4` sets minimap2's thread count. Bands are processed sequentially.

To start from the raw FAST5 files instead, basecall with a separately installed
Guppy 6.5.7 and its `dna_r9.4.1_450bps_sup` model, then run the command above:

```bash
guppy_basecaller -i /path/to/DNA_ladder -s basecalled_SUP \
  -c dna_r9.4.1_450bps_sup.cfg --device cuda:0 --recursive \
  --compress_fastq --num_callers 4 --gpu_runners_per_device 4
```

## Output

- `ladder_sequences.fasta`: per-band consensus, with nominal and actual sizes.
- `split_A.fasta`, `split_B.fasta`: consensuses from disjoint read subsets.
- `ladder_sequences_masked.fasta`: subset disagreements masked with N.
- `band_summary.tsv`: lengths, read counts, endpoint support and subset agreement.
- `uncertainty.bed`: disagreement intervals; zero-based, half-open coordinates.


## Method

Trim ligation adapters and retain reads of at least 50 bp with mean-error Q ≥12.
For each band size, select reads within ±6% or ±10 bp, whichever is larger.
Choose a read seed using up to 40 candidate seeds and 1,200 reads; supporters
must align across at least 88% of both sequences at ≥80% identity. Select up to
400 supporting reads independently from the full pool and each read half.
Combine these selections for the final consensus. Orient reads with
minimap2, align whole reads with abPOA, and set endpoints using modal
adapter-confirmed read starts supported by at least three molecules.
