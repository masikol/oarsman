# oarsman

Oarsman (**O**verlapping **A**mplicons, **R**eference) is a pipeline for obtaining consensus genome sequences from overlapping amplicon data using a reference genome sequence.

Current version is 0.2.b (2022-01-25 edition).

## Description

### Input data:

1. Reads (fastq): single-end or paired-end.
2. A CSV file of primers used to create amplicons for sequencing.
3. A reference sequence (fasta).

### Output data:

1. The consensus sequence (fasta) with variants applied according to the reads (output subdirectory `consensus`).
2. The consensus sequence with low-coverage regions annotated (files `consensus/*annotated_consensus.gbk`).
3. Various temporary files: reads cleaned by kromsatel, `.bam` mapping file, `.bcf` variant file etc.

### Dependencies

Oarsman is written in Python 3 and won't work on Python 2.

Oarsman is a pipeline, so it has an awfully massive bunch of dependencies.

1. kromsatel (version 2.0.a or later). Might be in PATH or specified with option `--kromsatel`.
2. consensus-highlighter (version 2.1.a or later). Might be in PATH or specified with option `--highlighter`.
3. samtools (version 1.12 or later). Might be in PATH or specified with option `--samtools`.
4. bcftools. Might be in PATH or specified with option `--bcftools`.
5. BLAST+. Must be in PATH.
6. bwa. Might be in PATH or specified with option `--bwa`.

## Usage

### Single-end reads, single sample

```
python3 oarsman.py \
    -1 20_S30_L001.fastq.gz \                           # unpaired reads
    -p nCov-2019-alt_primers.csv \                      # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \                # Reference sequence
    -t 4 \                                              # CPU threads
    -o my_outdir                                        # Output directory
```

### Single-end reads, multiple samples

```
python3 oarsman.py \
    -1 20_S30_L001.fastq.gz \                           # unpaired reads, sample 1
    -1 1_S7_L001.fastq.gz \                             # unpaired reads, sample 2
    -p nCov-2019-alt_primers.csv \                      # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \                # Reference sequence
    -t 4 \                                              # CPU threads
    -o my_outdir                                        # Output directory
```

### Paired-end reads, single sample

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \                    # forward (R1) reads
    -2 20_S30_L001_R2_001.fastq.gz \                    # reverse (R2) reads
    -p nCov-2019-alt_primers.csv \                      # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \                # Reference sequence
    -t 4 \                                              # CPU threads
    -o my_outdir                                        # Output directory
```

### Paired-end reads, multiple samples

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \                    # forward (R1) reads, sample 1
    -2 20_S30_L001_R2_001.fastq.gz \                    # reverse (R2) reads, sample 1
    -1 1_S7_L001_R1_001.fastq.gz \                      # forward (R1) reads, sample 2
    -2 1_S7_L001_R2_001.fastq.gz \                      # reverse (R2) reads, sample 2
    -p nCov-2019-alt_primers.csv \                      # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \                # Reference sequence
    -t 4 \                                              # CPU threads
    -o my_outdir                                        # Output directory
```

### If reads from a single sample is split into multiple files

Let the forward reads be stored in files `1_S7_L001_R1_001_part_one.fastq.gz` and `1_S7_L001_R1_001_part_two.fastq.gz`. And let the reverse reads be stored in files `1_S7_L001_R2_001_part_one.fastq.gz` and `1_S7_L001_R2_001_part_two.fastq.gz`. So, the command will be the following:

```
python3 oarsman.py \
    -1 1_S7_L001_R1_001_part_one.fastq.gz 1_S7_L001_R1_001_part_two.fastq.gz \
    -2 1_S7_L001_R2_001_part_one.fastq.gz 1_S7_L001_R2_001_part_two.fastq.gz \
    -p nCov-2019-alt_primers.csv \
    -r Wuhan-Hu-1-compele-genome.fasta \
    -t 4 \
    -o my_outdir
```

### With some dependencies

Here, some dependencies (kromsatel, consensus-highlighter, and bwa) are specified in command line instead of using PATH environment variable.

```
python3 oarsman.py \
    -1 20_S30_L001.fastq.gz \
    -p nCov-2019-alt_primers.csv \
    -r Wuhan-Hu-1-compele-genome.fasta \
    --kromsatel /some/dir/kromsatel/kromsatel.py \                            # Kromsatel executable file
    --highlighter /some/dir/consensus-highlighter/consensus-highlighter.py \  # consensus-highlighter executable script
    --bwa /some/dir/bwa \
    -t 4 \
    -o my_outdir
```
