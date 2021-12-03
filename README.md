# oarsman

Oarsman (**O**verlapping **A**mplicons, **R**eference) is a pipeline for obtaining consensus genome sequences from overlapping amplicon data using a reference genome sequence.

## Description

### Input data:

1. Reads (fastq): single-end or paired-end.
2. A CSV file of primers used to create amplicons for sequencing.
3. A reference sequence (fasta).

### Output data:

1. The consensus sequence (fasta) with variants applied according to the reads (output subdirectory `consensus`).
2. Various temporary files: reads cleaned by kromsatel, `.bam` mapping file, `.bcf` variant file etc.


## Usage

### Single-end reads, single sample

```
python3 oarsman.py \
    -1 20_S30_L001.fastq.gz \                          # unpaired reads
    -p nCov-2019-alt_primers.csv \                     # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \               # Reference sequence
    --kromsatel-dir /mnt/1.5_drive_0/repos/kromsatel \ # Kromsatel source directory
    -t 4 \                                             # CPU threads
    -o my_outdir                                       # Output directory
```

### Single-end reads, multiple samples

```
python3 oarsman.py \
    -1 20_S30_L001.fastq.gz \                          # unpaired reads, sample 1
    -1 1_S7_L001.fastq.gz \                            # unpaired reads, sample 2
    -p nCov-2019-alt_primers.csv \                     # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \               # Reference sequence
    --kromsatel-dir /mnt/1.5_drive_0/repos/kromsatel \ # Kromsatel source directory
    -t 4 \                                             # CPU threads
    -o my_outdir                                       # Output directory
```

### Paired-end reads, single sample

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \                   # forward (R1) reads
    -2 20_S30_L001_R2_001.fastq.gz \                   # reverse (R2) reads
    -p nCov-2019-alt_primers.csv \                     # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \               # Reference sequence
    --kromsatel-dir /mnt/1.5_drive_0/repos/kromsatel \ # Kromsatel source directory
    -t 4 \                                             # CPU threads
    -o my_outdir                                       # Output directory
```

### Paired-end reads, multiple samples

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \                   # forward (R1) reads, sample 1
    -2 20_S30_L001_R2_001.fastq.gz \                   # reverse (R2) reads, sample 1
    -1 1_S7_L001_R1_001.fastq.gz \                     # forward (R1) reads, sample 2
    -2 1_S7_L001_R2_001.fastq.gz \                     # reverse (R2) reads, sample 2
    -p nCov-2019-alt_primers.csv \                     # ARTIC primers
    -r Wuhan-Hu-1-compele-genome.fasta \               # Reference sequence
    --kromsatel-dir /mnt/1.5_drive_0/repos/kromsatel \ # Kromsatel source directory
    -t 4 \                                             # CPU threads
    -o my_outdir                                       # Output directory
```

### If reads from a single sample is split into multiple files

Let the forward reads be stored in files `1_S7_L001_R1_001_part_one.fastq.gz` and `1_S7_L001_R1_001_part_two.fastq.gz`. And let the reverse reads be stored in files `1_S7_L001_R2_001_part_one.fastq.gz` and `1_S7_L001_R2_001_part_two.fastq.gz`. So, the command will be the following:

```
python3 oarsman.py \
    -1 1_S7_L001_R1_001_part_one.fastq.gz 1_S7_L001_R1_001_part_two.fastq.gz \
    -2 1_S7_L001_R2_001_part_one.fastq.gz 1_S7_L001_R2_001_part_two.fastq.gz \
    -p nCov-2019-alt_primers.csv \
    -r Wuhan-Hu-1-compele-genome.fasta \
    --kromsatel-dir /mnt/1.5_drive_0/repos/kromsatel \
    -t 4 \
    -o my_outdir
```

## Dependencies

1. kromsatel
2. Python packages `pandas` and `numpy`.
3. samtools
4. bcftools
5. BLAST+
6. bwa
7. seqkit