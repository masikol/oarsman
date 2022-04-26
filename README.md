# Oarsman

Oarsman (**O**verlapping **A**mplicons, **R**eference) is a pipeline for obtaining consensus genome sequences from overlapping amplicon data using a reference genome sequence.

Current version is 0.2.f (2022-04-26 edition).

## Description

### Input data:

1. Reads (FASTQ). Oarsman currently can process only Illumina paired-end data. I plan to add the support of Nanopore data later.
2. A CSV file of primers used to create amplicons for sequencing. The file must be comma-separated, with two columns: primer name and primer sequence, without header. Example: [https://github.com/masikol/kromsatel/blob/main/primers/nCov-2019_primers.csv](https://github.com/masikol/kromsatel/blob/main/primers/nCov-2019_primers.csv).
3. A reference sequence in fasta format. Usually, the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) genome is used as a reference.

### Output data:

1. The consensus sequence (fasta) with variants applied according to the reads (output subdirectory `consensus`).
2. The consensus sequence with low-coverage regions annotated (files `consensus/*annotated_consensus.gbk`).
3. Various temporary files: reads cleaned by Kromsatel, `.bam` mapping file, `.bcf` variant file etc. They might be also useful if some loci are suspicious after Oarsman run.

### Dependencies

Oarsman is written in Python 3 and won't work on Python 2.

Oarsman is a pipeline, so it has an awfully massive bunch of dependencies.

1. Kromsatel (version 2.0.a or later). Might be in PATH or specified with option `--kromsatel`. Kromsatel might be obtained [here](https://github.com/masikol/kromsatel).
2. consensus-highlighter (version 2.1.a or later). Might be in PATH or specified with option `--highlighter`. Consensus-highlighter might be obtained [here](https://github.com/masikol/consensus-highlighter).
3. SAMTools (version 1.12 or later). Might be in PATH or specified with option `--samtools`. SAMTools might be obtained [here](https://www.htslib.org/download/).
4. BCFTools. Might be in PATH or specified with option `--bcftools`. BCFTools might be obtained [here](https://www.htslib.org/download/).
5. BLAST+. Must be in PATH. Oarsman tested on BLAST+ version 2.9.0. BLAST+ might be obtained [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
6. BWA. Might be in PATH or specified with option `--bwa`. Oarsman tested on BWA version 0.7.17-r1188. BLAST+ might be obtained [here](https://github.com/lh3/bwa).
7. LoFreq. Optional. If you want to call variants with LoFreq, specify path to its binary with option `--lofreq`. LoFreq might be obtained [here](https://sourceforge.net/projects/lofreq/files/). Oarsman tested on LoFreq version 2.1.2. See example on how to run LoFreq in ["Running Oarsman with LoFreq" section](#running-oarsman-with-lofreq) below.

### The pipeline steps

1. Pre-process reads (remove adapter and primer sequeences, split chimeric reads) with Kromsatel.
2. Map the pre-processed reads on the reference sequence with BWA.
3. Call variants with BCFTools or with LoFreq.
4. Make a consensus sequence, i.e. the sequence of the genome being actually sequenced, with BCFTools.
5. Annotate the consensus sequence with consensus-highlighter.

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

Here, some dependencies (Kromsatel, consensus-highlighter, and bwa) are specified in command line instead of using PATH environment variable.

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

### Running Oarsman with LoFreq


If you pass path to `lofreq` binary with `--lofreq` option, LoFreq will be used to call variants instead of default BCFTools. LoFreq seems to be more accurate ([Deng et al.](https://academic.oup.com/bib/article/22/3/bbaa123/5868070)), but it is slower than the BCFTools caller.

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -2 20_S30_L001_R2_001.fastq.gz \
    -p nCov-2019-alt_primers.csv \
    -r Wuhan-Hu-1-compele-genome.fasta \
    -t 4 \
    -o my_outdir \
    --lofreq /home/user/Sotf/lofreq_star-2.1.2/bin/lofreq
```

### Custom options for consensus-highlighter

You can pass additional options for [consensus-highlighter](https://github.com/masikol/consensus-highlighter) with option `--highlighter-args`. The options should be surrounded by quotation marks (see Example below).

**Beware:** Oarsman does not check these options, simply passes them to the highlighter.

With that said:

- you may specify the following options: `-c/--coverage-thresholds`, `-n/--no-zero-output`, `--circular`, `--organism`;
- do not specify the following highlighter options with `--highlighter-args`: `-f/--target-fasta`, `-b/--bam`, `-o/--outfile`.

#### Example

```
python3 oarsman.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -2 20_S30_L001_R2_001.fastq.gz \
    -p nCov-2019-alt_primers.csv \
    -r Wuhan-Hu-1-compele-genome.fasta \
    -t 4 \
    -o my_outdir \
    --highlighter-args '-c 5,8 --organism Sabaka -n --circular'
```

If your desired organism name contains spaces, just use double quotes:

```
... \
    --highlighter-args '-c 5,8 --organism "Sabaka dvorii" -n --circular'
```
