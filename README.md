# ChIP-seq-Pipeline-with-Spike-in-Normalization-S.-pombe-S.-cerevisiae-

**Organisms:** *Schizosaccharomyces pombe* (target) and *Saccharomyces cerevisiae* (spike-in)  
**Workflow:** Quality control → adapter trimming → alignment → duplicate marking → spike-in scaling → coverage tracks → peak calling → relative enrichment → metagene analysis

---

## 1. Overview

This Bash pipeline processes paired-end ChIP-seq FASTQ files. It aligns the reads separately against:

- ***S. pombe***: the target/reference genome used for the main ChIP-seq analysis.
- ***S. cerevisiae***: the spike-in reference used to estimate a sample-specific scaling factor.

The calculated spike-in scale factor is applied to the *S. pombe* BAM file to generate a spike-in-normalized bigWig coverage track. The pipeline also generates an RPM-normalized track, calls peaks with MACS2, filters peaks using a *S. pombe* blacklist and fold-enrichment threshold, calculates a ChIP/Input ratio track, and creates metagene profiles.

> **Important:** This README documents the commands and assumptions in the supplied script. It does not imply that the pipeline has been validated on every dataset. Review the checks and methodological notes in Sections 7–9 before interpreting results.

## 2. Workflow

1. Run FastQC on the raw paired-end FASTQ files.
2. Summarize pre-trimming QC reports with MultiQC.
3. Trim adapters and low-quality bases using Trim Galore.
4. Align trimmed reads to the *S. pombe* reference using BWA-MEM.
5. Sort, fix mate information, mark duplicates, and index the *S. pombe* BAM.
6. Align the same trimmed reads separately to the *S. cerevisiae* spike-in reference.
7. Count mapped spike-in alignments and calculate the scale factor:

   `Spike-in scale factor = 1,000,000 / mapped spike-in reads`

8. Generate a spike-in-scaled *S. pombe* bigWig track using deepTools `bamCoverage`.
9. Generate a separate RPM-normalized bedGraph and bigWig track.
10. Call peaks using MACS2 with an Input control BAM.
11. Remove peaks overlapping the *S. pombe* blacklist and retain peaks with fold enrichment greater than 2.
12. Calculate a ChIP/Input relative enrichment ratio (rER) bigWig using `bigwigCompare`.
13. Generate a metagene matrix and profile using deepTools.

## 3. Software requirements

The script expects the following command-line tools to be installed and available in `PATH`:

| Tool | Purpose |
|---|---|
| BWA | Read alignment |
| SAMtools | BAM conversion, sorting, mate fixing, duplicate marking, indexing and counting |
| FastQC | Raw-read quality control |
| MultiQC | QC report aggregation |
| Trim Galore | Adapter and quality trimming |
| deepTools | Coverage normalization, bigWig comparison and metagene analysis |
| BEDTools | Genome coverage and blacklist overlap filtering |
| MACS2 | ChIP-seq peak calling |
| UCSC `bedGraphToBigWig` | Convert bedGraph to bigWig |
| Bash, awk, cut | Pipeline execution and text processing |

Record the versions used in the project or analysis report. The script does not install software or manage environments.

## 4. Required input files and directory layout

The script assumes that the working directory contains the following structure:

```text
project/
├── raw_reads/
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── refs/
│   ├── pombe/
│   │   ├── ASM294v3.fa
│   │   ├── pombase_blacklist.bed
│   │   └── pombe_genes.bed
│   └── cerevisiae/
│       └── R64.fa
├── bams/
│   └── input_markdup.bam
├── tracks/
│   └── input_spike_norm.bw
└── run_chipseq_pipeline.sh
```

The pipeline creates these output directories automatically:

```text
qc/
trimmed/
alignments/
bams/
tracks/
peaks/
metagene/
```

### Input file descriptions

- `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`: paired-end raw ChIP-seq reads. The names must match the configured `PREFIX`.
- `ASM294v3.fa`: *S. pombe* reference genome FASTA.
- `R64.fa`: *S. cerevisiae* spike-in reference genome FASTA.
- `pombase_blacklist.bed`: blacklist regions for *S. pombe*, in BED format and compatible with the reference assembly.
- `pombe_genes.bed`: gene annotation regions for metagene analysis, in BED format and compatible with the *S. pombe* bigWig.
- `input_markdup.bam`: preprocessed Input control BAM used by MACS2. It should be aligned to the same *S. pombe* reference and use a compatible paired-end BAM format.
- `input_spike_norm.bw`: precomputed spike-in-scaled Input control bigWig used as the denominator in the rER calculation.

**Reference preparation:** BWA index files must be generated for both reference FASTA files before running the script (for example, with `bwa index`). Ensure the reference FASTA, blacklist, gene BED, BAM and bigWig all refer to the same *S. pombe* genome assembly. The Input BAM must also be indexed if required by downstream tools or your local workflow.

## 5. Configuration

At the beginning of the Bash script, edit the following variables for each sample:

```bash
THREADS=8
PREFIX="sample1"
R1_RAW="raw_reads/${PREFIX}_R1.fastq.gz"
R2_RAW="raw_reads/${PREFIX}_R2.fastq.gz"

REF_POMBE="refs/pombe/ASM294v3.fa"
REF_CEREVISIAE="refs/cerevisiae/R64.fa"
BLACKLIST_POMBE="refs/pombe/pombase_blacklist.bed"

INPUT_BAM="bams/input_markdup.bam"
INPUT_SPIKE_SCALED_BW="tracks/input_spike_norm.bw"
```

- Change `THREADS` according to the available CPU resources.
- Change `PREFIX` to the sample name. The FASTQ files must follow the naming convention in the script, or update `R1_RAW` and `R2_RAW`.
- Update reference and control paths if your project uses different files.
- `GENES_BED` is set later in the script to `refs/pombe/pombe_genes.bed`; update it if needed.

## 6. Running the pipeline

Save the supplied script as `run_chipseq_pipeline.sh` in the project directory.

```bash
chmod +x run_chipseq_pipeline.sh
./run_chipseq_pipeline.sh
```

The script uses:

```bash
set -euo pipefail
```

This makes the script exit when a command fails, an unset variable is used, or a command in a pipeline fails. Check the terminal output to identify the command that stopped execution.

## 7. Output files

| Directory | Main output | Description |
|---|---|---|
| `qc/` | FastQC reports | Raw-read quality reports |
| `qc/multiqc_pretrim/` | MultiQC report | Aggregated pre-trimming QC report |
| `trimmed/` | `*_R1_val_1.fq.gz`, `*_R2_val_2.fq.gz` | Trimmed paired-end reads |
| `bams/` | `*_pombe_markdup.bam` and `.bai` | Coordinate-sorted *S. pombe* BAM with duplicates marked |
| `alignments/` | `*_cerevisiae.bam` | Reads aligned to the spike-in reference |
| `tracks/` | `*_spike_norm.bw` | Spike-in-scaled target coverage track |
| `tracks/` | `*_rpm.bw` | RPM-scaled target coverage track |
| `tracks/` | `*_rER.bw` | ChIP/Input ratio track |
| `peaks/` | `*_peaks.narrowPeak` and `*_summits.bed` (MACS2 outputs) | MACS2 peak-calling outputs |
| `peaks/` | `*_filtered_enrichment_gt2.bed` | Peaks not overlapping the blacklist and with fold enrichment > 2 |
| `metagene/` | `*_matrix.gz` | deepTools computeMatrix output |
| `metagene/` | `*_profile.pdf` | Metagene profile plot |

Intermediate *S. pombe* name-sorted, fixmate and coordinate-sorted BAM files are deleted by the script after duplicate marking. The *S. cerevisiae* alignment BAM is retained.

## 8. Key parameters used

| Step | Parameter | Meaning |
|---|---|---|
| FastQC | `-t 8` | Uses 8 threads |
| Trim Galore | `--paired --phred20` | Paired-end trimming with Phred quality threshold 20 |
| BWA-MEM | `-t 8` | Uses 8 threads for alignment |
| `bamCoverage` | `--binSize 10` | Coverage bins of 10 bp |
| Spike-in scaling | `1e6 / SPIKE_READS` | Scale factor calculated from mapped spike-in alignments |
| RPM scaling | `1e6 / PROPER_READS` | Scale factor based on the script's target BAM read count |
| MACS2 | `-f BAMPE` | Paired-end BAM mode |
| MACS2 | `-g 1.25e7` | Effective genome size set to approximately 12.5 Mb |
| MACS2 | `-q 0.05` | q-value cutoff |
| Peak filtering | `$7 > 2` | Retains narrowPeak entries with fold enrichment > 2 |
| `bigwigCompare` | `--operation ratio` | Calculates ChIP/Input signal ratio |
| `bigwigCompare` | `--binSize 10 --smoothLength 100` | 10 bp bins and 100 bp smoothing |
| `computeMatrix` | 1,000 bp upstream; 3,000 bp body; 1,500 bp downstream | Scaled gene-region profile |
| `computeMatrix` | `--binSize 10` | 10 bp matrix bins |

## 9. Important methodological checks before interpreting results

These points are important because they affect the biological interpretation and reproducibility of the results.

1. **Separate versus competitive alignment:** The script aligns the same reads independently to the *S. pombe* and *S. cerevisiae* references. This is not a combined competitive alignment. Reads that can align to both genomes may therefore be counted in both alignments. Consider whether a combined-reference alignment and species-specific assignment are more appropriate for your experimental design.

2. **Spike-in read counting:** The command `samtools view -c -F 4` counts mapped alignment records; it does not, by itself, guarantee that reads are uniquely mapped or count only primary alignments. If the intended normalization is based on uniquely mapped spike-in fragments, define and apply the appropriate mapping-quality and alignment-flag filters. For paired-end data, confirm whether the intended unit is reads, read pairs or fragments.

3. **Duplicate handling:** `samtools markdup` marks duplicates; it does not remove them. Confirm whether duplicates should be retained or excluded for your assay and downstream analysis. The current `bamCoverage` command does not explicitly request duplicate exclusion.

4. **Spike-in scaling convention:** The script uses `1,000,000 / mapped spike-in alignments` as the scale factor. Confirm that this convention matches the spike-in amount, read-counting method and normalization strategy used in the experiment. Use a consistent method for all ChIP and Input samples.

5. **Input normalization:** `input_spike_norm.bw` must have been generated using a compatible reference, bin size, spike-in counting method and scaling convention. Otherwise, the ChIP/Input ratio may not be directly comparable.

6. **rER interpretation:** `bigwigCompare --operation ratio` calculates a ratio of the two bigWig signals. It is not a general background-subtraction procedure. Consider how zero or very low Input signal bins are handled and whether a pseudocount or another comparison method is appropriate.

7. **Peak calling normalization:** MACS2 receives the target BAM and Input BAM, but the command does not apply the spike-in scale factor to peak calling. The spike-in scaling is applied to the coverage track, not directly to MACS2 peak calling in this script.

8. **RPM count definition:** `samtools view -c -F 260` excludes unmapped and secondary alignments, but the exact interpretation of the resulting count should be checked for paired-end data and the desired fragment-counting convention. It is not necessarily a count of unique fragments.

9. **Peak filtering:** In a MACS2 narrowPeak file, column 7 is fold enrichment. The pipeline retains values strictly greater than 2 and removes peaks overlapping the supplied blacklist. Confirm the blacklist assembly and coordinate compatibility.

10. **Metagene input:** The gene BED file must use coordinates and chromosome names matching the *S. pombe* bigWig. Check BED coordinate conventions, gene orientation and whether the regions represent the intended gene set.

11. **Zero spike-in reads:** If no spike-in alignments are counted, the scale-factor calculation is invalid. Add a check to stop the pipeline when the spike-in count is zero.

These are review points, not claims that the current script has already implemented these safeguards.

## 10. Suggested quality-control review

Before using the output for biological conclusions, review:

- FastQC and MultiQC reports for read quality, adapter content and sequence duplication.
- The number and proportion of reads aligning to each reference.
- The spike-in mapped-read count and calculated scale factor for every sample.
- Duplicate rates and library complexity.
- ChIP and Input coverage tracks in IGV.
- MACS2 peak numbers, enrichment values and blacklist-filtered peak numbers.
- Consistency of the Input and ChIP bigWig generation methods.
- Metagene profiles for expected signal patterns and potential sample-specific artefacts.

Keep a record of software versions, reference genome assembly/release, sample name, control files, parameters and the date of each run.

## 11. Reproducibility record

For each sample, record at least:

| Field | Record |
|---|---|
| Sample ID | |
| ChIP target / antibody | |
| Biological replicate | |
| Input control ID | |
| *S. pombe* reference assembly | |
| *S. cerevisiae* reference assembly | |
| Raw FASTQ file names | |
| Spike-in mapped-read count | |
| Spike-in scale factor | |
| Software versions | |
| Date of analysis | |
| Notes / deviations | |

---

**Project note:** This pipeline is intended to generate *S. pombe* ChIP-seq coverage and peak-analysis outputs with a *S. cerevisiae* spike-in-derived scaling factor. Document any changes to filtering, read counting, normalization or peak-calling parameters before comparing samples or biological replicates.
