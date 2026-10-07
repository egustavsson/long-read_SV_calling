# Long-read structural variant calling

A Snakemake pipeline for analysing long-read DNA sequencing data, with one sample
per configuration file.

The workflow generates read statistics with **SeqKit**, aligns reads with
**minimap2** or **ngmlr**, creates a sorted and indexed BAM, calculates depth with
**samtools**, calls structural variants with **Sniffles2**, and detects tandem
repeat expansions with **Straglr**.

## Installation

Requires Linux and [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```bash
git clone https://github.com/egustavsson/long-read_SV_calling.git
cd long-read_SV_calling

conda env create -f envs/environment.yml
conda activate long-read_SV_calling
```

Snakemake installs the required tool environments automatically when the pipeline
is run with `--sdm conda`.

## Input and configuration

Required inputs:

- Long-read DNA reads in FASTQ format.
- The matching reference genome in uncompressed FASTA format.

Edit `config.yaml` before running. The main settings are:

```yaml
# Working directory for results
workdir: "/path/to/results"

# Prefix of output files
sample_name: "sample01"

# Reference genome
genome: "/path/to/reference.fa"

# Input FASTQ file or directory
fastq: "/path/to/fastq_directory/"

# Aligner to use
aligner: "minimap2"

# Maximum threads per rule
threads: 32
```

A FASTQ directory is searched recursively for `.fastq`, `.fq`, `.fastq.gz`, and
`.fq.gz` files. The files are sorted and tracked individually by Snakemake.

Relative paths are resolved against the pipeline directory. Keep input FASTQs
outside the results directory, and use a separate working directory for each
sample. Start the pipeline once sequencing and file transfer are complete.

The remaining settings control mapping presets, extra tool options, repeat
annotations and sorting memory. Comments in `config.yaml` describe these options.

## Run the pipeline

Run from the pipeline directory with the Conda environment activated:

```bash
snakemake --cores 32 --sdm conda
```

Replace `32` with the total number of cores available to the run. The `threads`
setting in `config.yaml` controls the maximum requested by individual rules.

### Optional commands

These are useful when checking a configuration or preparing a run. They are not
required before every execution.

```bash
# Preview the jobs without running them
snakemake --cores 32 --sdm conda --dry-run

# Display the shell commands during execution
snakemake --cores 32 --sdm conda --printshellcmds

# Install the tool environments without processing data
snakemake --cores 32 --sdm conda --conda-create-envs-only

# Resume after an interrupted run, rebuilding incomplete outputs
snakemake --cores 32 --sdm conda --rerun-incomplete
```

To use a different configuration file:

```bash
snakemake --configfile configs/my_sample.yaml --cores 32 --sdm conda
```

## Repeat annotation files

### Sniffles tandem-repeat annotations

Sniffles can use an optional BED file of tandem-repeat regions. Human reference
annotations are available from the
[Sniffles annotations directory](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations).

For the GRCh38 file used in the example configuration, download it from the
pipeline directory:

```bash
mkdir -p data

curl --fail --location \
    https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed \
    --output data/human_GRCh38_no_alt_analysis_set.trf.bed
```

Set its location in `config.yaml`:

```yaml
# Optional tandem-repeat annotations for Sniffles
tandem_repeat_region: "data/human_GRCh38_no_alt_analysis_set.trf.bed"
```

Choose annotations matching your reference assembly and chromosome names. To
omit them, set `tandem_repeat_region: ""`. The workflow supplies the
`--tandem-repeats` argument; do not also add it to `sniffles_opts`.

```yaml
# Optional targeted repeat genotyping
straglr_opts: "--loci /absolute/path/to/simple_repeats.bed"
```

Use an absolute path here because Straglr runs inside the results directory.
This loci BED is separate from the Sniffles annotation file.

## Output

Results are written under the configured `workdir`. `<sample>` is the value of
`sample_name`.

| File or directory | Contents |
| --- | --- |
| `qc/<sample>_seqkit_stats.tsv` | Read counts, bases, lengths, N50, quality and GC statistics |
| `qc/<sample>_fastq_inputs.tsv` | Ordered input file list, sizes and modification times |
| `mapping/<sample>.bam` | Coordinate-sorted alignments |
| `mapping/<sample>.bam.bai` | BAM index |
| `coverage/<sample>_depth.tsv` | Per-base depth at covered positions |
| `sniffles/<sample>.vcf` | Structural variant calls |
| `sniffles/<sample>.snf` | Data for subsequent Sniffles cohort calling |
| `straglr/<sample>.straglr.vcf` | Tandem repeat variants |
| `straglr/<sample>.straglr.tsv` | Read-level repeat results |
| `straglr/<sample>.straglr.bed` | Summarized repeat genotypes |
| `reference/` | Reference symlink and FASTA index |
| `logs/` | Logs from each processing step |

The concatenated FASTQ is temporary and removed after alignment and read QC.
Add `--notemp` to the run command to retain it. No intermediate SAM is written.

## Analysis options

- `minimap2_preset` selects the mapping preset. Additional options are supplied
  through `minimap2_opts`; `"--eqx"` is a simple choice. An explicit `-k` setting
  overrides the preset's k-mer length. The workflow supplies uppercase `-Y` for
  soft clipping of supplementary alignments.
- `sort_mem` sets samtools sorting memory **per thread**, rather than total job
  memory.
- `sniffles_opts`, `ngmlr_opts` and `straglr_opts` accept additional analysis
  options. Leave input, output and thread arguments to the workflow. Files named
  within these option strings are not separately tracked by Snakemake.
- Depth uses the default samtools filters and omits zero-depth positions; the
  output is not a genome-wide mean coverage summary.

## Software versions

The environment files pin the main tools to the following versions:

| Tool | Version |
| --- | --- |
| Snakemake | 9.27.0 |
| SeqKit | 2.10.1 |
| minimap2 | 2.31 |
| ngmlr | 0.2.7 |
| samtools | 1.24 |
| Sniffles | 2.8.1 |
| Straglr | 1.5.6 |

## Tool documentation

- [SeqKit](https://bioinf.shenwei.me/seqkit/usage/)
- [Minimap2](https://github.com/lh3/minimap2)
- [NGMLR](https://github.com/philres/ngmlr)
- [Samtools](https://www.htslib.org/doc/samtools.html)
- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [Straglr](https://github.com/BirolLab/straglr)
- [Snakemake](https://snakemake.readthedocs.io/)
