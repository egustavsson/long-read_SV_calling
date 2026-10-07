# Long-read structural variant calling

Development toward **v0.2.0**. The historical workflow remains available as the
GitHub release **v0.1.0**.

This single-sample Snakemake workflow concatenates long-read DNA FASTQs, computes
SeqKit read statistics, maps with minimap2 or ngmlr directly into a sorted BAM,
indexes the BAM, calculates per-base depth, calls structural variants with
Sniffles2, and detects tandem repeat expansions with Straglr.

## Installation

Requires Linux and Conda. Run commands from the repository root.

```bash
git clone https://github.com/egustavsson/long-read_SV_calling.git
cd long-read_SV_calling
conda env create -f envs/environment.yml
conda activate long-read_SV_calling
```

The base environment contains Snakemake. Tool environments are installed per
rule using `--sdm conda`. These YAML files pin the primary tools, but are not
complete locks of every transitive dependency and build.

## Configuration and execution

Edit `config.yaml`; its paths are examples from the previous workflow, not
portable defaults. Use **a new workdir** when validating this update. Paths that
are relative are resolved against the repository root, including paths supplied
in another config file. One sample belongs to each workdir.

`fastq` accepts a single file or a directory searched recursively. Supported
extensions are `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`. Files are sorted by
absolute path and tracked individually by Snakemake. Run after sequencing and
file transfer have finished. Keep input FASTQs outside workdir to avoid collecting
the workflow's own processed reads on the next invocation.

`genome` must be an uncompressed FASTA. A symlink and FASTA index are created in
workdir under `reference/`; the original reference directory need not be writable.
`tandem_repeat_region` accepts a BED matching the reference assembly, or an empty
string to disable tandem-repeat annotations in Sniffles. Do not put
`--tandem-repeats` in `sniffles_opts`.

`threads` sets the maximum requested per computational rule; `--cores` sets the
total scheduler budget. SeqKit uses up to eight threads, indexing/depth up to
four, and concatenation one. Alignment divides its allocation between the
aligner and sorting. At a one-core allocation both piped processes still have a
main thread. `sort_mem` is the memory limit **per sorting thread**, not total
workflow memory. Reference indexing and mapping also need memory; monitor a
representative sample before running many jobs concurrently.

```bash
# Inspect the planned work and rendered commands.
snakemake --cores 32 --sdm conda --dry-run --printshellcmds

# Install the rule environments before running data.
snakemake --cores 32 --sdm conda --conda-create-envs-only

# Run.
snakemake --cores 32 --sdm conda --printshellcmds
```

To use a separate config:

```bash
snakemake --configfile configs/my_sample.yaml --cores 32 --sdm conda
```

Extra tool options remain shell fragments: quote argument values when needed.
Use these fields for analysis options, not for overriding the workflow's input,
output, reference, or thread arguments. Files embedded in extra option strings
are not tracked automatically as dependencies.

## Mapping and repeat calling

The supplied configuration keeps the original minimap2 `map-ont`, `-k 17`,
`-K 5g`, `--eqx`, and `-y` options. The workflow additionally supplies `-Y` for
soft-clipped supplementary alignments, as recommended by Straglr. Presets are
configurable; if selecting another preset, review/remove the explicit `-k`
override. Updated tool versions and `-Y` can change alignments and calls, so
compare a known sample against v0.1.0 before using this as a production release.

Straglr now uses its upstream `straglr.py BAM FASTA PREFIX --nprocs N` interface
and tracks its VCF, TSV, and BED outputs. The former rule declared filtered VCF,
TRF BED, genotype, insertion and BAM statistics files that this command does not
produce. Those targets are removed rather than fabricated. This workflow uses
Straglr's default genome-scan mode; options such as `--loci` can be supplied in
`straglr_opts` for targeted genotyping (such files are not tracked separately).

## Outputs

| Location within workdir | Contents |
| --- | --- |
| `qc/<sample>_seqkit_stats.tsv` | Read count, bases, lengths, N50, quality and GC statistics |
| `qc/<sample>_fastq_inputs.tsv` | Ordered input paths, byte sizes and modification times |
| `mapping/<sample>.bam` | Coordinate-sorted alignments |
| `mapping/<sample>.bam.bai` | Tracked BAM index, rebuilt separately if missing |
| `coverage/<sample>_depth.tsv` | Samtools per-base depth at covered positions |
| `sniffles/<sample>.vcf` | Structural variant calls |
| `sniffles/<sample>.snf` | Candidate data for later Sniffles cohort calling |
| `straglr/<sample>.straglr.vcf` | Tandem repeat variants |
| `straglr/<sample>.straglr.tsv` | Detailed read-level results |
| `straglr/<sample>.straglr.bed` | Summarized locus genotypes |
| `reference/genome.fa` and `.fai` | Local reference symlink and FASTA index |
| `logs/` | Logs for each processing step |

The concatenated uncompressed FASTQ is temporary and removed once mapping and
SeqKit finish. Use `--notemp` if you want to retain it. No intermediate SAM is
written. Concatenation accepts mixed compressed/plain inputs and adds a trailing
newline between files when needed; it does not repair malformed FASTQ records.
Depth retains the previous `samtools depth` default filters and omission of
zero-depth positions; this is not a genome-wide mean coverage summary.

## Pinned primary tools

| Tool | Version |
| --- | --- |
| Snakemake | 9.27.0 |
| SeqKit | 2.10.1 |
| minimap2 | 2.31 |
| ngmlr | 0.2.7 |
| samtools | 1.24 |
| Sniffles | 2.8.1 |
| Straglr | 1.5.6 |

The environments still need to be solved and exercised on your machine before
tagging v0.2.0. See `VALIDATION.md` for the checks performed on this draft.

## References

- [SeqKit usage](https://bioinf.shenwei.me/seqkit/usage/)
- [Snakemake software deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
- [Minimap2](https://github.com/lh3/minimap2)
- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [Straglr interface and outputs](https://github.com/BirolLab/straglr)
