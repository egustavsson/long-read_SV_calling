from pathlib import Path
import shlex
from snakemake.utils import min_version, validate
from snakemake.exceptions import WorkflowError

min_version("9.0")

# Define and validate the config file
configfile: "config.yaml"
validate(config, "config_schema.yml")

# Resolve relative input and workdir paths against the repository, not workdir.
REPO = Path(workflow.snakefile).resolve().parent


def resolve_path(value):
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = REPO / p

    return str(p.resolve())


workdir: resolve_path(config["workdir"])

sample = config["sample_name"]
ref = resolve_path(config["genome"])

if ref.endswith(".gz"):
    raise WorkflowError("genome must be an uncompressed FASTA file")

# Discover the input FASTQ files
fq_path = Path(resolve_path(config["fastq"]))
suffixes = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

if fq_path.is_dir():
    fastqs = sorted(
        str(p.resolve())
        for p in fq_path.rglob("*")
        if p.is_file() and p.name.endswith(suffixes)
    )
elif fq_path.is_file() and fq_path.name.endswith(suffixes):
    fastqs = [str(fq_path)]
else:
    raise WorkflowError(f"fastq must be a FASTQ file or directory: {fq_path}")

if not fastqs:
    raise WorkflowError(f"No FASTQ files found in {fq_path}")

# Do not silently concatenate duplicate paths reached through symlinks.
if len(set(fastqs)) != len(fastqs):
    raise WorkflowError("FASTQ directory contains duplicate symlink targets")

# Optional tandem-repeat annotations for Sniffles
tr_bed = None
if config["tandem_repeat_region"]:
    tr_bed = resolve_path(config["tandem_repeat_region"])

if "--tandem-repeats" in shlex.split(config["sniffles_opts"]):
    raise WorkflowError("Remove --tandem-repeats from sniffles_opts; set tandem_repeat_region instead")


# Final outputs ---------------------------------------------------------

rule all:
    input:
        f"qc/{sample}_seqkit_stats.tsv",
        f"qc/{sample}_fastq_inputs.tsv",
        f"mapping/{sample}.bam",
        f"mapping/{sample}.bam.bai",
        f"coverage/{sample}_depth.tsv",
        f"sniffles/{sample}.vcf",
        f"sniffles/{sample}.snf",
        f"straglr/{sample}.straglr.vcf",
        f"straglr/{sample}.straglr.tsv",
        f"straglr/{sample}.straglr.bed"


# Concatenate reads -----------------------------------------------------

rule concatenate_reads:
    input:
        fastqs

    output:
        fq = temp(f"processed_reads/{sample}_reads.fq"),
        manifest = f"qc/{sample}_fastq_inputs.tsv"

    threads: 1

    log:
        f"logs/{sample}.concatenate.log"

    script:
        "scripts/concatenate_reads.py"


# Read stats ------------------------------------------------------------

rule seqkit_stats:
    input:
        fq = rules.concatenate_reads.output.fq

    output:
        stats = f"qc/{sample}_seqkit_stats.tsv"

    threads: min(8, config["threads"])

    log:
        f"logs/{sample}.seqkit.log"

    conda:
        "envs/seqkit.yml"

    shell:
        """
        seqkit stats --all --tabular -j {threads} {input.fq:q} > {output.stats:q} 2> {log:q}
        """

# Prepare reference -----------------------------------------------------

# Index a local symlink so the original reference directory can be read-only.
rule prepare_reference:
    input:
        ref = ref

    output:
        ref = "reference/genome.fa",
        fai = "reference/genome.fa.fai"

    threads: 1

    log:
        f"logs/{sample}.reference.log"

    conda:
        "envs/samtools.yml"

    shell:
        """
        ln -sf {input.ref:q} {output.ref:q}; samtools faidx {output.ref:q} 2> {log:q}
        """


# Align reads and sort BAM ----------------------------------------------

rule align:
    input:
        fq = rules.concatenate_reads.output.fq,
        ref = rules.prepare_reference.output.ref,
        fai = rules.prepare_reference.output.fai

    output:
        bam = f"mapping/{sample}.bam"

    threads: config["threads"]

    params:
        # samtools -@ counts additional workers; leave room for its main thread.
        sort_workers = lambda wildcards, threads: max(0, min(4, threads // 4) - 1),
        align_threads = lambda wildcards, threads: max(1, threads - max(1, min(4, threads // 4))),
        sort_mem = config["sort_mem"],
        aligner = config["aligner"],
        command = (
            f"minimap2 -a -x {shlex.quote(config['minimap2_preset'])} -Y {config['minimap2_opts']}"
            if config["aligner"] == "minimap2"
            else f"ngmlr -x ont {config['ngmlr_opts']}"
        )

    log:
        align = f"logs/{sample}.align.log",
        sort = f"logs/{sample}.sort.log"

    conda:
        "envs/mapping.yml"

    shell:
        """
        if [ {params.aligner:q} = minimap2 ]; then
            {params.command} -t {params.align_threads} {input.ref:q} {input.fq:q} 2> {log.align:q} |
                samtools sort -@ {params.sort_workers} -m {params.sort_mem:q} -T {output.bam:q}.tmp -o {output.bam:q} - 2> {log.sort:q}
        else
            {params.command} -t {params.align_threads} -r {input.ref:q} -q {input.fq:q} 2> {log.align:q} |
                samtools sort -@ {params.sort_workers} -m {params.sort_mem:q} -T {output.bam:q}.tmp -o {output.bam:q} - 2> {log.sort:q}
        fi
        """


# Index BAM -------------------------------------------------------------

rule index_bam:
    input:
        bam = rules.align.output.bam

    output:
        bai = f"mapping/{sample}.bam.bai"

    threads: min(4, config["threads"])

    params:
        workers = lambda wildcards, threads: max(0, threads - 1)

    log:
        f"logs/{sample}.index.log"

    conda:
        "envs/samtools.yml"

    shell:
        """
        samtools index -@ {params.workers} {input.bam:q} {output.bai:q} 2> {log:q}
        """


# Depth and coverage ----------------------------------------------------

rule depth:
    input:
        bam = rules.align.output.bam,
        bai = rules.index_bam.output.bai

    output:
        depth = f"coverage/{sample}_depth.tsv"

    threads: min(4, config["threads"])

    params:
        workers = lambda wildcards, threads: max(0, threads - 1)

    log:
        f"logs/{sample}.depth.log"

    conda:
        "envs/samtools.yml"

    shell:
        """
        samtools depth -@ {params.workers} {input.bam:q} > {output.depth:q} 2> {log:q}
        """


# Call SVs --------------------------------------------------------------

rule sniffles:
    input:
        bam = rules.align.output.bam,
        bai = rules.index_bam.output.bai,
        ref = rules.prepare_reference.output.ref,
        fai = rules.prepare_reference.output.fai,
        tandem_repeats = [tr_bed] if tr_bed else []

    output:
        vcf = f"sniffles/{sample}.vcf",
        snf = f"sniffles/{sample}.snf"

    threads: config["threads"]

    params:
        opts = config["sniffles_opts"],
        tr_arg = f"--tandem-repeats {shlex.quote(tr_bed)}" if tr_bed else "",
        sample = sample

    log:
        f"logs/{sample}.sniffles.log"

    conda:
        "envs/sniffles.yml"

    shell:
        """
        sniffles --input {input.bam:q} --vcf {output.vcf:q} --snf {output.snf:q} \
            --reference {input.ref:q} --sample-id {params.sample:q} \
            {params.tr_arg} {params.opts} --threads {threads} > {log:q} 2>&1
        """


# Call STRs -------------------------------------------------------------

rule straglr:
    input:
        bam = rules.align.output.bam,
        bai = rules.index_bam.output.bai,
        ref = rules.prepare_reference.output.ref,
        fai = rules.prepare_reference.output.fai

    output:
        vcf = f"straglr/{sample}.straglr.vcf",
        tsv = f"straglr/{sample}.straglr.tsv",
        bed = f"straglr/{sample}.straglr.bed"

    threads: config["threads"]

    params:
        prefix = lambda wildcards, output: str(Path(output.vcf).with_suffix("")),
        opts = config["straglr_opts"]

    log:
        f"logs/{sample}.straglr.log"

    conda:
        "envs/straglr.yml"

    shell:
        """
        straglr.py {input.bam:q} {input.ref:q} {params.prefix:q} --nprocs {threads} {params.opts} > {log:q} 2>&1
        """
