from snakemake.utils import min_version, validate
from pathlib import Path
from os import path
import yaml

min_version("7.18")

# Define the config file
configfile:"config.yaml"
workdir: config["workdir"]

# Validate config file
validate(config, schema="config_schema.yml")

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]

target_list = [
    f"sniffles/{sample}.vcf",
    f"coverage/{sample}_depth.tsv",
    f"Nanostat/{sample}_stat_out.txt",
    f"straglr/{sample}.straglr.vcf",
    f"straglr/{sample}.straglr.filtered.vcf",
    f"straglr/{sample}.straglr.trf.bed",
    f"straglr/{sample}.straglr.genotype.txt",
    f"straglr/{sample}.straglr.insertions.txt",
    f"straglr/{sample}.straglr.bamstats.txt"
]


rule all:
    input: 
        target_list

# Concatenate reads -------------------------------------------------------

in_fastq = config["fastq"]
if not path.isabs(in_fastq):
    in_fastq = path.join(SNAKEDIR, in_fastq)
    assert os.path.exists(in_fastq)

rule concatenate_reads:
	input:
		fq = in_fastq
	
	output:
		fq_concat = path.join("processed_reads", f"{sample}_reads.fq")
	
	threads: config["threads"]
	
	shell:"""
	find {input.fq}  -regextype posix-extended -regex '.*\.(fastq.gz|fq.gz)$' -exec zcat {{}} \\; > {output.fq_concat}
	"""
# Read stats --------------------------------------------------------------

rule nanostat:
    input:
        rules.concatenate_reads.output.fq_concat

    output: 
        ns = path.join("Nanostat", f"{sample}_stat_out.txt")

    threads: config["threads"]

    shell: """
        NanoStat -n {output.ns} -t {threads} --tsv --fastq {input.fq}
        """

# Align reads -------------------------------------------------------------

# Config file input on which aligner to use

print("Aligner selected:", config["aligner"])

rule align:
    input:
        fq = rules.concatenate_reads.output.fq_concat,
        ref = config["genome"]
    
    output:
        sam = path.join("mapping", f"{sample}.sam")
    
    threads: config["threads"]
    
    params:
        minimap2_opts = config["minimap2_opts"],
        ngmlr_opts = config["ngmlr_opts"]
    
    run:
        if config["aligner"] == "minimap2":
            shell(
                "minimap2 {params.minimap2_opts} -y -x map-ont -t {threads} -a --eqx -k 17 -K 5g {input.ref} {input.fq} -o {output.sam}"
            )
            
        else:
            shell(
                "ngmlr -r {input.ref} -q {input.fq} -t {threads} {params.ngmlr_opts} -x ont -o {output.sam}"
            )

# Sorted BAM --------------------------------------------------------------

rule sam_to_bam: 
    input:
        sam = lambda wildcards: path.join("mapping", f"{sample}.sam")
    output:
        bam = path.join("mapping", f"{sample}.bam")
    threads: config["threads"]

    shell:"""
        samtools sort -@ {threads} -O BAM -o {output.bam} {input.sam};
        samtools index {output.bam} 
    """
# Depth and coverage ------------------------------------------------------

rule depth:
    input: 
        bam = rules.sam_to_bam.output.bam

    output:
        depth_tsv = path.join("coverage", f"{sample}_depth.tsv")

    threads: config["threads"]

    conda: "envs/depth.yml"

    shell:"""
        samtools depth -@ {threads} {input.bam} > {output.depth_tsv}
        """

# Call SVs ----------------------------------------------------------------

rule sniffles:
    input: 
        bam = rules.sam_to_bam.output.bam

    output:
        vcf = path.join("sniffles", f"{sample}.vcf")

    params:
        sn_opts = config["sniffles_opts"],
	STR_bed = config["tandem_repeat_region"]

    threads: config["threads"]

    shell:"""
        sniffles -i {input.bam} -v {output.vcf} {params.sn_opts} {params.STR_bed} --threads {threads}
        """

# Call STRs --------------------------------------------------------------

rule straglr:
    input:
        bam=rules.sam_to_bam.output.bam,
        ref=config["genome"]
    
    output:
        vcf="straglr/{sample}.straglr.vcf",
        vcf_filt="straglr/{sample}.straglr.filtered.vcf",
        trf="straglr/{sample}.straglr.trf.bed",
        genotype="straglr/{sample}.straglr.genotype.txt",
        insertions="straglr/{sample}.straglr.insertions.txt",
        stats="straglr/{sample}.straglr.bamstats.txt"
    
    params:
        str_opts = config["straglr_opts"]

    conda:
        "envs/straglr.yml"
    
    threads: config["threads"]
    
    shell:"""
        straglr -t {threads} -i {input.bam} {params.str_opts} -r {input.ref} -o straglr/{wildcards.sample}.straglr.vcf
        """