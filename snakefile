from snakemake.utils import min_version, validate
from pathlib import Path
from os import path

min_version("7.18")

# Define the config file
configfile: "config.yaml"

# Read the aligner and its options from the config file
configfile: "config.yaml"
config = yaml.expand(config)

# Validate config file
validate(config, schema="config_schema.yml")

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]

target_list = [
    "sniffles/" + sample + ".vcf"
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

# minimap2
rule align_minimap2: 
    input:
        fq = rules.concatenate_reads.output.fq_concat,
        ref = config["genome"]

    output:
        sam = path.join("mapping", f"{sample}.sam")

    threads: config["threads"]

    params:
        opts = config["minimap2_opts"]
    
    shell:"""
    minimap2 {params.opts} -y -x map-ont -t {threads} -a --eqx -k 17 -K 5g {input.ref} {input.fq} -o {output.sam}
    """
# ngmlr
rule align_ngmlr: 
    input:
        fq = rules.concatenate_reads.output.fq_concat,
        ref = config["genome"]

    output:
        sam = path.join("mapping", f"{sample}.sam")
    
    threads: config["threads"]

    params:
        opts = config["nglmr_opts"]
    
    shell:"""
    ngmlr -t 4 -r {input.ref} -q {input.fq} -o {output.sam} -t {threads} -x ont -o {output.sam}
    """

# Conditionally invoke the appropriate aligner rule
if config["aligner"] == "minimap2":
    include: align_minimap2
else:
    include: align_ngmlr

# ------------------------------------------------------------

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
# ------------------------------------------------------------

rule sniffles:
    input: 
        bam = rules.sam_to_bam.output.bam

    output:
        vcf = path.join("sniffles", f"{sample}.vcf")

    params:
        sn_opts = config["sniffles_opts"]

    threads: config["threads"]

    shell:"""
        sniffles -i {input.bam} -v {output.vcf} {params.sn_opts} --threads {threads}
        """
