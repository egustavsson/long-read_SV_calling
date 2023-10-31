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

# rule depth:
#     input: 
    
#     output:
    
#     conda: "envs/depth.yml"

#     threads:

# Call SVs ----------------------------------------------------------------

rule sniffles:
    input: 
        bam = rules.sam_to_bam.output.bam
        STR_bed = config["tandem_repeat_region"]

    output:
        vcf = path.join("sniffles", f"{sample}.vcf")

    params:
        sn_opts = config["sniffles_opts"]

    threads: config["threads"]

    shell:"""
        sniffles -i {input.bam} -v {output.vcf} {params.sn_opts} {input.STR_bed} --threads {threads}
        """
