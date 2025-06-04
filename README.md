# Long-read structural variant calling

This is a `snakemake` pipeline that takes long-read DNA sequencing data (fastq) as input, generates fastq stats using `nanostat`, maps the reads to the genome using either `minimap2` or `ngmlr`, calculates depth and coverage using `samtools`, calls structural variants using `sniffles2`, and detects tandem repeat expansions using `straglr`.

# Getting Started

## Input

- fastq reads
- Reference genome assembly in fasta format

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (including snakemake) are installed via conda through the `environment.yml` file

Clone the directory:

```bash
git clone --recursive https://github.com/egustavsson/long-read_SV_calling.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd long-read_SV_calling
conda env create -f ./envs/environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. There is also the option of which aligner to use, default is `minimap2`. 

`snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate long-read_SV_calling`.

```bash
cd long-read_SV_calling
conda activate long-read_SV_calling
snakemake --use-conda -j <num_cores> all
```
It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

## Output
```
working directory  
|--- config.yml                      # parameters used  
|--- processed_reads/  
     |-- <sample>_reads.fq          # concatenated reads  
|--- Nanostat/  
     |-- <sample>_stat_out.txt      # output of NanoStat  
|--- mapping/  
     |-- <sample>.bam               # aligned reads  
     |-- <sample>.bam.bai           # BAM index  
|--- coverage/  
     |-- <sample>_depth.tsv         # coverage info from samtools  
|--- sniffles/  
     |-- <sample>.vcf               # structural variant calls  
|--- straglr/  
     |-- <sample>.straglr.vcf                 # tandem repeat calls  
     |-- <sample>.straglr.filtered.vcf       # filtered TRs  
     |-- <sample>.straglr.trf.bed            # regions detected by TRF  
     |-- <sample>.straglr.genotype.txt       # per-locus genotype  
     |-- <sample>.straglr.insertions.txt     # repeat insertions  
     |-- <sample>.straglr.bamstats.txt       # coverage stats
```

## Versions

```
snakemake  : 7.18.2
nanostat   : 1.6
minimap2   : 2.26
ngmlr      : 0.2.7
samtools   : 1.09
sniffles   : 2.0.7
straglr    : 1.5.3
```
