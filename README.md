# Long-read structural variant calling

This is a `snakemake` pipeline that takes long-read DNA sequencing data (fastq) as input, generates fastq stats using `nanostat`, map the reads to the genome using `minimap2` and uses `sniffles` for calling structural variants.

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
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate long-read_SV_calling`.

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
|--- config.yml           # a copy of the parameters used in the pipeline  
|--- Nanostat/  
     |-- # output of nanostat - fastq stats   
|--- Mapping/  
     |-- # output of minimap2 - aligned reads  
|--- sniffles/  
     |-- # output of sniffles - vcf with SVs
```
