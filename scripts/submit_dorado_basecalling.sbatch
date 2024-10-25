#!/bin/bash
#SBATCH --job-name=dorado_basecalling  # Job name
#SBATCH --cpus-per-task=200             # Number of CPU cores
#SBATCH --mem=400G                      # Memory

# Activate the Conda environment
conda activate long-read_SV_calling

# Run dorado_basecalling.sh with the specified parameters
/home/JHEmilGustavsson/MinaRyten/Emil/Snakemakes/long-read_SV_calling/scripts/dorado_basecalling.sh -i "/home/MinaRyten/Emil/projects/SCA4/data/230918_SCA4_run1/SCA4/20230918_1006_3D_PAQ85586_4d32fc77/pod5_pass/barcode11,/home/MinaRyten/Emil/projects/SCA4/data/230918_SCA4_run1/SCA4/20230918_1006_3D_PAQ85586_4d32fc77/pod5_pass/barcode15" -o "/home/MinaRyten/Emil/projects/SCA4/results/dorado_bams/" -s 23712
