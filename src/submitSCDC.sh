#!/bin/bash -login
#SBATCH --mem=100GB
#SBATCH --job-name=SCDC
#SBATCH --output=%x-%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --constrain=intel18

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript /mnt/research/bioinformaticsCore/projects/Yuzbasiyan-Gurkan/BCC119_vyg_deconvolution/src/SCDC.R