#!/bin/bash -login
#SBATCH --mem=100GB
#SBATCH --job-name=omnideconv
#SBATCH --output=%x-%j.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --constrain=intel18

################################################################## 
# Print parameters
##################################################################

start=`date +%s`
echo "omnideconv.R path:" $1
echo "ref.seurat.path:" $2
echo "bulk.mtx.path:" $3
echo "tool:" $4
echo "annotation.col:" $5
echo "batch.col:" $6
echo "results.dir:" $7

################################################################## 
# Run omnideconv
##################################################################

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $1 $2 $3 $4 $5 $6 $7

################################################################## 
# Finish
##################################################################

echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 
