#!/bin/bash -login
#SBATCH --mem=100GB
#SBATCH --job-name=scaden
#SBATCH --output=%x-%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --constrain=intel18

################################################################## 
# Print parameters
##################################################################

start=`date +%s`

echo "bulk.mtx.path:" $1
echo "ref.data.dir:" $2
echo "ref.data.pattern:" $3
echo "results.dir:" $4

bulk=$1
ref_dir=$2
ref_pattern=$3
res_dir=$4

################################################################## 
# Run scaden
##################################################################

module purge
module load Conda/3

# scaden simulate \
#   --out $ref_dir \
#   --data $ref_dir \
#   --pattern $ref_pattern 

cd $ref_dir

ref_file=$(ls $ref_pattern)
ref_prefix=${ref_file%$ref_pattern}

scaden process ${ref_prefix}.h5ad $bulk

scaden train processed.h5ad

scaden predict \
  --model_dir $ref_dir \
  --outname ${res_dir}/scaden_proportions.txt \
  $bulk




