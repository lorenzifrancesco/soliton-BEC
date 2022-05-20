#!/bin/bash
#SBATCH --job-name=BECdyn1.0
#SBATCH --error=BECdyn1.%j.err
#SBATCH --output=BECdyn1.%j.out
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --mail-user francesco.lorenzi98@gmail.com
#SBATCH --mail-type ALL

cd $WORKING_DIR
#your working directory
cd soliton-BEC

srun singularity exec bec.sif julia --threads=auto --project="." test/test_transmission.jl
