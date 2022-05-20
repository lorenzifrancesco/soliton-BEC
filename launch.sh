#!/bin/bash
#SBATCH --job-name=BEC3.0
#SBATCH --error=BECdyn2.err
#SBATCH --output=BECdyn2.out
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
