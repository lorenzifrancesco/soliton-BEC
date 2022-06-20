#!/bin/bash
#SBATCH --job-name=BECdyn1.0
#SBATCH --error=BEC3D.err
#SBATCH --output=BEC3D.out
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user francesco.lorenzi98@gmail.com
#SBATCH --mail-type ALL
#SBATCH --gres=gpu:1

cd $WORKING_DIR
#your working directory
cd soliton-BEC

srun singularity exec --nv bec.sif julia --threads=auto --project=. test/test_3d_transmission.jl