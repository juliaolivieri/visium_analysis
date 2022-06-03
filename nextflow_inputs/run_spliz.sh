#!/bin/bash
#
#SBATCH --job-name=nf-spliz
#SBATCH --output=nf-spliz.%j.out
#SBATCH --error=nf-spliz.%j.err
#SBATCH --time=48:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=4G
date

source ~/.bashrc
conda deactivate
conda activate nf_spliz_env

nextflow run salzmanlab/spliz -c ../visium.config  -latest -r main -resume --dataname V1_Mouse_Brain_Sagittal_Posterior
date
