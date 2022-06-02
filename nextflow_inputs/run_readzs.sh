#!/bin/bash
#
#SBATCH --job-name=nf-readzs
#SBATCH --output=nf-readzs.%j.out
#SBATCH --error=nf-readzs.%j.err
#SBATCH --time=48:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4G
date
#source activate /share/PI/horence/applications/anaconda3/envs/spliz_envR2
#a="make mini-slurm-nosing"
#ml load java
#ml load biology
#ml load samtools
#ml load R

source ~/.bashrc
conda deactivate
conda activate readzs_env

nextflow run /oak/stanford/groups/horence/JuliaO/ReadZS/main.nf \
	-c visium_readzs.config \
	--runName V1_Mouse_Brain_Sagittal_Posterior \
	-resume


date
