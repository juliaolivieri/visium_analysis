#!/bin/bash
#
#SBATCH --job-name=ge_readzs
#SBATCH --output=../job_output/ge_readzs.%j.out
#SBATCH --error=../job_output/ge_readzs.%j.err
#SBATCH --time=1:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=20G

date
source ~/.bashrc
conda activate jup_env

#DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior_Section_2"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior_Section_2"
#DATANAME="V1_Mouse_Kidney"
#DATANAME="10X_P1_1"

DATANAME="p20190_s003_3_BrainMetastasis"


THRESH=100


a="python ../readzs_ge.py --dataname ${DATANAME}  --thresh ${THRESH}"
echo $a
eval $a
date
