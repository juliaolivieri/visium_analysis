#!/bin/bash
#
#SBATCH --job-name=ising
#SBATCH --output=../job_output/ising.%j.out
#SBATCH --error=../job_output/ising.%j.err
#SBATCH --time=24:00:00
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

DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190-s004_4_BrainMetastasis"
#DATANAME="p20218-s001_L1"
#DATANAME="p20218-s002_L2"
#DATANAME="p20218-s003_L3"
#DATANAME="p20218-s004_L4"
#DATANAME="cta_ucsf-1-5_liver"


#SCORE="SpliZ_norm"
#SCORE="ge_norm"
#SCORE="SpliZ_resid"
SCORE="ReadZS"
#SCORE="SpliZ"
THRESH=5

#SCORE="ReadZS_norm"
#SCORE="ReadZS_ge_norm"
#SCORE="ReadZS"
#SCORE="ReadZS_resid"
#THRESH=1000

NUMPERMS=10



a="python ../ising.py --dataname ${DATANAME} --score ${SCORE}  --thresh ${THRESH} --num_perms ${NUMPERMS}"
echo $a
eval $a
date
