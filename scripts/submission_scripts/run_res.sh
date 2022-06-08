#!/bin/bash
#
#SBATCH --job-name=res
#SBATCH --output=../job_output/res.%j.out
#SBATCH --error=../job_output/res.%j.err
#SBATCH --time=24:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=20G

date
source ~/.bashrc
conda activate jup_env

#DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior_Section_2"
DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior_Section_2"
#DATANAME="V1_Mouse_Kidney"
#DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190-s004_4_BrainMetastasis"
#DATANAME="p20218-s001_L1"
#DATANAME="p20218-s002_L2"
#DATANAME="p20218-s003_L3"
#DATANAME="p20218-s004_L4"
#DATANAME="cta_ucsf-1-5_liver"

SCORE="SpliZ"
SCORE2="ge"
THRESH=100

#SCORE="ReadZS"
#SCORE2="ReadZS_ge"
#THRESH=100


OUTNAME="../output/save_residuals/${DATANAME}_${SCORE}_${SCORE2}_${THRESH}.tsv"



a="python ../save_residuals.py --dataname ${DATANAME} --score ${SCORE} --score2 ${SCORE2}  --thresh ${THRESH} --outname ${OUTNAME}"
echo $a
eval $a
date
