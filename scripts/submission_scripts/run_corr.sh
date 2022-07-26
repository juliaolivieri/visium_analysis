#!/bin/bash
#
#SBATCH --job-name=corr
#SBATCH --output=../job_output/corr.%j.out
#SBATCH --error=../job_output/corr.%j.err
#SBATCH --time=5:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence,quake
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=40G
date
source ~/.bashrc
conda deactivate
conda activate jup_env

#DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior_Section_2"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior_Section_2"
#DATANAME="V1_Mouse_Kidney"
#DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190_s004_4_BrainMetastasis"
#DATANAME="p20218_s001_L1"
#DATANAME="p20218_s002_L2"
#DATANAME="p20218_s003_L3"
DATANAME="p20218_s004_L4"
#DATANAME="cta_ucsf-1-5_liver"

#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#COL="SpliZ_resid"
#THRESH=10

COL1="ReadZS"
COL2="ReadZS_norm"
COL3="ReadZS_ge_norm"
COL4="ReadZS_resid"

THRESH=10


a="python ../pixel_correlation.py --dataname ${DATANAME}  --thresh ${THRESH} --score ${COL1} --outpath ../output/pixel_correlation/"
b="python ../pixel_correlation.py --dataname ${DATANAME}  --thresh ${THRESH} --score ${COL2} --outpath ../output/pixel_correlation/"
c="python ../pixel_correlation.py --dataname ${DATANAME}  --thresh ${THRESH} --score ${COL3} --outpath ../output/pixel_correlation/"
d="python ../pixel_correlation.py --dataname ${DATANAME}  --thresh ${THRESH} --score ${COL4} --outpath ../output/pixel_correlation/"


echo $a
eval $a

echo $b
eval $b

echo $c
eval $c

echo $d
eval $d


date
