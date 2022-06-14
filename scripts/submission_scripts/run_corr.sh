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
#DATANAME="cta_ucsf-1-5_liver"
#DATANAME="p20190-s003_3_BrainMetastasis"
#DATANAME="p20190-s004_4_BrainMetastasis"
#DATANAME="p20218-s001_L1"
#DATANAME="p20218-s002_L2"
#DATANAME="p20218-s003_L3"
#DATANAME="p20218-s004_L4"


#DATANAME="p20190-s003_3_BrainMetastasis_b${BOUND}"

#SPLIZPATH="/oak/stanford/groups/horence/JuliaO/data/visium/${DATANAME}/${DATANAME}_all_5000_v2.zscore"
#COL="z_scaled"
#THRESH=100
#GENECOL="window"

#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#COL="SpliZ_resid"
#THRESH=10

DATANAME="p20190_s003_3_BrainMetastasis"
COL="ReadZS"
THRESH=10


a="python ../pixel_correlation.py --dataname ${DATANAME}  --thresh ${THRESH} --score ${COL} --outpath ../output/pixel_correlation/"
echo $a
eval $a
date
