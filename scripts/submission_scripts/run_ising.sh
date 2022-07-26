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

#DATANAME="Visium_FFPE_Human_Breast_Cancer"
#DATANAME="Visium_FFPE_Human_Normal_Prostate"
#DATANAME="Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma"
#DATANAME="Visium_FFPE_Human_Prostate_Cancer"
#DATANAME="Visium_FFPE_Human_Prostate_IF"
#DATANAME="Visium_FFPE_Mouse_Brain"
#DATANAME="Visium_FFPE_Mouse_Brain_IF"
#DATANAME="Visium_FFPE_Mouse_Kidney"

#DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190_s004_4_BrainMetastasis"
#DATANAME="p20218_s001_L1"
#DATANAME="p20218_s002_L2"
#DATANAME="p20218_s003_L3"
DATANAME="p20218_s004_L4"
#DATANAME="cta_ucsf-1-5_liver"


#SCORE="SpliZ_norm"
#SCORE="ge_norm"
#SCORE="SpliZ_resid"
#SCORE="ReadZS"
#SCORE="SpliZ"
#THRESH=5

#SCORE="ReadZS_norm"
#SCORE="ReadZS_ge_norm"
#SCORE="ReadZS"
#SCORE="ReadZS_resid"
THRESH=100
NUMPERMS=100

#COL1="ReadZS"
#COL2="ReadZS_norm"
#COL3="ReadZS_ge_norm"
#COL4="ReadZS_resid"

COL1="SpliZ"
COL2="SpliZ_norm"
COL3="ge_norm"
COL4="SpliZ_resid"

#COL5="ReadZS_ge"


a="python ../ising.py --dataname ${DATANAME} --score ${COL1}  --thresh ${THRESH} --num_perms ${NUMPERMS} --suff _b0"
b="python ../ising.py --dataname ${DATANAME} --score ${COL2}  --thresh ${THRESH} --num_perms ${NUMPERMS} --suff _b0"
c="python ../ising.py --dataname ${DATANAME} --score ${COL3}  --thresh ${THRESH} --num_perms ${NUMPERMS}"
d="python ../ising.py --dataname ${DATANAME} --score ${COL4}  --thresh ${THRESH} --num_perms ${NUMPERMS} --suff"
#e="python ../ising.py --dataname ${DATANAME} --score ${COL5}  --thresh ${THRESH} --num_perms ${NUMPERMS}"


echo $a
eval $a

echo $b
eval $b

#echo $c
#eval $c
#
echo $d
eval $d

#echo $e
#eval $e

date
