#!/bin/bash
#
#SBATCH --job-name=spat_plot
#SBATCH --output=../job_output/spat_plot.%j.out
#SBATCH --error=../job_output/spat_plot.%j.err
#SBATCH --time=1:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=60G

DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior_Section_2"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior_Section_2"
#DATANAME="V1_Mouse_Kidney"
#DATANAME="Visium_FFPE_Mouse_Brain"
#DATANAME="Visium_FFPE_Mouse_Brain_IF"
#DATANAME="Visium_FFPE_Mouse_Kidney"

#DATANAME="Visium_FFPE_Human_Breast_Cancer"
#DATANAME="Visium_FFPE_Human_Normal_Prostate"
#DATANAME="Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma"
#DATANAME="Visium_FFPE_Human_Prostate_Cancer"
#DATANAME="Visium_FFPE_Human_Prostate_IF"
#DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190_s004_4_BrainMetastasis"
#DATANAME="p20218_s001_L1"
#DATANAME="p20218_s002_L2"
#DATANAME="p20218_s003_L3"
#DATANAME="p20218_s004_L4"


#SCORE="ReadZS"
#SCORE2="ReadZS_ge"
#WINDOWFILE="/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/ising/${DATANAME}_ReadZS_norm_100_100_plot.txt"

SCORE="SpliZ_norm"
SCORE2="ge"
WINDOWFILE="/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/ising/${DATANAME}_SpliZ_norm_100_100_b0_plot.txt"


#WINDOWFILE="/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/ising/${DATANAME}_ReadZS_norm_100_100_plot.txt"

date

a="python ../plot_gene_val.py --dataname ${DATANAME} --window_file ${WINDOWFILE} --score ${SCORE} --score2 ${SCORE2}"

echo $a
eval $a

date

