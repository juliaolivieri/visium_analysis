#!/bin/bash
#
#SBATCH --job-name=ge
#SBATCH --output=../job_output/ge.%j.out
#SBATCH --error=../job_output/ge.%j.err
#SBATCH --time=2:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence,quake
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=120G
date
source activate jup_env

#DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Anterior_Section_2"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior_Section_2"
#DATANAME="V1_Mouse_Kidney"
#MATPATH="/oak/stanford/groups/horence/JuliaO/data/visium/${DATANAME}/filtered_feature_bc_matrix/"

#DATANAME="cta_ucsf-1-5_liver"
#MATPATH="/oak/stanford/groups/horence/JuliaO/data/visium/Biohub_covid_liver/cta_ucsf-1-5_liver/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
#FEATPATH="/oak/stanford/groups/horence/JuliaO/data/visium/Biohub_covid_liver/cta_ucsf-1-5_liver/outs/filtered_feature_bc_matrix/features.tsv.gz"
#BARPATH="/oak/stanford/groups/horence/JuliaO/data/visium/Biohub_covid_liver/cta_ucsf-1-5_liver/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

#DATANAME="p20190_s003_3_BrainMetastasis"
#DATANAME="p20190_s004_4_BrainMetastasis"
#DATANAME="p20218_s003_L3"
#DATANAME="p20218_s002_L2"
#DATANAME="p20218_s004_L4"
#DATANAME="p20218_s001_L1"
#DATANAME="Visium_FFPE_Human_Breast_Cancer"
#DATANAME="Visium_FFPE_Human_Normal_Prostate"
#DATANAME="Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma"
#DATANAME="Visium_FFPE_Human_Prostate_Cancer"
#DATANAME="Visium_FFPE_Human_Prostate_IF"
#DATANAME="Visium_FFPE_Mouse_Brain"
#DATANAME="Visium_FFPE_Mouse_Brain_IF"
DATANAME="Visium_FFPE_Mouse_Kidney"

#MATPATH="/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/${DATANAME}/outs/filtered_feature_bc_matrix/matrix.mtx.gz

#"
#FEATPATH="/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/${DATANAME}/outs/filtered_feature_bc_matrix/features.tsv.gz"
#BARPATH="/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/${DATANAME}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"



#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"

THRESH=1000
#THRESH=0
#GENES="--genes Usp9y Kdm5d"
GENES=""

a="python ../parse_gene_expression.py --dataname ${DATANAME}  --thresh ${THRESH} ${GENES} --norm "
echo $a
eval $a
date
