#!/usr/bin/env bash

#ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20190_s004_4_BrainMetastasis/zscore/"
#DATANAME="p20190_s004_4_BrainMetastasis"

#ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20218_s001_L1/zscore/"
#DATANAME="p20218_s001_L1"

#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"

DATANAME="Visium_FFPE_Human_Breast_Cancer"
#DATANAME="Visium_FFPE_Human_Normal_Prostate"
#DATANAME="Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma"
#DATANAME="Visium_FFPE_Human_Prostate_Cancer"
#DATANAME="Visium_FFPE_Human_Prostate_IF"
#DATANAME="Visium_FFPE_Mouse_Brain"
#DATANAME="Visium_FFPE_Mouse_Brain_IF"
#DATANAME="Visium_FFPE_Mouse_Kidney"

ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/${DATANAME}/zscore/"
awk 'FNR>1 || NR==1' ${ZDIR}*.zscore > ${ZDIR}${DATANAME}.zscore
