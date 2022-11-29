#!/bin/bash
#
#SBATCH --job-name=png
#SBATCH --output=../job_output/png.%j.out
#SBATCH --error=../job_output/png.%j.err
#SBATCH --time=1:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=20G

DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
WINDOW="chr17_5159_minus"
PLOTTERFILE="../output/prep_plot/pre_plotter/${DATANAME}_${WINDOW}.tsv"
PATH="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/${DATANAME}"
a="Rscript ../make_plotter_files.R  ${PLOTTERFILE} 5000 quant 20  ../output/prep_plot  ${DATANAME} ${PATH} ../output/peak_plot/"
echo $a
eval $a

#Rscript ../make_plot.R \
#  rank_1_chr17_5159_minus_Antkmt,Metrn_chi2__pval_.plotterFile \
#  quant \
#  /oak/stanford/groups/horence/kaitlin/ref_files/gencode.vM26.annotation.gff3
