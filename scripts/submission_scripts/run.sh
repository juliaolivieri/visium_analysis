#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=test.%j.out
#SBATCH --error=test.%j.err
#SBATCH --time=1:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=20G

#source ~/.bashrc
#conda deactivate
#conda activate readzs_env
#date

DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
WINDOW="chr17_5159_minus"

a="Rscript ../make_plotter_files.R ../output/prep_plot/pre_plotter/${DATANAME}_${WINDOW}.tsv 5000 quant 20 ../output/prep_plot ${DATANAME} /oak/stanford/groups/horence/JuliaO/nf-readzs/visium/${DATANAME} ../output/peak_plot/"
echo $a
eval $a

