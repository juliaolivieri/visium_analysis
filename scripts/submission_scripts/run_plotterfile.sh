#!/bin/bash
#
#SBATCH --job-name=plotterfile
#SBATCH --output=../job_output/plotterfile.%j.out
#SBATCH --error=../job_output/plotterfile.%j.err
#SBATCH --time=4:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=20G

source ~/.bashrc
conda deactivate
conda activate jup_env
date

#DATANAME="V1_Mouse_Kidney"
#DATANAME="V1_Mouse_Brain_Sagittal_Posterior"
DATANAME="Visium_FFPE_Human_Breast_Cancer"
#DATANAME="Visium_FFPE_Mouse_Kidney"
WINDOWFILE="windows_${DATANAME}.txt"

#GFF="/oak/stanford/groups/horence/kaitlin/ref_files/gencode.vM26.annotation.gff3"
GFF="/oak/stanford/groups/horence/kaitlin/ref_files/gencode.v37.annotation.gff3 "

a="python -u ../prep_plot.py --dataname ${DATANAME} --window_file ${WINDOWFILE}"
echo $a
eval $a

conda deactivate
conda activate readzs_env
date


while read w; do
	echo "$w"
	a="Rscript ../make_plotter_files.R ../output/prep_plot/pre_plotter/${DATANAME}_${w}.tsv 5000 quant 20 ../output/prep_plot ${DATANAME} /oak/stanford/groups/horence/JuliaO/nf-readzs/visium/${DATANAME} ../output/peak_plot/"
	echo $a
	eval $a

	b="Rscript ../make_plot.R ../output/peak_plot/${DATANAME}_${w}.plotterFile quant ${GFF}"
	echo $b
	eval $b

done < $WINDOWFILE

#echo $a
#
#b="Rscript ../make_plot.R rank_1_chr17_5159_minus_Antkmt,Metrn_chi2__pval_.plotterFile quant /oak/stanford/groups/horence/kaitlin/ref_files/gencode.vM26.annotation.gff3"
date
