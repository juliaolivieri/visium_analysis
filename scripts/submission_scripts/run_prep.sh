#!/bin/bash
#
#SBATCH --job-name=prep
#SBATCH --output=../job_output/prep.%j.out
#SBATCH --error=../job_output/prep.%j.err
#SBATCH --time=1:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=60G


date

DATANAME="V1_Mouse_Brain_Sagittal_Anterior"
WINDOWFILE="windows.txt"
a="python -u ../prep_plot.py --dataname ${DATANAME} --window_file ${WINDOWFILE}"
echo $a
eval $a
 
date
