#!/bin/bash
#
#SBATCH --job-name=overview
#SBATCH --output=../job_output/overview.%j.out
#SBATCH --error=../job_output/overview.%j.err
#SBATCH --time=4:00:00
#SBATCH -p owners,horence,quake
#SBATCH --nodes=1
#SBATCH --mem=60G

date
source ~/.bashrc
conda activate jup_env

a="python ../make_overview_table.py"
echo $a
eval $a
date
