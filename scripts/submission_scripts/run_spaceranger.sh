#!/bin/bash
#
#SBATCH --job-name=spaceranger
#SBATCH --output=/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/job_output/spaceranger.%j.out
#SBATCH --error=/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/job_output/spaceranger.%j.err
#SBATCH --time=24:00:00
##SBATCH --qos=normal
#SBATCH -p owners,horence,quake
##SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=64G

date
export PATH="/home/groups/horence/applications/spaceranger-1.3.1/:$PATH"
#SAMPLE="p20190-s003_3_BrainMetastasis"
#PATIENT="pt15"
#SLIDE="V10A21-107"
#AREA="C1"



#SAMPLE="p20190-s004_4_BrainMetastasis"
#PATIENT="pt19"
#SLIDE="V10A21-106"
#AREA="D1"

#SAMPLE="p20218-s003_L3"
#PATIENT="pt24"
#SLIDE="V10A21-106"
#AREA="C1"

SAMPLE="p20218-s001_L1"
PATIENT="pt26"
SLIDE="V10A21-106"
AREA="A1"

#SAMPLE="p20218-s002_L2"
#PATIENT="pt27"
#SLIDE="V10A21-106"
#AREA="B1"

#SAMPLE="p20218-s004_L4"
#PATIENT="pt16"
#SLIDE="V10A21-107"
#AREA="D1"

spaceranger count --id=${SAMPLE} \
	--transcriptome=/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/refdata-gex-GRCh38-2020-A \
	--fastqs=/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases \
	--sample=${SAMPLE} \
	--image=/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases/all_spatial_images/${PATIENT}.tif  \
	--slide=${SLIDE} \
	--area=${AREA} \
	--localmem=64
date
