#!/usr/bin/env bash

#ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20190_s004_4_BrainMetastasis/zscore/"
#DATANAME="p20190_s004_4_BrainMetastasis"

#ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20218_s001_L1/zscore/"
#DATANAME="p20218_s001_L1"

ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20190_s003_3_BrainMetastasis/zscore/"
DATANAME="p20190_s003_3_BrainMetastasis"
awk 'FNR>1 || NR==1' ${ZDIR}*.zscore > ${ZDIR}${DATANAME}.zscore
