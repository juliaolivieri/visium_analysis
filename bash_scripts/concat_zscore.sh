#!/usr/bin/env bash

ZDIR="/oak/stanford/groups/horence/JuliaO/nf-readzs/visium/p20190_s004_4_BrainMetastasis/zscore/"
DATANAME="p20190_s004_4_BrainMetastasis"
awk 'FNR>1 || NR==1' ${ZDIR}*.zscore > ${ZDIR}${DATANAME}.zscore
