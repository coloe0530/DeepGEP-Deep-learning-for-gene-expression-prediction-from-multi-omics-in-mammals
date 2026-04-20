#!/bin/sh

species=pig
num=p348
tissueList=(muscle)
#tissueList=(adipose lung cortex hypothalamus liver muscle spleen cerebellum)
#tissueList=(adipose lung cortex hypothalamus liver muscle spleen)

for tissue in "${tissueList[@]}"
do
yhbatch -p rhenv -J cjl_tpm /BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/run_kf.sh \
		-b 100 \
		-f 5 \
		-n ${num} \
		-t ${tissue} \
		-p ${species}
done
