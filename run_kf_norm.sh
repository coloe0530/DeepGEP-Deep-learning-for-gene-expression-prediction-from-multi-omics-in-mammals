#!/bin/sh

species=cattle
num=M22
tissueList=(adipose lung muscle cortex hypothalamus liver muscle spleen)
#tissueList=(muscle)

for tissue in "${tissueList[@]}"
do
 yhbatch -N 1 -p rhenv -J cjl_2${tissue} /BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/kf_norm.sh -t ${tissue} -n ${num} -s ${species}
done

