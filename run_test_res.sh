#!/bin/sh
module load anaconda3/2020.07
species=pig
#tissueList=(adipose lung cerebellum cortex hypothalamus liver muscle spleen)
hmList=(model_hm1 model_hm2 model_hm3 model_hm4)
tissue=muscle
num=p348
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
for hm in "${hmList[@]}"
do
python /BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/test.res.py \
	--data_root ${dataPath}/HMexplore/${hm} \
	--species ${species} \
	--tissue ${tissue} \
	--num ${num}
done
