#!/bin/sh

while getopts "t:" opt; do
  case $opt in
    t) tissue=$OPTARG   ;;   # path+fa file name prefix,such as: path/susScr11
    *) echo 'error' >&2
       exit 1
  esac
done

numList=(p348 p350)
species=pig
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
log_norm=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/log_norm.py

module load bedtools2/2.26.0-gcc-4.8.5
module load anaconda3/2020.07
for num in "${numList[@]}"
do
python ${log_norm} \
	--data_root ${resPath}/${species}/${tissue}/${num} \
	--species ${species} \
	--tissue ${tissue} \
	--num ${num}
done
