#!/bin/sh

species=pig
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results

tissueList=(adipose lung cerebellum cortex hypothalamus liver muscle spleen)
tissueList=(muscle)
numList=(p348 p350)
for tissue in "${tissueList[@]}"
do
 for num in "${numList[@]}"
 do
# 计算训练集、测试集、验证集的基因数目
 total_genes=$(wc -l < ${dataPath}/${species}/${tissue}/${num}/${species}_${tissue}_${num}_last_exp_count.csv)
 train_genes=$((total_genes * 7 / 10))
 test_genes=$((total_genes * 1 / 10))
 valid_genes=$((total_genes * 2 / 10))
 TSS_region_bed=${dataPath}/${species}/${tissue}/${num}/${species}_${tissue}_${num}_TSS_region_last.bed
 resPath=${dataPath}/${species}/${tissue}/${num}
# 划分训练集
 awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${resPath}/Train_TSS_region.bed

# 去除训练集中的行，生成剩余行的bed文件
 awk 'NR==FNR{a[$0];next} !($0 in a)' $resPath/Train_TSS_region.bed $TSS_region_bed > $resPath/TSS_remaining.bed

# 划分测试集
 awk -v n=$test_genes 'BEGIN{srand()} {print rand() "\t" $0}' $resPath/TSS_remaining.bed | sort -n | cut -f2- | head -n $((test_genes)) > $resPath/Test_TSS_region.bed

# 去除测试集中的行，生成剩余行的bed文件
 awk 'NR==FNR{a[$0];next} !($0 in a)' $resPath/Test_TSS_region.bed $resPath/TSS_remaining.bed > $resPath/TSS_remaining2.bed

# 划分验证集
 awk -v n=$valid_genes 'BEGIN{srand()} {print rand() "\t" $0}' $resPath/TSS_remaining2.bed | sort -n | cut -f2- | head -n $((valid_genes)) > $resPath/Valid_TSS_region.bed
 done
done
