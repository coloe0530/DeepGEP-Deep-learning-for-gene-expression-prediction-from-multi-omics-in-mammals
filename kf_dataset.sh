#!/bin/sh

dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/
species=cattle
numList=(M08 M22)
#tissueList=(muscle adipose lung cerebellum cortex hypothalamus liver spleen)
tissueList=(muscle)
for num in "${numList[@]}"
 do
 for tissue in "${tissueList[@]}"
 do 
 cd ${dataPath}/${species}/${num}
 # 计算训练集、测试集、验证集的基因数目
 total_genes=$(wc -l < ${dataPath}/${species}/${num}/${tissue}/${species}_${tissue}_${num}_last_exp_count.csv)
 train_genes=$((total_genes * 9 / 10))
 test_genes=$((total_genes * 1 / 10))
 TSS_region_bed=${dataPath}/${species}/${num}/${tissue}/${species}_${tissue}_${num}_TSS_region_last.bed
 resPath=${dataPath}/${species}/${num}/${tissue}
 #划分大训练集
 awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${resPath}/Train_TSS_region.bed
 # 得到测试集
 awk 'NR==FNR{a[$0];next} !($0 in a)' $resPath/Train_TSS_region.bed $TSS_region_bed > $resPath/Test_TSS_region.bed
 #将训练集划分为五份
 dataset_line=$(wc -l ${resPath}/Train_TSS_region.bed | awk '{print $1}')
 split_lines=$((dataset_line / 5))
 split -l $split_lines -d ${resPath}/Train_TSS_region.bed ${resPath}/Train_TSS_region_ --additional-suffix=.bed
 
 if [ -f "${resPath}/Train_TSS_region_05.bed" ]; then
  # 合并第五个文件和第六个文件
  cat ${resPath}/Train_TSS_region_04.bed ${resPath}/Train_TSS_region_05.bed >> ${resPath}/Train_TSS_region_04.bed

  # 删除第六个文件
  rm ${resPath}/Train_TSS_region_05.bed
fi
done
done
