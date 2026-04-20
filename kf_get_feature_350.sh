#!/bin/sh

num=p350
tissue=spleen
species=pig
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results

#lung data
#H3K4me3="SRR12697675.bam"
#H3K4me1="SRR12697567.bam"
#H3K27me3="SRR12697615.bam"
#H3K27ac="SRR12697299.bam"
#ATAC="SRR12697174.bam"
#muscle data
#H3K4me3="SRR12697681.bam"
#H3K4me1="SRR12697573.bam"
#H3K27me3="SRR12697621.bam"
#H3K27ac="SRR12697305.bam"
#ATAC="SRR12697178.bam"
#adipose data
#H3K4me3="SRR12697645.bam"
#H3K4me1="SRR12697537.bam"
#H3K27me3="SRR12697585.bam"
#H3K27ac="SRR12697269.bam"
#ATAC="SRR12697155.bam"
#cerebellum data
#H3K4me3="SRR12697651.bam"
#H3K4me1="SRR12697543.bam"
#H3K27me3="SRR12697591.bam"
#H3K27ac="SRR12697275.bam"
#ATAC="SRR12697158.bam"
#cortex data
#H3K4me3="SRR12697657.bam"
#H3K4me1="SRR12697549.bam"
#H3K27me3="SRR12697597.bam"
#H3K27ac="SRR12697281.bam"
#ATAC="SRR12697162.bam"
#hypothalamus data
#H3K4me3="SRR12697663.bam"
#H3K4me1="SRR12697555.bam"
#H3K27me3="SRR12697603.bam"
#H3K27ac="SRR12697287.bam"
#ATAC="SRR12697166.bam"
#liver data
#H3K4me3="SRR12697669.bam"
#H3K4me1="SRR12697561.bam"
#H3K27me3="SRR12697609.bam"
#H3K27ac="SRR12697293.bam"
#ATAC="SRR12697170.bam"
#spleen data
H3K4me3="SRR12697687.bam"
H3K4me1="SRR12697579.bam"
H3K27me3="SRR12697627.bam"
H3K27ac="SRR12697311.bam"
ATAC="SRR12697182.bam"

module load anaconda3/2020.07
export PATH=/BIGDATA2/scau_xlyuan_1/.conda/envs/gwt/bin:$PATH

dataType=(Test_TSS_region Train_TSS_region_00 Train_TSS_region_01 Train_TSS_region_02 Train_TSS_region_03 Train_TSS_region_04)
#/////////////////////////////////////////////////////#
#并发数
threadTask=6  # <<<--- 设置并行数 ！！！
#创建fifo管道
fifoFile="test_fifo"
rm -f ${fifoFile}
mkfifo ${fifoFile}
# 建立文件描述符关联
exec 9<> ${fifoFile}
rm -f ${fifoFile}
# 预先向管道写入数据
for ((i=0;i<${threadTask};i++))
do
    echo "" >&9
done
echo "wait all task finish,then exit!!!"
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
for data in "${dataType[@]}"
do
read -u9              # <<<--- 添加该行
{                     # <<<--- 添加该行
cd ${resPath}/${species}/${num}/${tissue}

bedtools makewindows -b ${data}.bed -n 100 -i srcwinnum > ${data}_100bin.bed #200 bin of every region

bedtools multicov \
-bams \
${dataPath}/${species}/${tissue}/${H3K4me3} \
${dataPath}/${species}/${tissue}/${H3K4me1} \
${dataPath}/${species}/${tissue}/${H3K27me3} \
${dataPath}/${species}/${tissue}/${H3K27ac} \
${dataPath}/${species}/${tissue}/${ATAC} \
-bed \
${data}_100bin.bed > ${data}_feature_reads.csv
awk '{print $4,$5,$6,$7,$8,$9}' ${data}_feature_reads.csv|tr " " "," > ${data}_features.csv
     echo "" >&9      # <<<--- 添加该行
} &                  # <<<--- 添加该行
done
wait                 # <<<--- 添加该行

#/////////////////////////////////////////////////////#
# 关闭管道
exec 9>&-
echo
echo "success"
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
