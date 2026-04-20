#!/bin/sh

num=M08
tissue=spleen
species=cattle
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
log_norm=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/log_norm.py
#lung data
#H3K4me3="SRR12697672.bam"
#H3K4me1="SRR12697564.bam"
#H3K27me3="SRR12697612.bam"
#H3K27ac="SRR12697296.bam"
#ATAC="SRR12697171.bam"
#muscle data
#H3K4me3="SRR12697678.bam"
#H3K4me1="SRR12697570.bam"
#H3K27me3="SRR12697618.bam"
#H3K27ac="SRR12697302.bam"
#ATAC="SRR12697175.bam"
#adipose data
#H3K4me3="SRR12697642.bam"
#H3K4me1="SRR12697534.bam"
#H3K27me3="SRR12697582.bam"
#H3K27ac="SRR12697266.bam"
#ATAC="SRR12697152.bam"
#cerebellum data
#H3K4me3="SRR12697650.bam"
#H3K4me1="SRR12697542.bam"
#H3K27me3="SRR12697590.bam"
#H3K27ac="SRR12697274.bam"
#ATAC=".bam"
#cortex data
#H3K4me3="SRR12697654.bam"
#H3K4me1="SRR12697546.bam"
#H3K27me3="SRR12697594.bam"
#H3K27ac="SRR12697278.bam"
#ATAC="SRR12697159.bam"
#hypothalamus data
#H3K4me3="SRR12697660.bam"
#H3K4me1="SRR12697552.bam"
#H3K27me3="SRR12697600.bam"
#H3K27ac="SRR12697284.bam"
#ATAC="SRR12697163.bam"
#liver data
#H3K4me3="SRR12697666.bam"
#H3K4me1="SRR12697558.bam"
#H3K27me3="SRR12697606.bam"
#H3K27ac="SRR12697290.bam"
#ATAC="SRR12697167.bam"
#spleen data
H3K4me3="SRR12697684.bam"
H3K4me1="SRR12697576.bam"
H3K27me3="SRR12697624.bam"
H3K27ac="SRR12697308.bam"
ATAC="SRR12697179.bam"

module load bedtools2/2.26.0-gcc-4.8.5
module load anaconda3/2020.07

dataType=(Train Test Valid)
#/////////////////////////////////////////////////////#
#并发数
threadTask=3  # <<<--- 设置并行数 ！！！
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
cd ${resPath}/${species}/${tissue}/${num}

bedtools makewindows -b ${data}_TSS_region.bed -n 200 -i srcwinnum > ${data}_region_200bin.bed #200 bin of every region

bedtools multicov \
-bams \
${dataPath}/${species}/${tissue}/${H3K4me3} \
${dataPath}/${species}/${tissue}/${H3K4me1} \
${dataPath}/${species}/${tissue}/${H3K27me3} \
${dataPath}/${species}/${tissue}/${H3K27ac} \
${dataPath}/${species}/${tissue}/${ATAC} \
-bed \
${data}_region_200bin.bed > ${data}_feature_reads.csv
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

