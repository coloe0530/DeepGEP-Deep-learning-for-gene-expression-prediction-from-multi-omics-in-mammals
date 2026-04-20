#!/bin/sh

num=M22
tissue=lung
species=cattle
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
log_norm=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/log_norm.py
#lung data
H3K4me3="SRR12697673.bam"
H3K4me1="SRR12697565.bam"
H3K27me3="SRR12697613.bam"
H3K27ac="SRR12697297.bam"
ATAC="SRR12697172.bam"
#muscle data
#H3K4me3="SRR12697679.bam"
#H3K4me1="SRR12697571.bam"
#H3K27me3="SRR12697619.bam"
#H3K27ac="SRR12697303.bam"
#ATAC="SRR12697176.bam"
#adipose data
#H3K4me3="SRR12697643.bam"
#H3K4me1="SRR12697535.bam"
#H3K27me3="SRR12697583.bam"
#H3K27ac="SRR12697267.bam"
#ATAC="SRR12697153.bam"
#cerebellum data
#H3K4me3="SRR12697649.bam"
#H3K4me1="SRR12697541.bam"
#H3K27me3="SRR12697589.bam"
#H3K27ac="SRR12697273.bam"
#ATAC="SRR12697156.bam"
#cortex data
#H3K4me3="SRR12697655.bam"
#H3K4me1="SRR12697547.bam"
#H3K27me3="SRR12697595.bam"
#H3K27ac="SRR12697279.bam"
#ATAC="SRR12697160.bam"
#hypothalamus data
#H3K4me3="SRR12697661.bam"
#H3K4me1="SRR12697553.bam"
#H3K27me3="SRR12697601.bam"
#H3K27ac="SRR12697285.bam"
#ATAC="SRR12697164.bam"
#liver data
#H3K4me3="SRR12697667.bam"
#H3K4me1="SRR12697559.bam"
#H3K27me3="SRR12697607.bam"
#H3K27ac="SRR12697291.bam"
#ATAC="SRR12697168.bam"
#spleen data
#H3K4me3="SRR12697685.bam"
#H3K4me1="SRR12697577.bam"
#H3K27me3="SRR12697625.bam"
#H3K27ac="SRR12697309.bam"
#ATAC="SRR12697180.bam"

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

