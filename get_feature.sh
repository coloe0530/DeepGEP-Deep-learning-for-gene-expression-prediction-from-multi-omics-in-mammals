#!/bin/sh

num=p348
tissue=muscle
species=pig
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
log_norm=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/log_norm.py
#lung data
#H3K4me3="SRR12697674.bam"
#H3K4me1="SRR12697566.bam"
#H3K27me3="SRR12697614.bam"
#H3K27ac="SRR12697298.bam"
#ATAC="SRR12697173.bam"
#muscle data
H3K4me3="SRR12697680.bam"
H3K4me1="SRR12697572.bam"
H3K27me3="SRR12697620.bam"
H3K27ac="SRR12697304.bam"
ATAC="SRR12697177.bam"
#adipose data
#H3K4me3="SRR12697644.bam"
#H3K4me1="SRR12697536.bam"
#H3K27me3="SRR12697584.bam"
#H3K27ac="SRR12697268.bam"
#ATAC="SRR12697154.bam"
#cerebellum data
#H3K4me3="SRR12697650.bam"
#H3K4me1="SRR12697542.bam"
#H3K27me3="SRR12697590.bam"
#H3K27ac="SRR12697274.bam"
#ATAC="SRR12697157.bam"
#cortex data
#H3K4me3="SRR12697656.bam"
#H3K4me1="SRR12697548.bam"
#H3K27me3="SRR12697596.bam"
#H3K27ac="SRR12697280.bam"
#ATAC="SRR12697161.bam"
#hypothalamus data
#H3K4me3="SRR12697662.bam"
#H3K4me1="SRR12697554.bam"
#H3K27me3="SRR12697602.bam"
#H3K27ac="SRR12697286.bam"
#ATAC="SRR12697165.bam"
#liver data
#H3K4me3="SRR12697668.bam"
#H3K4me1="SRR12697560.bam"
#H3K27me3="SRR12697608.bam"
#H3K27ac="SRR12697292.bam"
#ATAC="SRR12697169.bam"
#spleen data
#H3K4me3="SRR12697686.bam"
#H3K4me1="SRR12697578.bam"
#H3K27me3="SRR12697626.bam"
#H3K27ac="SRR12697310.bam"
#ATAC="SRR12697181.bam"

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

python ${log_norm} \
	--data_root ${resPath}/${species}/${tissue}/${num} \
	--species ${species}
	--tissue ${tissue}
	--num ${num}
