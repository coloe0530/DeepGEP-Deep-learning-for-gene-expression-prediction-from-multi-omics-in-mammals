#!/bin/sh

num=M22
tissue=adipose
species=cattle
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
:<<EOF
#lung data
H3K4me3="SRR12697673.bam"
H3K4me1="SRR12697565.bam"
H3K27me3="SRR12697613.bam"
H3K27ac="SRR12697297.bam"
ATAC="SRR12697172.bam"
:<<EOF
#muscle data
H3K4me3="SRR12697679.bam"
H3K4me1="SRR12697571.bam"
H3K27me3="SRR12697619.bam"
H3K27ac="SRR12697303.bam"
ATAC="SRR12697176.bam"
EOF
#adipose data
H3K4me3="SRR12697643.bam"
H3K4me1="SRR12697535.bam"
H3K27me3="SRR12697583.bam"
H3K27ac="SRR12697267.bam"
ATAC="SRR12697153.bam"
:<<EOF
#cerebellum data
#H3K4me3="SRR12697649.bam"
#H3K4me1="SRR12697541.bam"
#H3K27me3="SRR12697589.bam"
#H3K27ac="SRR12697273.bam"
#ATAC="SRR12697156.bam"
:<<EOF
#cortex data
H3K4me3="SRR12697655.bam"
H3K4me1="SRR12697547.bam"
H3K27me3="SRR12697595.bam"
H3K27ac="SRR12697279.bam"
ATAC="SRR12697160.bam"
:<<EOF
#hypothalamus data
H3K4me3="SRR12697661.bam"
H3K4me1="SRR12697553.bam"
H3K27me3="SRR12697601.bam"
H3K27ac="SRR12697285.bam"
ATAC="SRR12697164.bam"
:<<EOF
#liver data
H3K4me3="SRR12697667.bam"
H3K4me1="SRR12697559.bam"
H3K27me3="SRR12697607.bam"
H3K27ac="SRR12697291.bam"
ATAC="SRR12697168.bam"
:<<EOF
#spleen data
H3K4me3="SRR12697685.bam"
H3K4me1="SRR12697577.bam"
H3K27me3="SRR12697625.bam"
H3K27ac="SRR12697309.bam"
ATAC="SRR12697180.bam"
EOF

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
