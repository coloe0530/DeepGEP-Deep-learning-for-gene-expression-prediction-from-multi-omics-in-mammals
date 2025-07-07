#!/bin/sh

#####################用户输入部分###################################################
GTF=$dataPath/GTF.gtf
#GTF路径
TSS_path=$dataPath/TSS
Fasta=$dataPath/dna_sm.toplevel.fa
Fai=$dataPath/dna_sm.toplevel.fa.fai
#TSS文件夹路径，参考基因组Fasta路径，Fai路径
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
#可选参数 length：选中TSS区域后，向两端各延伸多少bp
length = 5000

subread_bin = /BIGDATA2/scau_xlyuan_1/software/subread-2.0.6-Linux-x86_64/bin
conda_env = /BIGDATA2/scau_xlyuan_1/.conda/envs/gwt/bin
H3K4me3="SRR12697686.bam"
H3K4me1="SRR12697578.bam"
H3K27me3="SRR12697626.bam"
H3K27ac="SRR12697310.bam"
ATAC="SRR12697181.bam"

##############第一部分##########1.1获取TSS数据########################################

awk '$3=="gene" {print}'  $GTF > $TSS_path/gene.gtf
grep "protein_coding" $TSS_path/gene.gtf|less -S > $TSS_path/gene.pcg.gtf

awk '$7=="+" {print}' $TSS_path/gene.pcg.gtf > $TSS_path/gene.pcg_1.gtf
awk '$7=="-" {print}' $TSS_path/gene.pcg.gtf > $TSS_path/gene.pcg_2.gtf

#筛选PCG的TSS 加geneid（后面用来匹配）
awk '{print $1,$4,$10}' $TSS_path/gene.pcg_1.gtf > $TSS_path/pcgTSS_1.bed
awk '{print $1,$5,$10}' $TSS_path/gene.pcg_2.gtf > $TSS_path/pcgTSS_2.bed
cat $TSS_path/pcgTSS_1.bed $TSS_path/pcgTSS_2.bed > $TSS_path/pcgTSS.bed
sed -i 's/[;"]//g' $TSS_path/pcgTSS.bed
#提取常染色体TSS
grep -v KI $TSS_path/pcgTSS.bed|grep -v GL|grep -v MT|grep -v X|grep -v Y > $TSS_path/TSS.grep.bed
#transcript's start-1 as TSS's start
awk '{print $1,$2-1,$2,$3}' $TSS_path/TSS.grep.bed > $TSS_path/positive_TSS.bed
sed -i 's/[ ]/\t/g' $TSS_path/positive_TSS.bed
awk '{if($2>=5000) print $0}' $TSS_path/positive_TSS.bed > $TSS_path/positive_TSS_1.bed 
#positive_TSS_1.bed为TSS文件中的基因ID
#获得chom_size
module load samtools/1.11-gcc-4.8.5
samtools faidx $dataPath/$Fasta
cut -f1,2 $Fai > $dataPath/chrom.sizes
#region 
module load bedtools2/2.26.0-gcc-4.8.5
bedtools slop -i $TSS_path/positive_TSS_1.bed -g $dataPath/chrom.sizes -l length -r length > $TSS_path/region.bed

#####################################################1.2 run_expression表达量#################################################

module load stringtie
export PATH=$subread_bin:$PATH

stringtie -e -B -p 23 -G ${GTF} -o $dataPath/stringtie.gtf -A $dataPath/SAR.tsv ${dataPath}/SAR.bam
#输入RNA表达量${dataPath}/SAR.bam
featureCounts -T 23 -p -t exon -g gene_id -a ${GTF} -o $dataPath/featureCounts.txt ${dataPath}/SAR.bam
awk -v OFS="," 'NR>1{print $1,$7}' $dataPath/featureCounts.txt > $dataPath/exp_count.csv
#$dataPath/count.csv为TSS文件中的基因的表达量

#####################################1.3match tss gene#############################

module load R/4.0.3-gcc-4.8.5

cd ${dataPath}
sed -i 's/[ ]/\t/g' ${dataPath}/exp_count.txt
#此文件为步骤1.2的输出文件

TPM= ${dataPath}/exp_count.txt
TSS= $TSS_path/region.bed
#此文件为步骤1.1的输出文件

cd ${dataPath}
awk 'NR==FNR{a[$4]=$0;next}NR>FNR{if($1 in a )print a[$1],$0}' ${TSS} ${TPM} > region_count.bed
sort -k1,1V -k2,2n region_count.bed > region_count.sort.bed
awk '{print $1,$2,$3,$4}' region_count.sort.bed > TSS_region_last.bed
sed -i 's/[ ]/\t/g' TSS_region_last.bed
awk '{print $5,$6}' region_count.sort.bed > last_exp_count.csv
sed -i 's/[ ]/,/g' last_exp_count.csv


#################################################第2部分#################2.1对基因数据进行训练集和测试集的划分。将训练集进一步划分为五份，用于五折交叉验证。################################

cd ${dataPath}
# 计算训练集、测试集、验证集的基因数目
total_genes=$(wc -l < ${dataPath}/last_exp_count.csv)
train_genes=$((total_genes * 9 / 10))
test_genes=$((total_genes * 1 / 10))
TSS_region_bed=${dataPath}/TSS_region_last.bed

#划分大训练集
awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${dataPath}/Train_TSS_region.bed

# 得到测试集
awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Train_TSS_region.bed $TSS_region_bed > $dataPath/Test_TSS_region.bed

#将训练集划分为五份
dataset_line=$(wc -l $dataPath/Train_TSS_region.bed | awk '{print $1}')
split_lines=$((dataset_line / 5))
split -l $split_lines -d ${dataPath}/Train_TSS_region.bed ${dataPath}/Train_TSS_region_ --additional-suffix=.bed
 
 if [ -f "${dataPath}/Train_TSS_region_05.bed" ]; then
  # 合并第五个文件和第六个文件
  cat ${dataPath}/Train_TSS_region_04.bed ${dataPath}/Train_TSS_region_05.bed >> ${dataPath}/Train_TSS_region_04.bed

  # 删除第六个文件
  rm ${dataPath}/Train_TSS_region_05.bed
fi

########################################################################第2.2部分#################get feature###############

module load anaconda3/2020.07
export PATH=$conda_env:$PATH

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
cd ${dataPath}

bedtools makewindows -b ${data}.bed -n 100 -i srcwinnum > ${data}_100bin.bed #200 bin of every region

bedtools multicov \
-bams \
${dataPath}/${H3K4me3} \
${dataPath}/${H3K4me1} \
${dataPath}/${H3K27me3} \
${dataPath}/${H3K27ac} \
${dataPath}/${ATAC} \
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

##################第3部分#################进行五折交叉验证，生成训练集、验证集和测试集文件###############划分完训练集和测试集后，将训练集进一步分为5份#################

# Set the number of folds and the prefix for the output directory names
n_folds=5
prefix="fold_"

# Create the output directories 建 fold0, fold1, fold2, fold3, fold4 这五个文件夹

cd ${dataPath}
 for i in $(seq 0 $((n_folds-1))); do
 mkdir -p "${prefix}${i}"
done


# Perform 5-fold cross-validation

for i in $(seq 0 $((n_folds-1))); do
  # Rename the i-th data file as the validation set
  validation_file="${dataPath}/Train_TSS_region_0${i}_features.csv"
  cp "${validation_file}" "${dataPath}/${prefix}${i}/validation.csv"

  # Concatenate the remaining data files as the training set
  training_files=$(ls ${dataPath}/Train_TSS_region_0[0-4]_features.csv | grep -v "${validation_file}")
  cat ${training_files} > "${dataPath}/${prefix}${i}/training.csv"
  cp ${dataPath}/Test_TSS_region_features.csv ${dataPath}/${prefix}${i}/test.csv
  cp ${dataPath}/last_exp_count.csv ${dataPath}/${prefix}${i}/last_exp_count.csv
  cd ${dataPath}/${prefix}${i}
  cat training.csv validation.csv test.csv > all_gene_feature.csv
done

##################第4部分##########################先划分训练集，再从剩余数据中划分测试集，最后从剩余数据中划分验证集#######################

# 计算训练集、测试集、验证集的基因数目
total_genes=$(wc -l < ${dataPath}/last_exp_count.csv)
train_genes=$((total_genes * 7 / 10))
test_genes=$((total_genes * 1 / 10))
valid_genes=$((total_genes * 2 / 10))
TSS_region_bed=${dataPath}/TSS_region_last.bed

# 划分训练集
awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${dataPath}/Train_TSS_region.bed

# 去除训练集中的行，生成剩余行的bed文件
 awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Train_TSS_region.bed $TSS_region_bed > $dataPath/TSS_remaining.bed

# 划分测试集
 awk -v n=$test_genes 'BEGIN{srand()} {print rand() "\t" $0}' $dataPath/TSS_remaining.bed | sort -n | cut -f2- | head -n $((test_genes)) > $dataPath/Test_TSS_region.bed

# 去除测试集中的行，生成剩余行的bed文件
 awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Test_TSS_region.bed $dataPath/TSS_remaining.bed > $dataPath/TSS_remaining2.bed

# 划分验证集
 awk -v n=$valid_genes 'BEGIN{srand()} {print rand() "\t" $0}' $dataPath/TSS_remaining2.bed | sort -n | cut -f2- | head -n $((valid_genes)) > $dataPath/Valid_TSS_region.bed
