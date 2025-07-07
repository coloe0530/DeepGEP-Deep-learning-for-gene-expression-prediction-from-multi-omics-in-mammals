#!/bin/sh

#####################�û����벿��###################################################
GTF=$dataPath/GTF.gtf
#GTF·��
TSS_path=$dataPath/TSS
Fasta=$dataPath/dna_sm.toplevel.fa
Fai=$dataPath/dna_sm.toplevel.fa.fai
#TSS�ļ���·�����ο�������Fasta·����Fai·��
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
#��ѡ���� length��ѡ��TSS����������˸��������bp
length = 5000

subread_bin = /BIGDATA2/scau_xlyuan_1/software/subread-2.0.6-Linux-x86_64/bin
conda_env = /BIGDATA2/scau_xlyuan_1/.conda/envs/gwt/bin
H3K4me3="SRR12697686.bam"
H3K4me1="SRR12697578.bam"
H3K27me3="SRR12697626.bam"
H3K27ac="SRR12697310.bam"
ATAC="SRR12697181.bam"

##############��һ����##########1.1��ȡTSS����########################################

awk '$3=="gene" {print}'  $GTF > $TSS_path/gene.gtf
grep "protein_coding" $TSS_path/gene.gtf|less -S > $TSS_path/gene.pcg.gtf

awk '$7=="+" {print}' $TSS_path/gene.pcg.gtf > $TSS_path/gene.pcg_1.gtf
awk '$7=="-" {print}' $TSS_path/gene.pcg.gtf > $TSS_path/gene.pcg_2.gtf

#ɸѡPCG��TSS ��geneid����������ƥ�䣩
awk '{print $1,$4,$10}' $TSS_path/gene.pcg_1.gtf > $TSS_path/pcgTSS_1.bed
awk '{print $1,$5,$10}' $TSS_path/gene.pcg_2.gtf > $TSS_path/pcgTSS_2.bed
cat $TSS_path/pcgTSS_1.bed $TSS_path/pcgTSS_2.bed > $TSS_path/pcgTSS.bed
sed -i 's/[;"]//g' $TSS_path/pcgTSS.bed
#��ȡ��Ⱦɫ��TSS
grep -v KI $TSS_path/pcgTSS.bed|grep -v GL|grep -v MT|grep -v X|grep -v Y > $TSS_path/TSS.grep.bed
#transcript's start-1 as TSS's start
awk '{print $1,$2-1,$2,$3}' $TSS_path/TSS.grep.bed > $TSS_path/positive_TSS.bed
sed -i 's/[ ]/\t/g' $TSS_path/positive_TSS.bed
awk '{if($2>=5000) print $0}' $TSS_path/positive_TSS.bed > $TSS_path/positive_TSS_1.bed 
#positive_TSS_1.bedΪTSS�ļ��еĻ���ID
#���chom_size
module load samtools/1.11-gcc-4.8.5
samtools faidx $dataPath/$Fasta
cut -f1,2 $Fai > $dataPath/chrom.sizes
#region 
module load bedtools2/2.26.0-gcc-4.8.5
bedtools slop -i $TSS_path/positive_TSS_1.bed -g $dataPath/chrom.sizes -l length -r length > $TSS_path/region.bed

#####################################################1.2 run_expression�����#################################################

module load stringtie
export PATH=$subread_bin:$PATH

stringtie -e -B -p 23 -G ${GTF} -o $dataPath/stringtie.gtf -A $dataPath/SAR.tsv ${dataPath}/SAR.bam
#����RNA�����${dataPath}/SAR.bam
featureCounts -T 23 -p -t exon -g gene_id -a ${GTF} -o $dataPath/featureCounts.txt ${dataPath}/SAR.bam
awk -v OFS="," 'NR>1{print $1,$7}' $dataPath/featureCounts.txt > $dataPath/exp_count.csv
#$dataPath/count.csvΪTSS�ļ��еĻ���ı����

#####################################1.3match tss gene#############################

module load R/4.0.3-gcc-4.8.5

cd ${dataPath}
sed -i 's/[ ]/\t/g' ${dataPath}/exp_count.txt
#���ļ�Ϊ����1.2������ļ�

TPM= ${dataPath}/exp_count.txt
TSS= $TSS_path/region.bed
#���ļ�Ϊ����1.1������ļ�

cd ${dataPath}
awk 'NR==FNR{a[$4]=$0;next}NR>FNR{if($1 in a )print a[$1],$0}' ${TSS} ${TPM} > region_count.bed
sort -k1,1V -k2,2n region_count.bed > region_count.sort.bed
awk '{print $1,$2,$3,$4}' region_count.sort.bed > TSS_region_last.bed
sed -i 's/[ ]/\t/g' TSS_region_last.bed
awk '{print $5,$6}' region_count.sort.bed > last_exp_count.csv
sed -i 's/[ ]/,/g' last_exp_count.csv


#################################################��2����#################2.1�Ի������ݽ���ѵ�����Ͳ��Լ��Ļ��֡���ѵ������һ������Ϊ��ݣ��������۽�����֤��################################

cd ${dataPath}
# ����ѵ���������Լ�����֤���Ļ�����Ŀ
total_genes=$(wc -l < ${dataPath}/last_exp_count.csv)
train_genes=$((total_genes * 9 / 10))
test_genes=$((total_genes * 1 / 10))
TSS_region_bed=${dataPath}/TSS_region_last.bed

#���ִ�ѵ����
awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${dataPath}/Train_TSS_region.bed

# �õ����Լ�
awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Train_TSS_region.bed $TSS_region_bed > $dataPath/Test_TSS_region.bed

#��ѵ��������Ϊ���
dataset_line=$(wc -l $dataPath/Train_TSS_region.bed | awk '{print $1}')
split_lines=$((dataset_line / 5))
split -l $split_lines -d ${dataPath}/Train_TSS_region.bed ${dataPath}/Train_TSS_region_ --additional-suffix=.bed
 
 if [ -f "${dataPath}/Train_TSS_region_05.bed" ]; then
  # �ϲ�������ļ��͵������ļ�
  cat ${dataPath}/Train_TSS_region_04.bed ${dataPath}/Train_TSS_region_05.bed >> ${dataPath}/Train_TSS_region_04.bed

  # ɾ���������ļ�
  rm ${dataPath}/Train_TSS_region_05.bed
fi

########################################################################��2.2����#################get feature###############

module load anaconda3/2020.07
export PATH=$conda_env:$PATH

dataType=(Test_TSS_region Train_TSS_region_00 Train_TSS_region_01 Train_TSS_region_02 Train_TSS_region_03 Train_TSS_region_04)
#/////////////////////////////////////////////////////#
#������
threadTask=6  # <<<--- ���ò����� ������
#����fifo�ܵ�
fifoFile="test_fifo"
rm -f ${fifoFile}
mkfifo ${fifoFile}
# �����ļ�����������
exec 9<> ${fifoFile}
rm -f ${fifoFile}
# Ԥ����ܵ�д������
for ((i=0;i<${threadTask};i++))
do
    echo "" >&9
done
echo "wait all task finish,then exit!!!"
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
for data in "${dataType[@]}"
do
read -u9              # <<<--- ��Ӹ���
{                     # <<<--- ��Ӹ���
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
     echo "" >&9      # <<<--- ��Ӹ���
} &                  # <<<--- ��Ӹ���
done
wait                 # <<<--- ��Ӹ���

#/////////////////////////////////////////////////////#
# �رչܵ�
exec 9>&-
echo
echo "success" 

##################��3����#################�������۽�����֤������ѵ��������֤���Ͳ��Լ��ļ�###############������ѵ�����Ͳ��Լ��󣬽�ѵ������һ����Ϊ5��#################

# Set the number of folds and the prefix for the output directory names
n_folds=5
prefix="fold_"

# Create the output directories �� fold0, fold1, fold2, fold3, fold4 ������ļ���

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

##################��4����##########################�Ȼ���ѵ�������ٴ�ʣ�������л��ֲ��Լ�������ʣ�������л�����֤��#######################

# ����ѵ���������Լ�����֤���Ļ�����Ŀ
total_genes=$(wc -l < ${dataPath}/last_exp_count.csv)
train_genes=$((total_genes * 7 / 10))
test_genes=$((total_genes * 1 / 10))
valid_genes=$((total_genes * 2 / 10))
TSS_region_bed=${dataPath}/TSS_region_last.bed

# ����ѵ����
awk -v n=$train_genes 'BEGIN{srand()} {print rand() "\t" $0}' ${TSS_region_bed} | sort -n | cut -f2- | head -n $((train_genes)) > ${dataPath}/Train_TSS_region.bed

# ȥ��ѵ�����е��У�����ʣ���е�bed�ļ�
 awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Train_TSS_region.bed $TSS_region_bed > $dataPath/TSS_remaining.bed

# ���ֲ��Լ�
 awk -v n=$test_genes 'BEGIN{srand()} {print rand() "\t" $0}' $dataPath/TSS_remaining.bed | sort -n | cut -f2- | head -n $((test_genes)) > $dataPath/Test_TSS_region.bed

# ȥ�����Լ��е��У�����ʣ���е�bed�ļ�
 awk 'NR==FNR{a[$0];next} !($0 in a)' $dataPath/Test_TSS_region.bed $dataPath/TSS_remaining.bed > $dataPath/TSS_remaining2.bed

# ������֤��
 awk -v n=$valid_genes 'BEGIN{srand()} {print rand() "\t" $0}' $dataPath/TSS_remaining2.bed | sort -n | cut -f2- | head -n $((valid_genes)) > $dataPath/Valid_TSS_region.bed
