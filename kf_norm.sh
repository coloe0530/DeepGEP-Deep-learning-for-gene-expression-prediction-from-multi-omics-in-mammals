#!/bin/sh

while getopts "t:n:s:" opt; do
  case $opt in
    t) tissue=$OPTARG   ;;
    n) num=$OPTARG   ;;
    s) species=$OPTARG   ;;   # path+fa file name prefix,such as: path/susScr11
    *) echo 'error' >&2
       exit 1
  esac
done

prefix="fold_"
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/
log_norm=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/log_norm.py

module load anaconda3/2020.07
export PATH=/BIGDATA2/scau_xlyuan_1/.conda/envs/gwt/bin:$PATH

#/////////////////////////////////////////////////////#
#并发数
threadTask=5  # <<<--- 设置并行数 ！！！
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

for i in $(seq 0 4)
do
read -u9              # <<<--- 添加该行
{                     # <<<--- 添加该行
python ${log_norm} \
	--data_root ${resPath}/${species}/${num}/${tissue}/${prefix}${i} \
	--species ${species} \
	--tissue ${tissue} \
	--num ${num}
cd ${resPath}/${species}/${num}/${tissue}/${prefix}${i}
cat train_norm_feature.csv test_norm_feature.csv valid_norm_feature.csv > all_gene_norm_feature.csv
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

