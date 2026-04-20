#!/bin/sh

while getopts "f:n:t:p:b:u:" opt; do
  case $opt in
	n) num=$OPTARG   ;;
        f) nFeatures=$OPTARG ;;
        t) tissue=$OPTARG   ;;
        p) species=$OPTARG   ;;
        b) n_bins=$OPTARG  ;;
        *) echo 'error' >&2
       exit 1
  esac
done

module load anaconda3/2020.07
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
scriptPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code
prefix="fold_"
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

#for i in $(seq 0 1)
iList=(0)
for i in "${iList[@]}"
do
read -u9              # <<<--- 添加该行
{                     # <<<--- 添加该行

cd ${resPath}/${species}/${num}/${tissue}/${prefix}${i}
mkdir -p model_tpm
	python ${scriptPath}/kf_model_tpm.py \
		--model_name LSTM_tpm_${i} \
		--save_root ${resPath}/${species}/${num}/${tissue}/${prefix}${i}/model_tpm \
		--data_root ${resPath}/${species}/${num}/${tissue}/${prefix}${i} \
		--n_features ${nFeatures} \
		--tissue ${tissue} \
		--species ${species} \
                --n_bins ${n_bins} \
		--num ${num}
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


