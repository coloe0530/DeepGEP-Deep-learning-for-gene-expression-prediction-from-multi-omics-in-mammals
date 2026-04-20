#!/bin/sh

#tissueList=(adipose lung cerebellum cortex hypothalamus liver muscle spleen)
tissueList=(muscle)
#/////////////////////////////////////////////////////#
#并发数
threadTask=2  # <<<--- 设置并行数 ！！！
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
for tissue in "${tissueList[@]}"
do
read -u9              # <<<--- 添加该行
{                     # <<<--- 添加该行
 yhbatch -N 1 -p bigdata -J cjl_${tissue} /BIGDATA2/scau_xlyuan_1/CJL/keti/results/code/tissun_log_norm.sh -t ${tissue}
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
