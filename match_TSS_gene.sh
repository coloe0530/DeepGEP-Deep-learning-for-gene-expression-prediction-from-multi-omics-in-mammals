#!/bin/sh

module load R/4.0.3-gcc-4.8.5

species=cattle
TSS=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/${species}/TSS/region_10000.bed
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/${species}
datapath=/BIGDATA2/scau_xlyuan_1/CJL/keti/data/histone/${species}

#tissueList=(adipose lung cerebellum cortex hypothalamus liver spleen)
tissueList=(muscle)
numList=(M08 M22)
for num in "${numList[@]}"
 do
 cd ${resPath}
 mkdir ${num}
 for tissue in "${tissueList[@]}"
 do
   cd ${resPath}/${num}
   mkdir ${tissue}
   cd ${datapath}/${tissue}
   Rscript /BIGDATA2/scau_xlyuan_1/CJL/keti/2022.11test/code/macth_expression_geneID.R ensembl_name_ensembl_name_string_uniq.txt gene_count_matrix.csv ${resPath}/${species}_${tissue}_count_match_res.txt
   if [ ${num} == 'M08']
    then
    cut -f 1,3 ${resPath}/${species}_${tissue}_count_match_res.txt > ${resPath}/${num}/${tissue}/${species}_${tissue}_${num}_exp_count.txt
   else
    cut -f 1,4 ${resPath}/${species}_${tissue}_count_match_res.txt > ${resPath}/${num}/${tissue}/${species}_${tissue}_${num}_exp_count.txt
   fi
   sed -i 's/[ ]/\t/g' ${resPath}/${num}/${tissue}/${species}_${tissue}_${num}_exp_count.txt
   TPM=${resPath}/${num}/${tissue}/${species}_${tissue}_${num}_exp_count.txt
   cd ${resPath}/${num}/${tissue}
   awk 'NR==FNR{a[$4]=$0;next}NR>FNR{if($1 in a )print a[$1],$0}' ${TSS} ${TPM} > ${species}_${tissue}_${num}_region_10000_count.bed
   sort -k1,1V -k2,2n ${species}_${tissue}_${num}_region_10000_count.bed > ${species}_${tissue}_${num}_region_10000_count.sort.bed
   awk '{print $1,$2,$3,$4}' ${species}_${tissue}_${num}_region_10000_count.sort.bed > ${species}_${tissue}_${num}_TSS_region_last.bed
   sed -i 's/[ ]/\t/g' ${species}_${tissue}_${num}_TSS_region_last.bed
   awk '{print $5,$6}' ${species}_${tissue}_${num}_region_10000_count.sort.bed > ${species}_${tissue}_${num}_last_exp_count.csv
   sed -i 's/[ ]/,/g' ${species}_${tissue}_${num}_last_exp_count.csv
  done
done


