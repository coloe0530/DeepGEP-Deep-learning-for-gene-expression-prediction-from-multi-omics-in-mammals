#!/bin/sh

species=cattle
resPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results
dataPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/cattle-norm
numList=(M08 M22)
tissueList=(muscle adipose lung cortex hypothalamus liver spleen)
# Set the number of folds and the prefix for the output directory names
n_folds=5
prefix="fold_"

for num in "${numList[@]}"
 do
 for tissue in "${tissueList[@]}"
 do
# Create the output directories
 outPath=${resPath}/${species}/${num}/${tissue}
 cd ${outPath}
 for i in $(seq 0 $((n_folds-1))); do
 mkdir -p "${prefix}${i}"
done
done
done

# Perform 5-fold cross-validation
for num in "${numList[@]}"
 do
 for tissue in "${tissueList[@]}"
 do
 for i in $(seq 0 $((n_folds-1))); do
  # Rename the i-th data file as the validation set
  Path=${dataPath}/${num}/${tissue}
  outPath=${resPath}/${species}/${num}/${tissue}
  validation_file="${Path}/Train_TSS_region_0${i}_features.csv"
  cp "${validation_file}" "${outPath}/${prefix}${i}/validation.csv"

  # Concatenate the remaining data files as the training set
  training_files=$(ls ${Path}/Train_TSS_region_0[0-4]_features.csv | grep -v "${validation_file}")
  cat ${training_files} > "${outPath}/${prefix}${i}/training.csv"
  cp ${Path}/Test_TSS_region_features.csv ${outPath}/${prefix}${i}/test.csv
  cp ${Path}/${species}_${tissue}_${num}_last_exp_count.csv ${outPath}/${prefix}${i}/last_exp_count.csv
  cd ${outPath}/${prefix}${i}
  cat training.csv validation.csv test.csv > all_gene_feature.csv
done
done
done


