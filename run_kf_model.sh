#!/bin/sh

while getopts "m:s:d:n:t:p:b:u:" opt; do
  case $opt in
  	m) modelName=$OPTARG   ;;   
	s) saveRoot=$OPTARG   ;;
	d) dataRoot=$OPTARG   ;;
	n) nFeatures=$OPTARG   ;;
	t) tissue=$OPTARG   ;;
        p) species=$OPTARG   ;;
	b) n_bins=$OPTARG   ;;
        u) num=$OPTARG   ;;
	*) echo 'error' >&2
       exit 1
  esac
done


scriptPath=/BIGDATA2/scau_xlyuan_1/CJL/keti/results/code
module load anaconda3/2020.07

/usr/bin/time -f "%E\t%S\t%U\t%P\t%M\t/kf_model.py\t${tissue}\t${species}" -o ${saveRoot}/computeResources.log -a \
python ${scriptPath}/kf_model_tpm.py \
	--model_name $modelName \
	--save_root $saveRoot \
	--data_root $dataRoot \
	--n_features $nFeatures \
	--tissue ${tissue} \
	--species ${species} \
	--n_bins ${n_bins} \
	--num ${num}




