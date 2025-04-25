#!/bin/bash

threads=10
pepsiscript=run_pepsi_constantparams_parallel.py

for traj in PPARg_rep1 PPARg_rep2
do

dir=${traj}/SAXS
#rm -r $dir
mkdir $dir
cp $pepsiscript $dir
cd $dir

max_job=$((${threads}-1))

for k in $(seq 0 $max_job)
do

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N SAXS_${traj}" >> temp
cat ../../run_pepsi_parallel.sh >> temp
mv temp run_pepsi_parallel.sh

qsub run_pepsi_parallel.sh -v start_frame=$k,skip_frame=$threads

done

cd ../..

done

