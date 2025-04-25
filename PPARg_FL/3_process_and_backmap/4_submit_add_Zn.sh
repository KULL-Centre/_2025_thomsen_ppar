#!/bin/bash

for traj in PPARg_rep1 PPARg_rep2
do

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N PPARg_addZn_${traj}" >> temp
cat add_Zn.sh >> temp
mv temp temp.sh

qsub temp.sh -v traj=$traj
rm temp.sh

done
