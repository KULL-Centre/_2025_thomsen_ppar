#!/bin/bash

for traj in PPARg_rep1 PPARg_rep2
do

cd ${traj}

mkdir Backmapping
cd Backmapping

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N PPARg_${traj}_backmap" >> temp
cat ../../Backmap.sh >> temp
mv temp Backmap.sh

qsub Backmap.sh -v traj=$traj

cd ../..

done
