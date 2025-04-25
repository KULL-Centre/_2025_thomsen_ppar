#!/bin/bash

mkdir calc_Rg

for traj in PPARg_rep1 PPARg_rep2
do

cd $traj
cp ../Rg_gyrate.sh .
nohup sh Rg_gyrate.sh $traj &
cd ..

done
