#!/bin/bash

for rep in 1 2
do

cd PPARg_rep${rep}

cp ../prodrun_grompp.sh .
sbatch prodrun_grompp.sh
cd ..

done
