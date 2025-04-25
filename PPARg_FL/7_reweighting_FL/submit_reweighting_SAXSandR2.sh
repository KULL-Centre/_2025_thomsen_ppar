#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

for R2_err in $(seq 0.00 0.05 0.50) #$(seq 0.0 0.01 0.30)
do
	echo $R2_err
	rm -r R2_err_fmod_${R2_err}
	mkdir R2_err_fmod_${R2_err}
	cd R2_err_fmod_${R2_err}
	nohup $python ../reweighting_SAXSandR2.py ${R2_err} &
	cd ..
done
