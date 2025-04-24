#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

cd backmapping_SAXS

$python ../run_pepsi.py ../../SAXS_expt/PPARg_ABonly_SAXSexpt.dat traj_backmapped.xtc top_backmapped.pdb calc_SAXS.dat

cd ..
