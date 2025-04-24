#!/bin/bash

python=/projects/prism/people/hzr104/envs/md_standard/bin/python3.12

cd backmapping_SAXS

$python ../run_pepsi.py ../SAXS_expt/PPARa_AB.dat traj_backmapped.xtc top_backmapped.pdb calc_SAXS.dat

cd ..
