#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

mkdir backmapping_SAXS
cd backmapping_SAXS

$python ../backmap.py ../PPARg_AB/traj.xtc ../PPARg_AB/top.pdb traj_backmapped.xtc top_backmapped.pdb top_backmapped.gro

cd ..
