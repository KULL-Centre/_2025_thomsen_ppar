#!/bin/bash

python=/projects/prism/people/hzr104/envs/pulchra/bin/python3.12

mkdir backmapping_SAXS
cd backmapping_SAXS

$python ../backmap.py ../PPARa_AB/traj.xtc ../PPARa_AB/top.pdb traj_backmapped.xtc top_backmapped.pdb top_backmapped.gro

cd ..
