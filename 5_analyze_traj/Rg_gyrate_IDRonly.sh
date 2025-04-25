#!/bin/bash
gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

traj=$1

#$gmx grompp -f ../md.mdp -p all_PRO_lambda1.06.top -c relax.gro -o prodrun_BINF.tpr -maxwarn 1 -v

$gmx make_ndx -f prodrun_BINF.tpr -o IDR.ndx <<EOF
ri 1-137
q
EOF

$gmx trjconv -s prodrun_BINF.tpr -f prodrun_nopbc.xtc -o prodrun_IDRonly.xtc -n IDR.ndx <<EOF
r_1-137
EOF

$gmx editconf -f PRO_CG.gro -o PRO_CG_IDRonly.gro -n IDR.ndx <<EOF
r_1-137
EOF

$gmx gyrate -f prodrun_nopbc.xtc -s prodrun_BINF.tpr -o ../calc_Rg/Rg_gyrate_IDRonly_${traj}.xvg -n IDR.ndx <<EOF
r_1-137
EOF
