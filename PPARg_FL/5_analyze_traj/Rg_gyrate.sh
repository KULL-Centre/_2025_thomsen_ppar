#!/bin/bash
gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

traj=$1

$gmx grompp -f ../md.mdp -p all_PRO_lambda1.06.top -c relax.gro -o prodrun_BINF.tpr -maxwarn 1 -v

$gmx gyrate -f prodrun_nopbc.xtc -s prodrun_BINF.tpr -o ../calc_Rg/Rg_gyrate_${traj}.xvg <<EOF
Protein
EOF
