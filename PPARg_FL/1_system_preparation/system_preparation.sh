#!/bin/bash
#


export PATH="/lindorffgrp-isilon/thomasen/software/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lindorffgrp-isilon/wyong/software/openmpi401/lib

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi
python=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/python3.7
python2=/lindorffgrp-isilon/thomasen/software/miniconda2/bin/python2.7

martinize=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/martinize2

wget http://cgmartini.nl/images/tools/insane/insane.py
insane=insane.py

minmdp=minimization.mdp
FF=martini3001
ffdir=/storage1/thomasen/software/force-fields/Martini/martini_v300

dssp=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/mkdssp

PWrescaling_script=PW_rescaling_martini3.py
EN_script=choose_rubber_bands_PPARg.py

salt=0.15

mkdir PPARg
cd PPARg
pdb=../Structures/PPARg_FL_modeller_IDRmoved.pdb

#Martinize
$python $martinize -f $pdb -o PRO_topol.top -x PRO_CG.pdb -ff $FF -elastic -ef 700.0 -el 0.5 -eu 1.2 -ea 0 -ep 0 -ff-dir $ffdir/martini_v3.0.0_proteins/force_fields/ -map-dir $ffdir/martini_v3.0.0_proteins/mappings/ -ss CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCSCCBCSSSCSBCCEEETTEEECHHHHHHHHHHHHHTCCCCCCCSCCCCCTTGGGGCHHHHHHHHHHHCCCTTTCCCSCCSSTTTTTTCSSSSTHHHHSCSSHHHHHHHHHHHHHHHHHHCSSCHHHHHHHHHTCSSSCCCEEECSHHHHHHHHHSSCCCCCCSSCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTTSCHHHHHHHHHHHHHHHHHHHHHHHCCSSEEEEGGGTEEEEHHHHHTSSGGGGGGTHHHHHHHHHHHHHCCCHHHHHHHHHHHHSCTTSTTCSCHHHHHHHHHHHHHHHHHHHHHHCSSCTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCSCCCCCHHHHHHHTTCC

#Put protein in box
$gmx editconf -f PRO_CG.pdb -o PRO_CG.gro -bt dodecahedron -d 0.1 <<EOF
1
EOF

#Solvate using insane.py
$python2 ../$insane -f PRO_CG.gro -o PRO_SOL_IONS.gro -pbc keep -salt ${salt} -sol W -center -p PRO_topol_SOL_IONS.top

#The next few blocks modify the toplogy file and molecule_0.itp file:

#Remove #include include martini.itp and substitute ion names in topology file
perl -pi -e's/#include "martini.itp"//g' PRO_topol_SOL_IONS.top
perl -pi -e's/NA\+/NA/g' PRO_topol_SOL_IONS.top
perl -pi -e's/CL-/CL/g' PRO_topol_SOL_IONS.top

#Rename molecule_0.itp to PRO.itp and rename "molecule_0" as "Protein" in PRO.itp file
mv molecule_0.itp PRO.itp
perl -pi -e's/molecule_0/Protein/g' PRO.itp

#Add "#include .itp" lines to PRO_topol_SOL_IONS.top
cat <<EOF > others.top
#include "$ffdir/martini_v3.0.0.itp"
#include "PRO.itp"
#include "$ffdir/martini_v3.0.0_ions_v1.itp"
#include "$ffdir/martini_v3.0.0_solvents_v1.itp"
EOF
cat others.top PRO_topol_SOL_IONS.top >a
mv a PRO_topol_SOL_IONS.top

#Run energy minimization
$gmx grompp -f ../$minmdp -p PRO_topol_SOL_IONS.top -c PRO_SOL_IONS.gro -o trash.tpr -pp all_PRO.top -maxwarn 3 -r PRO_SOL_IONS.gro
rm trash.tpr

$python ../$EN_script
$python ../$PWrescaling_script -i rubberbands_all_PRO.top -o all_PRO_lambda1.06.top -l 1.06

