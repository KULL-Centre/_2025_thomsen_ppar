#!/bin/sh
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N trjcat
#PBS -l nodes=1:ppn=1:thinnode
#PBS -l walltime=200:00:00
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
# Load all required modules for the job
module load tools
module load cuda/toolkit/10.2.89 openmpi/gcc/64/1.10.2 gcc/9.3.0
gmx=/home/projects/ku_10001/apps/GMX20203/bin/gmx_mpi

backmapped_dir=${traj}/Backmapping
output_dir=${traj}
tmp_workingdir=${traj}/trj_cat_tmp
last_frame=40000

#rm -r ${tmp_workingdir}
mkdir ${tmp_workingdir}

for i in $(seq 0 ${last_frame})
do

$gmx trjconv -f ${backmapped_dir}/AA_frame${i}.pdb -o ${tmp_workingdir}/AA_frame${i}.xtc

echo -n " ${tmp_workingdir}/AA_frame${i}.xtc" >> ${tmp_workingdir}/frame_input_tmp.txt
echo $(( ${i}*1000 )) >> ${tmp_workingdir}/time_input_tmp.txt

done


$gmx trjcat -f `cat ${tmp_workingdir}/frame_input_tmp.txt` -o ${output_dir}/prodrun_AAbackmapped.xtc -cat -settime < ${tmp_workingdir}/time_input_tmp.txt

$gmx trjconv -f ${backmapped_dir}/AA_frame0.pdb -s ${backmapped_dir}/AA_frame0.pdb -o ${output_dir}/prodrun_AAbackmapped.gro <<EOF
1
EOF

rm -r ${tmp_workingdir}
