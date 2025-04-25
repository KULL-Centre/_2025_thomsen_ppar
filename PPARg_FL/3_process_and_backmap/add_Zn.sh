#PBS -l nodes=1:ppn=1:thinnode
#PBS -l walltime=480:00:00
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes


python=/home/projects/ku_10001/people/fretho/miniconda3/bin/python3.8

$python add_Zn.py $traj
