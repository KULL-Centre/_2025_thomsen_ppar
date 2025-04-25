#!/bin/bash
#SBATCH --job-name=md_grompp_PPARg
#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
echo "========= Job started  at `date` =========="
cd $SLURM_SUBMIT_DIR
source /comm/specialstacks/gromacs-volta/bin/modules.sh
module load gromacs-gcc-8.2.0-openmpi-4.0.3-cuda-10.1

gmx_mpi grompp -f ../md.mdp -p all_PRO_lambda1.06.top -c relax.gro -t relax.cpt -o prodrun.tpr -maxwarn 1 -v

echo "========= Job finished at `date` =========="

