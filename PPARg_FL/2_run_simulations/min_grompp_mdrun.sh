#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:v100:1
echo "========= Job started  at `date` =========="
cd $SLURM_SUBMIT_DIR
source /comm/specialstacks/gromacs-volta/bin/modules.sh
module load gromacs-gcc-8.2.0-openmpi-4.0.3-cuda-10.1

export GMX_MAXCONSTRWARN=-1

gmx_mpi grompp -f ../minimization.mdp -c PRO_SOL_IONS.gro -p all_PRO_lambda1.06.top -o min.tpr -maxwarn 1
gmx_mpi mdrun -s min.tpr -deffnm min -ntomp 18 -maxh 23.9 -v


echo "========= Job finished at `date` =========="

