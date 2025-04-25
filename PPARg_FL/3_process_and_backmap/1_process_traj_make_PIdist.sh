#!/bin/bash
#SBATCH --job-name=process_PI_PPARg
#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
echo "========= Job started  at `date` =========="
cd $SLURM_SUBMIT_DIR
source /comm/specialstacks/gromacs-volta/bin/modules.sh
module load gromacs-gcc-8.2.0-openmpi-4.0.3-cuda-10.1

gmx_mpi trjconv -s prodrun.tpr -f prodrun.xtc -o prodrun_nopbc.xtc -pbc mol -center <<EOF
1
1
EOF

gmx_mpi mindist -f prodrun_nopbc.xtc -s prodrun.tpr -od pi_mindist.xvg -tu us -pi <<EOF
1
EOF

echo "========= Job finished at `date` =========="


