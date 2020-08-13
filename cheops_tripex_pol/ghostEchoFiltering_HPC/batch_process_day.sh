#! /bin/bash -l
#SBATCH --tasks-per-node=1
#SBATCH --nodes=1
#SBATCH --account=AG-Kneifel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lterzi@uni-koeln.de
#SBATCH --time=48:00:00
#SBATCH --mem=45gb

#srun /scratch2/lterzi/tripex_pol/ghostEchoFiltering_HPC/regridData.py
srun /scratch2/lterzi/tripex_pol/ghostEchoFiltering_HPC/ghostEchoFilteringLeonie.py ${SLURM_ARRAY_TASK_ID} 
