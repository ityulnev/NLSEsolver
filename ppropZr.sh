#!/bin/bash
#SBATCH --partition=all                                                   ## this is the default partition.
#SBATCH -t 55:20:00                                                         ## default is 1h. The maximum is partition dependent, have a look at sview or scontrol for details.
#SBATCH --nodes=1                                                           ## number of nodes

#SBATCH --ntasks-per-node=1
##SBATCH --array 0-3
##SBATCH --constraint GPU
##SBATCH -n 50                                                           ## Number of threads. 
#SBATCH --output    output_figures/%j-%N.out                                 ## File to which STDOUT will be written
#SBATCH --error     output_figures/%j-%N.err                                 ## File to which STDERR will be written
#SBATCH --mail-type END                                                   ## Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user igor.tyulnev@desy.de                                       ## Email to which notifications will be sennet 
##SBATCH --dependency singleton
##slurm stop reads batch command after the first normal commands "export ....."


##load mpi model
export LD_PRELOAD=""
source  /etc/profile.d/modules.sh
##module load mpi/mpich-3.2-x86_64

echo "nlse solver"
echo "SLURM_JOB_ID           $SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID     $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID    $SLURM_ARRAY_TASK_ID"
echo "SLURM_ARRAY_TASK_COUNT $SLURM_ARRAY_TASK_COUNT"
echo "SLURM_ARRAY_TASK_MAX   $SLURM_ARRAY_TASK_MAX"
echo "SLURM_ARRAY_TASK_MIN   $SLURM_ARRAY_TASK_MIN"
echo "SLURM_JOB_NODELIST     $SLURM_JOB_NODELIST"

module load matlab
##srun  matlab -nodisplay -nosplash -r  "tpf2D_main($SLURM_ARRAY_TASK_ID)" 
srun  matlab -nodisplay -nosplash -r  "main_NLSEsolver3D" 
exit