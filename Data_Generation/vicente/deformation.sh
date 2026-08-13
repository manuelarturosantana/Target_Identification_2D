#!/bin/bash

#SBATCH --job-name=deformation
#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=1000G
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --array=1-26%3

#SBATCH --chdir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente
#SBATCH --output=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/deformation_%A_%a.out
#SBATCH --error=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/deformation_%A_%a.err

module load matlab

status_dir=/home/vhojas/Code/Target_Identification_2D/Data_Generation/vicente/logs/status/deformation
mkdir -p "$status_dir"

echo "Starting deformation config ${SLURM_ARRAY_TASK_ID}"

matlab -nodisplay -nosplash -nodesktop -nojvm -singleCompThread \
    -batch "data_gen_deformation"

exit_code=$?

if [ "$exit_code" -eq 0 ]; then
    echo "DONE config ${SLURM_ARRAY_TASK_ID}" \
        > "${status_dir}/config_${SLURM_ARRAY_TASK_ID}.done"
else
    echo "FAILED config ${SLURM_ARRAY_TASK_ID} with exit code ${exit_code}" \
        > "${status_dir}/config_${SLURM_ARRAY_TASK_ID}.failed"
fi

exit "$exit_code"