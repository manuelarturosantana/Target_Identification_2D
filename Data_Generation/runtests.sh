#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --qos=normal

#SBATCH --time=5-00:00:00                 # run time limit (DD-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks=1		  # This for multiple workers for parallel distributed jobs
#SBATCH --cpus-per-task=128                # This for multi threading. How many CPUS to use


module load matlab

# matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_constellation 
# matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_example  
# matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_twinjet
matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_quadjet
#matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_cavities    # ellipses and circles
# matlab -nodisplay -nosplash -nodesktop -nojvm -r data_gen_planes        # twinjet and quadjet geometries
