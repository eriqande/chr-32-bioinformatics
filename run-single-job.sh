#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --output=run-single-out-%j
#SBATCH --error=run-single-error-%j 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=yourname@yourmail.edu


./map-N-files-from-K.sh 128 1

