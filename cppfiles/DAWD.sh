#!/bin/bash
#SBATCH --qos=tiger-short
#SBATCH --nodes=1
#SBATCH -c 10
#SBATCH --mem=80G
#SBATCH --time=23:59:59

module load anaconda3/5.3.1
conda activate project
./hydro.exe