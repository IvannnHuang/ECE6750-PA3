#!/bin/bash
#SBATCH --exclusive
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --job-name=my_job_name
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=1

./serial

