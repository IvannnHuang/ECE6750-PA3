#!/bin/bash
#SBATCH --exclusive
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --job-name=pa3_spgemm
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=1

echo "MODE=$MODE IMPL=$IMPL matrix=$matrix p=$p"

# Rebuild on compute node so -march=native targets the right CPU
make clean
if [ "$MODE" = "parallel" ]; then
    make IMPL="$IMPL" test="$test" pspgemm
    ./pspgemm "$matrix" "$p" | awk 'NR<=20{print} NR==21{print "... output truncated after 20 lines ..."}'
else
    make test="$test" serial
    ./serial "$matrix" | awk 'NR<=20{print} NR==21{print "... output truncated after 20 lines ..."}'
fi