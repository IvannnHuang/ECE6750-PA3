#!/bin/bash
#SBATCH --exclusive
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --job-name=pa3_spgemm
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=1

echo "MODE=$MODE IMPL=$IMPL matrix=$matrix p=$p"

make clean

if [ "$MODE" = "parallel" ]; then
    make IMPL="$IMPL" test="$test" pspgemm || exit 1
    ./pspgemm "$matrix" "$p" | awk 'NR<=20{print} NR==21{print "... output truncated after 20 lines ..."}'
else
    make test="$test" serial || exit 1
    ./serial "$matrix" | awk 'NR<=20{print} NR==21{print "... output truncated after 20 lines ..."}'
fi