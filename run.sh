#!/bin/bash
#SBATCH --exclusive
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --job-name=pa3_spgemm
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=1

echo "MODE=$MODE IMPL=$IMPL matrix=$matrix p=$p test=$test"

make clean

if [ "$MODE" = "parallel" ]; then
    if [ "$IMPL" = "pthreads" ]; then
        make IMPL="$IMPL" test="$test" pthreads
        ./pthread "$matrix" "$p"
    elif [ "$IMPL" = "omp" ]; then
        make IMPL="$IMPL" test="$test" omp
        ./omp "$matrix" "$p"
    elif [ "$IMPL" = "balance" ]; then
        make IMPL="$IMPL" test="$test" balance
        ./balance "$matrix" "$p"
    else
        echo "Unknown IMPL=$IMPL"
        exit 1
    fi
else
    make test="$test" serial
    ./serial "$matrix"
fi