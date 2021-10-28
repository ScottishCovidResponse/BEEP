#!/bin/bash

#SBATCH --job-name="BEEPmbp"
#SBATCH --export=ALL

#SBATCH --partition=medium
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32


# Usage: sbatch submit.sh --inputfile "examples/EX1.toml"  --seed 10  --ngeneration 10 --nparticle 128 --dir "Output_EX1"

# dos2unix test.sh


inputfile="examples/infRegiontest.toml"
seed=0
ngeneration=10
nparticle=128
dir="Output"

while [ $# -gt 0 ]; do
    if [ $1 == "--inputfile" ]; then
        inputfile=$2
        shift 2
    elif [ $1 == "--seed" ]; then
        seed=$2
        shift 2
    elif [ $1 == "--ngeneration" ]; then
        ngeneration=$2
				shift 2
    elif [ $1 == "--nparticle" ]; then
        nparticle=$2
        shift 2
		elif [ $1 == "--dir" ]; then
        dir=$2
        shift 2
    else
        echo "Unrecognised arguments: $*" >&2
        exit 1
    fi
done


echo "Input file: $inputfile    seed: $seed    nproc: $nprocs    ngeneration: $ngeneration   nparticle: $nparticle    output: $dir"

echo "CPUs available: $SLURM_CPUS_PER_TASK"

echo "Starting job on $HOSTNAME"

mpirun -n 32 ./beepmbp inputfile=$inputfile mode="abcmbp" seed=$seed ngeneration=$ngeneration nparticle=$nparticle  outputdir=$dir
 
echo "Job finished"
