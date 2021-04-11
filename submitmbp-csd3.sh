#!/bin/bash

# Usage: ./submit-csd3 --groups 65536 --seed 0 --samples 10 --nprocs 32 --nnodes 1 --walltime 1:00:00 --dir /path/to/run/dir


# Usage: ./submitmbp-csd3.sh --inputfile  "examples/infMSOA.toml" --seed 10 --cputime 440 --nprocs 64 --nnodes 2 --walltime 8:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/Output1

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/infMSOAtest.toml" --seed 10 --cputime 220 --nprocs 128 --nnodes 4 --walltime 4:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/Output5

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/infMSOAtest.toml" --seed 10 --cputime 110 --nprocs 32 --nnodes 1 --walltime 2:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/OutputVV0

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/infRegiontest.toml"  --seed 10  --cputime 1000 --nprocs 32 --nnodes 1 --walltime 2:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/OutputVVV0 

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/inf_UK_deaths.toml"  --seed 10  --cputime 1000 --nprocs 128 --nnodes 4 --walltime 4:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/OutputUK20

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/inf_UK_new.toml"  --seed 10  --cputime 1000 --nprocs 128 --nnodes 4 --walltime 4:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/OutputUK22

# Usage: ./submitmbp-csd3.sh --inputfile  "examples/inf_Scotland_firstwave.toml"  --seed 10  --cputime 1000 --nprocs 128 --nnodes 4 --walltime 1:00:00 --dir ~/rds/rds-dirac-dc003/dc-pool1/OutputUK41


# squeue -u dc-pool1
# scancel <jobid>

# mpirun -n 1 ./beepmbp inputfile= "examples/inf_UK.toml"  mode="sim" seed=10 outputdir=Output


# dos2unix submitmbp-csd3.sh

set -e
set -u

# Defaults
inputfile="examples/infRegiontest.toml"
seed=0
cputime=10.0
nprocs=32
nnodes=1
walltime=1:00:00

while [ $# -gt 0 ]; do
    if [ $1 == "--inputfile" ]; then
        inputfile=$2
        shift 2
    elif [ $1 == "--seed" ]; then
        seed=$2
        shift 2
    elif [ $1 == "--samples" ]; then
        samples=$2
        shift 2
	  elif [ $1 == "--cputime" ]; then
        cputime=$2
        shift 2
    elif [ $1 == "--nprocs" ]; then
        nprocs=$2
				shift 2
    elif [ $1 == "--nnodes" ]; then
        nnodes=$2
				shift 2
    elif [ $1 == "--walltime" ]; then
        walltime=$2
        shift 2
		elif [ $1 == "--dir" ]; then
        dir=$2
        shift 2
    else
        echo "Unrecognised arguments: $*" >&2
        exit 1
    fi
done

jobname=$(basename $dir)

echo $dir

mkdir $dir

cp -a ./beepmbp $dir

mkdir $dir/examples

#mkdir $dir/examples/Data_ScotlandMSOA
mkdir $dir/examples/Data_UKLSOA
#mkdir $dir/examples/Data_ScotlandRegion

cp -a $inputfile $dir/$inputfile 

cp -a ./examples/Data_UKLSOA/* $dir/examples/Data_UKLSOA/
#cp -a ./examples/Data_ScotlandRegion/* $dir/examples/Data_ScotlandRegion/
 
cat >$dir/submit.sh <<EOF
#!/bin/bash
#SBATCH --job-name $jobname
#SBATCH --account DIRAC-DC003-CPU
#SBATCH --nodes $nnodes
#SBATCH --ntasks $nprocs
#SBATCH --time $walltime
#SBATCH --mail-type ALL
#SBATCH --no-requeue
#SBATCH --partition skylake
#SBATCH --output log.txt

####################################################################
# Setup environment
####################################################################

. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4

####################################################################
# Run BEEPmbp
####################################################################

echo
echo "Running simulation"
echo "=================="
echo
#mpirun -n 1 ./beepmbp inputfile=$inputfile mode="sim" seed=$seed outputdir=Output

echo
echo "Running inference"
echo "================="

mpirun -n $nprocs ./beepmbp inputfile=$inputfile mode="abcmbp" seed=$seed cputime=$cputime ngeneration=161 nparticle=256 outputdir=Output
EOF
						
						
cd $dir
chmod u+x submit.sh
sbatch submit.sh
