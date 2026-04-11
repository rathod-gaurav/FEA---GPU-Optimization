#!/usr/bin/env zsh
#SBATCH -p instruction
#SBATCH --job-name=FEA_Parallel
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=00:10:00
#SBATCH --output=FEA.out
#SBATCH --error=FEA.err

module load cmake/4.1.0

rm -rf build/
mkdir build
cmake -B build
cmake --build build

BINARY="./build/main"

if [ "$(nproc --all)" -gt 8 ]; then # If the system has more than 8 logical cores, the script assumes it is on the 20-core server
    # cpu2: 40 logical CPUs detected — use physical cores only, NUMA interleave
    echo "Detected cpu — Euler : running with 20 threads, NUMA interleave"
    export OMP_NUM_THREADS=20 #physical cores only, not HT siblings
    export OMP_PROC_BIND=close   # pack threads onto same socket first
    export OMP_PLACES=cores      # bind to physical cores, not HT siblings
    numactl --interleave=all $BINARY "$@"
    # --interleave=all: distribute memory pages round-robin across node0
    # and node1 at allocation time. Combined with the first-touch changes
    # in the code, this ensures mesh data is spread evenly rather than
    # all landing on node0.
else
    # cpu1: 4 cores — straightforward
    echo "Detected cpu - Lab — running with 4 threads"
    export OMP_NUM_THREADS=4
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    $BINARY "$@"
fi

