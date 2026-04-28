#!/usr/bin/env zsh
#SBATCH -p compphys2026
#SBATCH --job-name=FEA160
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=8:00:00
#SBATCH --output=FEA160.out
#SBATCH --error=FEA160.err

module load cmake/3.27.9
# module load valgrind/3.25.1

# Load valgrind on chtc spark
export VALGRIND_DIR=/home/grathod/lib/valgrind-3.26.0-install
export PATH=$VALGRIND_DIR/bin:$PATH
export LD_LIBRARY_PATH=$VALGRIND_DIR/lib:$LD_LIBRARY_PATH

# Define your thread array
THREADS=(2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48)

# 1. Build ONCE before the loop starts
export NEL_X1=160
export NEL_X2=48
export NEL_X3=48
mkdir -p "build/${NEL_X1}"
cmake -B "build/${NEL_X1}" -S .
cmake --build "build/${NEL_X1}" --parallel 24

# 2. Now run the experiment loop
for t in "${THREADS[@]}"; do
    export OMP_NUM_THREADS=$t
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    # Execute the binary directly instead of calling run_experiment.sh repeatedly
    ./build/${NEL_X1}/main > "experiments/experiment${NEL_X1}/FEA${NEL_X1}_${t}.out" 2>&1
done