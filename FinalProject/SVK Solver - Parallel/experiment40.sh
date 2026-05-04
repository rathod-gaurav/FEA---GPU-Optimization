#!/usr/bin/env zsh
#SBATCH -p instruction
#SBATCH --job-name=FEA40
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:45:00
#SBATCH --output=FEA40.out
#SBATCH --error=FEA40.err

module load cmake
# module load valgrind/3.25.1

# Load valgrind on chtc spark
# export VALGRIND_DIR=/home/grathod/lib/valgrind-3.26.0-install
# export PATH=$VALGRIND_DIR/bin:$PATH
# export LD_LIBRARY_PATH=$VALGRIND_DIR/lib:$LD_LIBRARY_PATH

# Define your thread array
THREADS=(2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48)

# 1. Build ONCE before the loop starts
export NEL_X1=40
export NEL_X2=12
export NEL_X3=12
mkdir -p "build/${NEL_X1}"
cmake -B "build/${NEL_X1}" -S .
cmake --build "build/${NEL_X1}" --parallel 24

# 2. Now run the experiment loop
for t in "${THREADS[@]}"; do
    export OMP_NUM_THREADS=$t
    # Execute the binary directly instead of calling run_experiment.sh repeatedly
    ./build/${NEL_X1}/main > "experiments/experiment${NEL_X1}/FEA${NEL_X1}_${t}.out" 2>&1
done