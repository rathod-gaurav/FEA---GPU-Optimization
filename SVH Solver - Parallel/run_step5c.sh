#!/usr/bin/env zsh
#SBATCH -p compphys2026
#SBATCH --job-name=FEA_Parallel_step5c
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=2:00:00
#SBATCH --output=FEA5c.out
#SBATCH --error=FEA5c.err

module load cmake/3.27.9
# module load valgrind/3.25.1

# Load valgrind on chtc spark
export VALGRIND_DIR=/home/grathod/lib/valgrind-3.26.0-install
export PATH=$VALGRIND_DIR/bin:$PATH
export LD_LIBRARY_PATH=$VALGRIND_DIR/lib:$LD_LIBRARY_PATH

rm -rf build5c/
mkdir build5c
cmake -B build5c
cmake --build build5c

BINARY="./build5c/main"

if [ ! -f "$BINARY" ]; then
    echo "Error: $BINARY not found. Run cmake + cmake --build first."
    exit 1
fi

# ── Detect machine ─────────────────────────────────────────────────────────────
CPU_MODEL=$(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2 | xargs)

# if echo "$CPU_MODEL" | grep -q "EPYC 9254"; then
#     # ── euler: AMD EPYC 9254 — Zen 4, 48 physical cores, 2 NUMA nodes ─────────
#     MACHINE="euler"
#     # 48 physical cores (96 logical with HT).
#     # HT shares the FPU per physical core — no benefit for AVX-512 FP work.
#     OMP_NUM_THREADS=16
#     # close: pack threads onto the same socket before crossing to the other.
#     # Keeps threads near their data (NUMA node affinity).
#     OMP_PROC_BIND=close
#     # cores: bind each OpenMP thread to one physical core.
#     # Prevents the OS from migrating threads between cores mid-run,
#     # which would invalidate the first-touch NUMA placement done in the code.
#     OMP_PLACES=cores

if echo "$CPU_MODEL" | grep -q "i5-7500"; then
    # ── workstation: Intel i5-7500 — Kaby Lake, 4 cores, 1 NUMA node ──────────
    MACHINE="lab workstation"
    OMP_NUM_THREADS=4
    OMP_PROC_BIND=close
    OMP_PLACES=cores

elif echo "$CPU_MODEL" | grep -q "5625U"; then
    # ── laptop: AMD Ryzen 5 5625U — Zen 3, 6 physical cores, 1 NUMA node ──────
    # Running under Hyper-V: OMP_PROC_BIND/PLACES may be partially ignored
    # by the hypervisor's virtual CPU scheduler, but still worth setting.
    MACHINE="laptop"
    OMP_NUM_THREADS=6
    OMP_PROC_BIND=close
    OMP_PLACES=cores

else
    # ── unknown machine — safe fallback ────────────────────────────────────────
    MACHINE="unknown ($CPU_MODEL)"
    # Use half of logical CPUs as a conservative default.
    # OMP_NUM_THREADS=$(( $(nproc --all) / 2 ))
    OMP_NUM_THREADS=16 #16 outer threads with each spawning 3 inner threads for a total of 48 threads
    OMP_PROC_BIND=close
    OMP_PLACES=cores
    echo "Warning: unrecognized CPU. Defaulting to manually set OMP_NUM_THREADS=${OMP_NUM_THREADS}."
    echo "Set OMP_NUM_THREADS manually for best results."
fi

# ── Export and report ──────────────────────────────────────────────────────────
export OMP_NUM_THREADS
export OMP_PROC_BIND
export OMP_PLACES

echo "================================================"
echo "Machine          : ${MACHINE}"
echo "CPU              : ${CPU_MODEL}"
echo "OMP_NUM_THREADS  : ${OMP_NUM_THREADS}"
echo "OMP_PROC_BIND    : ${OMP_PROC_BIND}"
echo "OMP_PLACES       : ${OMP_PLACES}"
echo "Binary           : ${BINARY}"
echo "================================================"

# ── Launch ─────────────────────────────────────────────────────────────────────
# exec "$BINARY" "$@"

# ── Step 5c: Cache behaviour at full scale ─────────────────────────────────
OMP_NUM_THREADS=16 \
OMP_MAX_ACTIVE_LEVELS=2 \
  likwid-perfctr \
  -C 0-47 \
  -g L3 \
  -m \
  ./build5c/main 2>&1 | tee likwid_l3.log