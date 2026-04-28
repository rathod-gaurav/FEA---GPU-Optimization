cmake_minimum_required(VERSION 3.16)
# 3.16+ recommended: better target-based dependency propagation,
# precompiled header support (useful if you add PCH later),
# and more reliable OpenMP detection on Linux.

project(MyProject LANGUAGES CXX)

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
# Extensions OFF: disables GCC-specific extensions (-std=gnu++17 → -std=c++17).
# Keeps behaviour portable and avoids subtle ABI differences.

# ── Build type default ────────────────────────────────────────────────────────
# If the user doesn't specify -DCMAKE_BUILD_TYPE, default to Release.
# Debug builds without this will have no -O3 and your code will be very slow.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "No build type specified — defaulting to Release")
endif()


# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(OpenMP REQUIRED)
find_package(CURL REQUIRED)
# REQUIRED: hard-fail if OpenMP is missing rather than silently building
# without parallelism. OpenMP::OpenMP_CXX is an imported target that
# automatically adds -fopenmp to compile AND link flags.

# Eigen: use an INTERFACE imported target instead of include_directories().
# This scopes the include path to targets that actually need it,
# and propagates correctly if other targets depend on yours.
add_library(Eigen3 INTERFACE IMPORTED)
set_target_properties(Eigen3 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES
    # "/mnt/c/Users/ratho/Downloads/eigen-5.0.0/eigen-5.0.0" #for laptop
    "/home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0" #for lab workstation
    # "/srv/home/grathod/eigen-5.0.0" #for euler
    # "/home/grathod/eigen-5.0.0" #for spark
)
# If you ever install Eigen system-wide or via vcpkg/conan you can replace
# the two lines above with: find_package(Eigen3 REQUIRED)
# and link against Eigen3::Eigen instead.

# ── Sources ───────────────────────────────────────────────────────────────────
file(GLOB SOURCES "src/*.cpp")
# Note: file(GLOB) doesn't re-run when you add new .cpp files — you'd need
# to re-run cmake manually. An explicit list is safer for larger projects,
# but GLOB is fine here.

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(main main.cpp ${SOURCES})


#############################################################
###below part varies depending on the target cpu architecture
#############################################################

# -march=native: detected at compile time for whichever machine runs cmake.
#   laptop      → znver3  (AVX2 + FMA, no AVX-512)
#   workstation → kabylake (AVX2 + FMA, no AVX-512)
#   euler       → znver4  (AVX2 + FMA + AVX-512, 8 doubles/op)
# Rule: compile on the machine you run on. Never copy binaries between machines.
#
# Override at configure time if needed:
#   cmake -DTARGET_ARCH=x86-64-v3    ← portable across all 3, loses AVX-512 on euler
#   cmake -DTARGET_ARCH=znver4       ← euler only, crashes on laptop/workstation
set(TARGET_ARCH "native" CACHE STRING
    "CPU target: native (default) | x86-64-v3 | znver3 | kabylake | znver4")
message(STATUS "TARGET_ARCH = ${TARGET_ARCH}")

# Detect at configure time whether this machine supports AVX-512,
# so we can add -mprefer-vector-width=512 only when relevant.
# This check is skipped when TARGET_ARCH is overridden to a non-native value.
if(TARGET_ARCH STREQUAL "native" OR TARGET_ARCH STREQUAL "znver4")
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS "-mavx512f")
    check_cxx_source_runs("
        #include <immintrin.h>
        int main() {
            __m512d a = _mm512_set1_pd(1.0);
            __m512d b = _mm512_set1_pd(2.0);
            __m512d c = _mm512_add_pd(a, b);
            return 0;
        }
    " HAS_AVX512)
    unset(CMAKE_REQUIRED_FLAGS)
else()
    set(HAS_AVX512 FALSE)
endif()

if(HAS_AVX512)
    message(STATUS "AVX-512 detected — enabling 512-bit vector width preference")
else()
    message(STATUS "AVX-512 not detected — using AVX2 (256-bit)")
endif()

# ── Compile options ───────────────────────────────────────────────────────────
# Use target_compile_options instead of CMAKE_CXX_FLAGS.
# CMAKE_CXX_FLAGS is global and bleeds into ALL targets (including test
# executables, third-party subdirectories, etc.).
# target_compile_options scopes flags to this target only.
target_compile_options(main PRIVATE
    # Optimisation flags (only for this target, only in Release)
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-march=native>
    # -march=native: enables AVX2, FMA, etc. for your hardware automatically.
    # WARNING: binaries built with -march=native are not portable to other
    # machines. Fine for a research/local build, not for distribution.

    # On euler (znver4) with native: tell GCC to prefer 512-bit SIMD ops
    # over falling back to 256-bit. No effect on laptop/workstation (no AVX-512).
    $<$<AND:$<CONFIG:Release>,$<BOOL:${HAS_AVX512}>>:-mprefer-vector-width=512>

    $<$<CONFIG:Release>:-ffast-math>
    # -ffast-math: allows reordering of FP operations, enables vectorisation
    # of reductions, assumes no NaN/Inf. Gives ~10-30% speedup on FP-heavy
    # code like yours. Side effect: breaks std::isnan() checks and strict
    # IEEE 754 compliance. Acceptable for FEM where you control inputs.

    $<$<CONFIG:Release>:-funroll-loops>
    # Unrolls small loops (like your 3x3 Eigen operations) so the CPU's
    # out-of-order engine has more instructions to pipeline simultaneously.

    # Debug flags
    $<$<CONFIG:Debug>:-O2>
    $<$<CONFIG:Debug>:-g>
    # $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    # AddressSanitizer + UBSan: catches memory errors and undefined behaviour
    # (out-of-bounds Eigen access, signed overflow, etc.) at runtime in Debug.
    # Pairs with the linker flag below.
)

# ── Compile definitions ───────────────────────────────────────────────────────
target_compile_definitions(main PRIVATE
    EIGEN_DONT_PARALLELIZE
    # Disables Eigen's internal OpenMP threading. Without this, Eigen may
    # spawn its own threads inside your OpenMP threads, causing
    # oversubscription on your 4-core i5-7500. Your explicit parallelism
    # (the #pragma omp parallel for in computeElement) is the only
    # parallelism you want. Applied at compile time — equivalent to putting
    # #define EIGEN_DONT_PARALLELIZE before every #include <Eigen/...>.

    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
    # Disables Eigen's internal assertion checks (bounds checking on
    # matrix dimensions, etc.) in Release builds. Removes overhead from
    # hot paths like block<3,3>() calls inside your quadrature loop.
    # Keep it ON in Debug so assertion failures point you to bugs.

    $<$<CONFIG:Release>:NDEBUG>
    # Disables standard assert() calls in Release. Also disables some
    # Eigen checks that are gated on NDEBUG rather than EIGEN_NO_DEBUG.
)

# ── Include paths ─────────────────────────────────────────────────────────────
target_include_directories(main PRIVATE
    "include"
    # Your own project headers. PRIVATE means this path is not propagated
    # to targets that link against main (not relevant for an executable,
    # but good habit for when you split into libraries).
)

# ── Link libraries ────────────────────────────────────────────────────────────
target_link_libraries(main PRIVATE
    CURL::libcurl
    OpenMP::OpenMP_CXX
    # OpenMP::OpenMP_CXX is a modern imported target. It automatically adds:
    #   - -fopenmp to compile flags
    #   - -fopenmp to link flags (needed to link the OpenMP runtime libgomp)
    # Using this target is safer than manually appending to CMAKE_CXX_FLAGS
    # because it handles both compile and link in one place.

    Eigen3
    # Our INTERFACE imported target — propagates the include path only,
    # no actual library to link since Eigen is header-only.

    # $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    # Linker also needs the sanitizer flag to link the ASan/UBSan runtime.
)

# ── Status messages ───────────────────────────────────────────────────────────
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "TARGET_ARCH:       ${TARGET_ARCH}")
message(STATUS "AVX-512:           ${HAS_AVX512}")
message(STATUS "OpenMP found:      ${OpenMP_CXX_FOUND} ${OpenMP_CXX_VERSION}")




#Previous CmakeLists.txt content (before edits):
# cmake_minimum_required(VERSION 3.10)
# project(MyProject LANGUAGES CXX)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# # Find OpenMP — sets OpenMP_CXX_FOUND, provides imported target
# find_package(OpenMP REQUIRED)

# # Include headers
# include_directories("include")
# include_directories("/home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0") #include Eigen headers

# # Gather all source files
# file(GLOB SOURCES "src/*.cpp")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -funroll-loops") #add OpenMP flag for parallelization

# # Add executable with main.cpp + all src files
# add_executable(main main.cpp ${SOURCES})

# # Link curl
# target_link_libraries(main PRIVATE curl OpenMP::OpenMP_CXX)