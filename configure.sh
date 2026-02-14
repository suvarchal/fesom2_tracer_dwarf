#!/bin/bash
#
# Configure and build FESOM2 Tracer Dwarf
#
# Usage:
#   ./configure.sh [OPTIONS]
#
# Options:
#   --compiler COMPILER   gnu (default), intel, nvidia
#   --precision PREC      dp (default), sp, hp
#   --openmp              Enable OpenMP (default: off)
#   --openacc             Enable OpenACC (default: off)
#   --build-type TYPE     Release (default) or Debug
#   --clean               Remove build directory before configuring
#   --build               Also run make after configuring
#   --help                Show this help
#
# Build directories are named: build_<compiler>_<precision>
# e.g. build_gnu_dp, build_intel_sp, build_nvidia_dp
#
# Examples:
#   ./configure.sh                                  # GNU double precision
#   ./configure.sh --compiler intel --precision sp   # Intel single precision
#   ./configure.sh --compiler nvidia --build         # NVIDIA DP + build
#   ./configure.sh --compiler gnu --precision sp --build
#
# For convenience, a symlink 'build' -> latest build directory is created.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# ========================================
# Defaults
# ========================================
COMPILER="gnu"
PRECISION="dp"
BUILD_TYPE="Release"
ENABLE_OPENMP="OFF"
ENABLE_OPENACC="OFF"
DO_CLEAN=false
DO_BUILD=false

# ========================================
# Parse arguments
# ========================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --compiler)
            COMPILER="$2"
            shift 2
            ;;
        --precision)
            PRECISION="$2"
            shift 2
            ;;
        --build-type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        --openmp)
            ENABLE_OPENMP="ON"
            shift
            ;;
        --openacc)
            ENABLE_OPENACC="ON"
            shift
            ;;
        --clean)
            DO_CLEAN=true
            shift
            ;;
        --build)
            DO_BUILD=true
            shift
            ;;
        --help|-h)
            head -28 "$0" | tail -26
            exit 0
            ;;
        *)
            echo "Error: unknown option '$1'"
            echo "Run with --help for usage"
            exit 1
            ;;
    esac
done

# ========================================
# Validate inputs
# ========================================
case "$COMPILER" in
    gnu|intel|nvidia) ;;
    *)
        echo "Error: unknown compiler '$COMPILER' (use: gnu, intel, nvidia)"
        exit 1
        ;;
esac

case "$PRECISION" in
    dp|sp|hp) ;;
    *)
        echo "Error: unknown precision '$PRECISION' (use: dp, sp, hp)"
        exit 1
        ;;
esac

# ========================================
# Set compiler paths and CMake variables
# ========================================
CMAKE_EXTRA_ARGS=""
MPIRUN="mpirun"

case "$COMPILER" in
    gnu)
        # Use system default gfortran/gcc
        FC=$(which gfortran 2>/dev/null || true)
        CC=$(which gcc 2>/dev/null || true)
        if [ -z "$FC" ]; then
            echo "Error: gfortran not found in PATH"
            exit 1
        fi
        ;;
    intel)
        INTEL_ROOT="/opt/intel/oneapi"
        if [ -f "$INTEL_ROOT/setvars.sh" ]; then
            echo "Sourcing Intel oneAPI environment..."
            # Must not pipe source output (pipe runs in subshell, loses env)
            source "$INTEL_ROOT/setvars.sh" --force > /dev/null 2>&1
        fi
        FC=$(which ifx 2>/dev/null || echo "$INTEL_ROOT/compiler/latest/bin/ifx")
        CC=$(which icx 2>/dev/null || echo "$INTEL_ROOT/compiler/latest/bin/icx")
        if [ ! -x "$FC" ]; then
            echo "Error: ifx not found at $FC"
            exit 1
        fi
        # Intel MPI: use wrappers directly as compilers so CMake calls
        # 'mpiifx -show' at runtime (which resolves I_MPI_ROOT correctly)
        # instead of parsing the wrapper script file.
        MPI_FC=$(which mpiifx 2>/dev/null || echo "$INTEL_ROOT/mpi/latest/bin/mpiifx")
        MPI_CC=$(which mpiicx 2>/dev/null || echo "$INTEL_ROOT/mpi/latest/bin/mpiicx")
        if [ -x "$MPI_FC" ]; then
            # Use MPI wrappers AS the compilers — they embed ifx/icx + MPI flags
            FC="$MPI_FC"
            CC="$MPI_CC"
        fi
        ;;
    nvidia)
        NVHPC_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/2025"
        FC="$NVHPC_ROOT/compilers/bin/nvfortran"
        CC="$NVHPC_ROOT/compilers/bin/nvc"
        if [ ! -x "$FC" ]; then
            echo "Error: nvfortran not found at $FC"
            exit 1
        fi
        # System OpenMPI's Fortran modules are compiled with gfortran and
        # can't be used by nvfortran. Use NVIDIA's bundled MPI instead.
        NV_MPI_FC="$NVHPC_ROOT/comm_libs/mpi/bin/mpifort"
        NV_MPI_CC="$NVHPC_ROOT/comm_libs/mpi/bin/mpicc"
        NV_MPIRUN="$NVHPC_ROOT/comm_libs/mpi/bin/mpirun"
        if [ -x "$NV_MPI_FC" ]; then
            FC="$NV_MPI_FC"
            CC="$NV_MPI_CC"
            MPIRUN="$NV_MPIRUN"
        fi
        ;;
esac

# Precision flags
USE_SINGLE="OFF"
USE_HALF="OFF"
if [ "$PRECISION" = "sp" ]; then
    USE_SINGLE="ON"
elif [ "$PRECISION" = "hp" ]; then
    USE_HALF="ON"
fi

# Build directory
BUILD_DIR="build_${COMPILER}_${PRECISION}"

# ========================================
# Print configuration
# ========================================
echo ""
echo "========================================="
echo "FESOM2 Tracer Dwarf Configuration"
echo "========================================="
echo "  Compiler:    $COMPILER ($FC)"
echo "  Precision:   $PRECISION"
echo "  Build type:  $BUILD_TYPE"
echo "  OpenMP:      $ENABLE_OPENMP"
echo "  OpenACC:     $ENABLE_OPENACC"
echo "  Build dir:   $BUILD_DIR"
echo "========================================="
echo ""

# ========================================
# Clean if requested
# ========================================
if $DO_CLEAN && [ -d "$BUILD_DIR" ]; then
    echo "Removing existing $BUILD_DIR/"
    rm -rf "$BUILD_DIR"
fi

# ========================================
# Configure
# ========================================
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake .. \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_Fortran_COMPILER="$FC" \
    -DCMAKE_C_COMPILER="$CC" \
    -DUSE_SINGLE_PRECISION="$USE_SINGLE" \
    -DUSE_HALF_PRECISION="$USE_HALF" \
    -DENABLE_OPENMP="$ENABLE_OPENMP" \
    -DENABLE_OPENACC="$ENABLE_OPENACC" \
    $CMAKE_EXTRA_ARGS

echo ""
echo "========================================="
echo "Configuration complete!"
echo "========================================="

# ========================================
# Build if requested
# ========================================
if $DO_BUILD; then
    echo ""
    echo "Building..."
    make -j$(nproc)
    echo ""
    echo "========================================="
    echo "Build complete!"
    echo "  Executable: $BUILD_DIR/bin/fesom_tracer_analytic"
    echo "  Library:    $BUILD_DIR/lib/libfesom_tracer_Fortran.so"
    echo "========================================="
else
    echo ""
    echo "To build, run:"
    echo "  cd $BUILD_DIR && make -j\$(nproc)"
fi

# Create/update 'build' symlink to latest build
cd "$SCRIPT_DIR"
if [ -L "build" ] || [ ! -e "build" ]; then
    ln -sfn "$BUILD_DIR" build
    echo ""
    echo "Symlink: build -> $BUILD_DIR"
fi

# ========================================
# Generate run.sh wrapper with correct mpirun
# ========================================
cd "$SCRIPT_DIR/$BUILD_DIR"
cat > run.sh <<RUNEOF
#!/bin/bash
# Auto-generated run wrapper for $BUILD_DIR
# Uses the correct mpirun for the $COMPILER compiler
#
# Usage: ./run.sh NP [program-args...]
#   NP = number of MPI processes
#
# Environment:
#   OMP_NUM_THREADS  Number of OpenMP threads per MPI rank (default: 1)
#
# Examples:
#   ./run.sh 1 20 20 10 --periodic
#   OMP_NUM_THREADS=4 ./run.sh 1 20 20 10 --periodic
#   ./run.sh 1 50 50 10 --save-mesh --save-scalars --periodic

set -e
if [ \$# -lt 1 ]; then
    echo "Usage: ./run.sh NP [program-args...]"
    echo "  NP = number of MPI processes"
    exit 1
fi

NP="\$1"
shift

export OMP_NUM_THREADS="\${OMP_NUM_THREADS:-1}"

SCRIPT_DIR="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )"
cd "\$SCRIPT_DIR"
exec $MPIRUN -np "\$NP" ./bin/fesom_tracer_analytic "\$@"
RUNEOF
chmod +x run.sh

echo ""
echo "Run the analytic dwarf:"
echo "  cd $BUILD_DIR && ./run.sh 1 20 20 10 --periodic"
echo ""
echo "Or directly:"
echo "  cd $BUILD_DIR && $MPIRUN -np 1 ./bin/fesom_tracer_analytic 20 20 10 --periodic"
echo ""
