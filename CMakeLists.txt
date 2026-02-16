cmake_minimum_required(VERSION 3.16)
project(fesom_tracer C Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
include(fesom_setup)

# Option to enable C interface for Python binding
option(ENABLE_TRACER_C_INTERFACE "Build C interface for tracer library" OFF)

# Precision options (mutually exclusive)
option(USE_HALF_PRECISION "Build with half precision (WP=2, FP16)" OFF)
option(USE_SINGLE_PRECISION "Build with single precision (WP=4)" OFF)
# Double precision (WP=8) is the default when both options are OFF

# Option to enable OpenMP support
option(ENABLE_OPENMP "compile with OpenMP support" OFF)

# Option to enable OpenACC support
option(ENABLE_OPENACC "compile with OpenACC support" ON)
option(DISABLE_OPENACC_ATOMICS "disable kernels using atomic statement for reproducible results" OFF)
option(ENABLE_NVTX "enable use of nvtx tools" OFF)

# Find required packages
find_package(MPI REQUIRED)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP enabled: ${OpenMP_Fortran_VERSION}")
endif()

# Determine library directory - works for both standalone and integrated builds
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # This is a standalone build
    set(lib_home ${CMAKE_CURRENT_SOURCE_DIR}/lib)
    message(STATUS "Building as standalone tracer dwarf")

    # Set default install prefix if not specified
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../.. CACHE PATH "Install path prefix" FORCE)
    endif()

    # Add uninstall target
    if(NOT TARGET uninstall)
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
            IMMEDIATE @ONLY)
        add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    endif()
else()
    # This is part of the main FESOM build
    set(lib_home ${CMAKE_CURRENT_SOURCE_DIR}/lib)
    message(STATUS "Building as tracer build")
endif()

message(STATUS "Library directory: ${lib_home}")

# ========================================
# Precision validation
# ========================================
if(USE_HALF_PRECISION AND USE_SINGLE_PRECISION)
    message(FATAL_ERROR "Cannot enable both USE_HALF_PRECISION and USE_SINGLE_PRECISION")
endif()

if(USE_HALF_PRECISION)
    message(STATUS "Building with HALF precision (WP=2, FP16)")
    message(WARNING "Half precision support is experimental and requires compiler support for REAL(2)")
elseif(USE_SINGLE_PRECISION)
    message(STATUS "Building with SINGLE precision (WP=4)")
else()
    message(STATUS "Building with DOUBLE precision (WP=8, default)")
endif()

# ========================================
# Unified compiler flag function
# ========================================
# Applies correct compiler flags and precision definitions to a target.
# This eliminates the per-target flag duplication that previously caused
# precision mismatches between library and executables.
function(apply_fesom_compile_flags target_name)
    # Precision preprocessor definitions
    if(USE_HALF_PRECISION)
        target_compile_definitions(${target_name} PRIVATE USE_HALF_PRECISION)
    elseif(USE_SINGLE_PRECISION)
        target_compile_definitions(${target_name} PRIVATE USE_SINGLE_PRECISION)
    endif()

    # Determine precision-specific real-size flag
    if(USE_HALF_PRECISION)
        set(_intel_real_flag "-r2")
        set(_nvhpc_real_flag "-Mr2")
        set(_cray_real_flag "-s real16")
    elseif(USE_SINGLE_PRECISION)
        set(_intel_real_flag "-r4")
        set(_nvhpc_real_flag "-Mr4")
        set(_cray_real_flag "-s real32")
    else()
        set(_intel_real_flag "-r8")
        set(_nvhpc_real_flag "-Mr8")
        set(_cray_real_flag "-s real64")
    endif()

    # --- Intel Classic (ifort) or Intel LLVM (ifx) ---
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        if(USE_HALF_PRECISION)
            message(FATAL_ERROR "Intel Fortran does not support REAL(kind=2). Use NVIDIA (nvfortran) for half precision.")
        endif()
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
            # ifx flags (LLVM-based): no -ip, -no-prec-sqrt, -fimf-use-svml
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:DEBUG>:${_intel_real_flag} -i4 -fp-model precise -no-prec-div -xHost -g -traceback -check all,noarg_temp_created,bounds,uninit>
                $<$<CONFIG:RELEASE>:${_intel_real_flag} -i4 -fp-model precise -no-prec-div -g -traceback>
            )
        else()
            # ifort flags (classic)
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:DEBUG>:${_intel_real_flag} -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit>
                $<$<CONFIG:RELEASE>:${_intel_real_flag} -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -g -traceback>
            )
        endif()
        if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
            if($ENV{FESOM_PLATFORM_STRATEGY} STREQUAL "levante.dkrz.de")
                target_compile_options(${target_name} PRIVATE -march=core-avx2 -mtune=core-avx2)
            else()
                target_compile_options(${target_name} PRIVATE -xHost)
            endif()
        else()
            target_compile_options(${target_name} PRIVATE -xHost)
        endif()

    # --- GNU (gfortran) ---
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        if(USE_HALF_PRECISION)
            message(FATAL_ERROR "GNU Fortran does not support native REAL(2)/FP16.")
        endif()
        if(USE_SINGLE_PRECISION)
            # Single precision: no -fdefault-real flag, WP=4 handles it via preprocessor
            target_compile_options(${target_name} PRIVATE
                -O3 -finit-local-zero -finline-functions -march=native -fimplicit-none -ffree-line-length-none
            )
        else()
            target_compile_options(${target_name} PRIVATE
                -O3 -finit-local-zero -finline-functions -march=native -fimplicit-none -fdefault-real-8 -ffree-line-length-none
            )
        endif()
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
            target_compile_options(${target_name} PRIVATE -fallow-argument-mismatch)
        endif()

    # --- Cray ---
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        if(USE_HALF_PRECISION)
            message(FATAL_ERROR "Cray Fortran does not support native REAL(2)/FP16.")
        endif()
        target_compile_options(${target_name} PRIVATE
            -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am ${_cray_real_flag} -hnoomp
        )

    # --- NVIDIA HPC SDK (nvfortran) or legacy PGI ---
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        # nvfortran uses -r4/-r8, legacy PGI uses -Mr4/-Mr8
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
            if(USE_HALF_PRECISION)
                # No -r flag for HP; WP=2 is set via USE_HALF_PRECISION preprocessor define
                set(_nv_real_flag "")
            elseif(USE_SINGLE_PRECISION)
                set(_nv_real_flag "-r4")
            else()
                set(_nv_real_flag "-r8")
            endif()
        else()
            set(_nv_real_flag "${_nvhpc_real_flag}")
        endif()
        if(ENABLE_OPENACC)
            if(USE_HALF_PRECISION)
                message(WARNING "Half precision on NVHPC requires GPU with FP16 support (Volta/Turing/Ampere+)")
            endif()
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:DEBUG>:-Mallocatable=95 ${_nv_real_flag} -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=math_uniform,cuda11.7,cc80>
                $<$<CONFIG:RELEASE>:-Mallocatable=95 ${_nv_real_flag} -pgf90libs -Minfo=all -acc=verystrict -gpu=cc80>
            )
        else()
            if(USE_HALF_PRECISION)
                message(WARNING "Half precision on NVHPC without OpenACC uses software FP16 emulation.")
            endif()
            target_compile_options(${target_name} PRIVATE
                $<$<CONFIG:DEBUG>:-Mallocatable=95 ${_nv_real_flag} -pgf90libs -Mnofma>
                $<$<CONFIG:RELEASE>:-Mallocatable=95 ${_nv_real_flag} -pgf90libs>
            )
        endif()
        if(ENABLE_NVTX)
            target_link_libraries(${target_name} nvhpcwrapnvtx)
        endif()
    endif()
endfunction()

# ========================================
# Library sources
# ========================================
set(library_sources_Fortran
    ${lib_home}/hp_math_intrinsics.F90
    ${lib_home}/oce_modules.F90
    ${lib_home}/gen_modules_config.F90
    ${lib_home}/gen_modules_rotate_grid.F90
    ${lib_home}/gen_modules_partitioning.F90
    ${lib_home}/fortran_utils.F90
    ${lib_home}/MOD_MESH.F90
    ${lib_home}/MOD_PARTIT.F90
    ${lib_home}/MOD_TRACER.F90
    ${lib_home}/MOD_DYN.F90
    ${lib_home}/MOD_READ_BINARY_ARRAYS.F90
    ${lib_home}/MOD_WRITE_BINARY_ARRAYS.F90
    ${lib_home}/gen_halo_exchange.F90
    ${lib_home}/io_restart_derivedtype.F90
    ${lib_home}/oce_mesh.F90
    ${lib_home}/oce_edge_coloring.F90
    ${lib_home}/oce_muscl_adv.F90
    ${lib_home}/oce_adv_tra_fct.F90
    ${lib_home}/oce_adv_tra_hor.F90
    ${lib_home}/oce_adv_tra_ver.F90
    ${lib_home}/oce_adv_tra_driver.F90
    ${lib_home}/analytic_mesh.F90
    ${lib_home}/mesh_output.F90
    ${lib_home}/tracer_init_from_mesh.F90
    ${lib_home}/fesom_profiler.F90
)

# Add C interface source if enabled
if(ENABLE_TRACER_C_INTERFACE)
    list(APPEND library_sources_Fortran ${lib_home}/tracer_c_interface.F90)
    message(STATUS "C interface for tracer library enabled")
endif()

# ========================================
# Shared library
# ========================================
add_library(${PROJECT_NAME}_Fortran SHARED ${library_sources_Fortran})
apply_fesom_compile_flags(${PROJECT_NAME}_Fortran)
target_link_libraries(${PROJECT_NAME}_Fortran PUBLIC MPI::MPI_Fortran)
if(ENABLE_OPENMP)
    target_link_libraries(${PROJECT_NAME}_Fortran PUBLIC OpenMP::OpenMP_Fortran)
endif()
set_target_properties(${PROJECT_NAME}_Fortran PROPERTIES LINKER_LANGUAGE Fortran)

# OpenACC linker flags for NVHPC
if(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    if(ENABLE_OPENACC)
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "-acc=verystrict -Mnofma -gpu=math_uniform,cuda11.7,cc80")
        set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "-acc=verystrict -gpu=cc80")
    endif()
    if(DISABLE_OPENACC_ATOMICS)
        target_compile_definitions(${PROJECT_NAME}_Fortran PRIVATE DISABLE_OPENACC_ATOMICS)
    endif()
endif()

# ========================================
# Main executable (file-based mesh)
# ========================================
add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/fesom.F90)
apply_fesom_compile_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_Fortran)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# ========================================
# Mesh-init executable
# ========================================
set(mesh_init_executable ${PROJECT_NAME}_mesh_init)
add_executable(${mesh_init_executable} ${CMAKE_CURRENT_SOURCE_DIR}/src/fesom_mesh_init.F90)
apply_fesom_compile_flags(${mesh_init_executable})
target_link_libraries(${mesh_init_executable} PUBLIC ${PROJECT_NAME}_Fortran)
set_target_properties(${mesh_init_executable} PROPERTIES LINKER_LANGUAGE Fortran)

# ========================================
# Analytic mesh executable
# ========================================
set(analytic_executable ${PROJECT_NAME}_analytic)
add_executable(${analytic_executable} ${CMAKE_CURRENT_SOURCE_DIR}/src/fesom_analytic.F90)
apply_fesom_compile_flags(${analytic_executable})
target_link_libraries(${analytic_executable} PUBLIC ${PROJECT_NAME}_Fortran)
set_target_properties(${analytic_executable} PROPERTIES LINKER_LANGUAGE Fortran)

# OpenACC linker flags for executables (NVHPC)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    if(ENABLE_OPENACC)
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc=verystrict -Mnofma -gpu=math_uniform,cuda11.7,cc80")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc=verystrict -gpu=cc80")
    endif()
endif()

# ========================================
# Install targets
# ========================================
install(TARGETS ${PROJECT_NAME} ${mesh_init_executable} ${analytic_executable}
    RUNTIME DESTINATION bin
    COMPONENT tools
)

install(TARGETS ${PROJECT_NAME}_Fortran
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    COMPONENT tools
)
