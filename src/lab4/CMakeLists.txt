cmake_minimum_required(VERSION 3.20)
project(amss_ncku_lab LANGUAGES CXX Fortran)

# --- Options ----------------------------------------------------------------
option(AMSS_ENABLE_GPU "Enable CUDA GPU build (ABEGPU target)" OFF)
set(AMSS_OPT "-O3" CACHE STRING "C/C++/Fortran optimization level")
set(AMSS_ARCH_FLAGS "" CACHE STRING "Architecture-specific compiler flags")
option(AMSS_ENABLE_OPENMP "Link OpenMP to all executables" ON)
option(AMSS_TWOP_OPENMP "Link OpenMP to the standalone TwoPuncture solver" ON)
set(AMSS_TWOP_DEFAULT_NRELAX "10" CACHE STRING "Default TwoPuncture line-preconditioner sweeps")
option(AMSS_FAST_MATH "Enable -ffast-math for the CPU ABE target only" ON)
option(AMSS_POLINT_UNIFORM6 "Use a fixed six-point uniform interpolation kernel in global_interp" ON)
option(AMSS_SOFTWARE_PREFETCH "Enable compiler-generated software prefetch loops" OFF)
option(AMSS_RHS_NAN_CHECK "Scan RHS input arrays for NaN on every call" OFF)
option(AMSS_FINE_PREDICTOR_EVOLVE "Use the evolution-only RHS mode for fine-level predictors" ON)
option(AMSS_RHS_WORKSPACE_POOL "Pool legacy fused-RHS stencil arrays outside the POINTWISE hot path" ON)
option(AMSS_VACUUM_BSSN "Compile the CPU BSSN RHS without stress-energy loads" ON)
option(AMSS_FUSED_RHS_TAIL "Fuse lopsided and KO scans in the B=2 batch path" ON)
option(AMSS_BATCH_STENCIL "Enable batched lopsided/KO stencils in the CPU BSSN path" ON)
option(AMSS_BATCH_STENCIL_SINGLE_REGION "Use one OpenMP worksharing region per batched stencil" OFF)
option(AMSS_RHS_BULK_SIMD "Request SIMD on the fused RHS i-direction loop" OFF)
option(AMSS_RHS_SKIP_CONSTRAINT_STORES "Skip fused RHS constraint-only array stores when co=1" OFF)
option(AMSS_RHS_SKIP_CHIN1_SCAN "Skip the fused RHS chin1 auxiliary scan when co=1" ON)
option(AMSS_RHS_RAW_DIAG_LOPSIDED "Use dxx/dyy/dzz directly for batch lopsided derivatives" ON)
option(AMSS_FGFS_COLORED_SLAB "Use cache-colored independent CPU fgfs allocations" OFF)
option(AMSS_FGFS_HUGEPAGE_SLAB "Pool CPU fgfs into one 2MiB-aligned MADV_HUGEPAGE slab (THP for dTLB)" OFF)
set(AMSS_RHS_TILE_J "0" CACHE STRING "POINTWISE RHS j tile height (0 disables tiling)")
set(AMSS_RHS_TILE_K "0" CACHE STRING "POINTWISE RHS k tile depth (0 disables tiling)")
set(AMSS_RHS_LOCALITY "POINTWISE" CACHE STRING "RHS locality path: LEGACY, HALO, or POINTWISE")
set_property(CACHE AMSS_RHS_LOCALITY PROPERTY STRINGS LEGACY HALO POINTWISE)
option(AMSS_BATCH_OUTBD "Batch AMR coarse-to-fine OutBd transfers in the CPU ABE path" OFF)
option(AMSS_LEVEL_DIAG "Enable RecursiveStep/Constraint_Out timing diagnostics" OFF)
option(AMSS_CACHE_COMM_PLANS "Cache CPU AMR communication segment plans" OFF)
option(AMSS_ASYNC_LEVEL_SYNC "Defer final level Sync while child AMR levels compute" OFF)
option(AMSS_ENFORCE_GA_POINTWISE "Use a scalar pointwise enforce_ga kernel" OFF)
option(AMSS_RK4_FUSE_LOWERBOUND "Fuse the phi lower bound into its RK4 update" OFF)
option(AMSS_RK4_BATCHED "Use a tiled multi-field CPU RK4 update" OFF)
option(AMSS_SOMMERFELD_BATCHED "Batch CPU Sommerfeld boundary updates" OFF)
set(AMSS_MPI_CUDA_AWARE 0 CACHE STRING "Set MPI_CUDA_AWARE macro for ABEGPU (0/1)")
set_property(CACHE AMSS_MPI_CUDA_AWARE PROPERTY STRINGS 0 1)

# --- Compile flags ----------------------------------------------------------
# Legacy AMSS-NCKU has strict-aliasing UB; -fno-strict-aliasing is kept
# regardless of AMSS_OPT so changing -O does not lose the safety flag.
# `-Dfortran3 -Dnewc` select the modern C/Fortran source branches (without
# them the C++ side falls back to pre-ANSI <iostream.h> includes and the
# Fortran side skips the f90 modules). Apply to every supported compiler.
separate_arguments(AMSS_OPT_LIST UNIX_COMMAND "${AMSS_OPT}")
separate_arguments(AMSS_ARCH_FLAGS_LIST UNIX_COMMAND "${AMSS_ARCH_FLAGS}")

# Clang/IntelLLVM/ARMClang all accept the GCC-style flags used here.
# COMPILE_LANG_AND_ID takes comma-separated compiler IDs (not CMake list).
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,IntelLLVM,Intel,ARMClang>:${AMSS_OPT_LIST};${AMSS_ARCH_FLAGS_LIST};-fno-strict-aliasing;-Wno-deprecated;-Dfortran3;-Dnewc>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU,Flang,IntelLLVM,Intel,ARMClang>:${AMSS_OPT_LIST};${AMSS_ARCH_FLAGS_LIST};-fno-strict-aliasing;-cpp>"
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- MPI --------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS CXX)

# --- Shared Fortran sources -------------------------------------------------
set(AMSS_F90SOURCES
  src/bssn_rhs.f90
  src/batch_stencils.f90
  src/diff_new.f90 src/enforce_algebra.f90
  src/fadmquantites_bssn.f90 src/fmisc.f90
  src/getnp4.f90
  src/initial_puncture.f90 src/kodiss.f90
  src/lopsidediff.f90
  src/prolongrestrict_cell.f90 src/rungekutta4_rout.f90
  src/sommerfeld_rout.f90
)

# --- ABE: CPU main evolution executable -------------------------------------
set(ABE_SOURCES
  src/ABE.C src/Ansorg.C src/Block.C src/MPatch.C src/Parallel.C
  src/bssn_class.C src/cgh.C src/checkpoint.C
  src/misc.C src/monitor.C src/perf.C
  src/surface_integral.C src/transpbh.C src/var.C
  ${AMSS_F90SOURCES}
)
add_executable(ABE ${ABE_SOURCES})
target_link_libraries(ABE PRIVATE MPI::MPI_CXX)
set(AMSS_FORTRAN_MODULE_ROOT "${CMAKE_CURRENT_BINARY_DIR}/fortran-mod")
set_target_properties(ABE PROPERTIES
  Fortran_MODULE_DIRECTORY "${AMSS_FORTRAN_MODULE_ROOT}/ABE")
if(AMSS_RHS_LOCALITY STREQUAL "HALO")
  target_compile_definitions(ABE PRIVATE AMSS_RHS_HALO)
elseif(AMSS_RHS_LOCALITY STREQUAL "POINTWISE")
  target_compile_definitions(ABE PRIVATE AMSS_RHS_HALO AMSS_RHS_POINTWISE)
  target_compile_options(ABE PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none;-finline-functions>")
elseif(NOT AMSS_RHS_LOCALITY STREQUAL "LEGACY")
  message(FATAL_ERROR "AMSS_RHS_LOCALITY must be LEGACY, HALO, or POINTWISE")
endif()
if(AMSS_FAST_MATH)
  # Keep the experiment scoped to the evolution executable.  In particular,
  # do not change TwoPunctureABE's initial-data solver semantics or timing.
  target_compile_options(ABE PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:-ffast-math>"
    "$<$<COMPILE_LANGUAGE:Fortran>:-ffast-math>")
endif()
if(AMSS_POLINT_UNIFORM6)
  # The specialized kernel is deliberately limited to the CPU ABE target.
  # It is selected only by global_interp when ordn==6; AMR and Sommerfeld
  # interpolation retain the legacy polint/polin3 path.
  target_compile_definitions(ABE PRIVATE AMSS_POLINT_UNIFORM6)
endif()
if(AMSS_SOFTWARE_PREFETCH)
  target_compile_options(ABE PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-fprefetch-loop-arrays>")
endif()
if(AMSS_BATCH_OUTBD)
  target_compile_definitions(ABE PRIVATE AMSS_BATCH_OUTBD)
endif()
if(AMSS_LEVEL_DIAG)
  target_compile_definitions(ABE PRIVATE AMSS_LEVEL_DIAG)
endif()
if(AMSS_CACHE_COMM_PLANS)
  target_compile_definitions(ABE PRIVATE AMSS_CACHE_COMM_PLANS)
endif()
if(AMSS_ASYNC_LEVEL_SYNC)
  target_compile_definitions(ABE PRIVATE AMSS_ASYNC_LEVEL_SYNC)
endif()
if(AMSS_ENFORCE_GA_POINTWISE)
  target_compile_definitions(ABE PRIVATE AMSS_ENFORCE_GA_POINTWISE)
endif()
if(AMSS_RK4_FUSE_LOWERBOUND)
  target_compile_definitions(ABE PRIVATE AMSS_RK4_FUSE_LOWERBOUND)
endif()
if(AMSS_RK4_BATCHED)
  target_compile_definitions(ABE PRIVATE AMSS_RK4_BATCHED)
endif()
if(AMSS_SOMMERFELD_BATCHED)
  target_compile_definitions(ABE PRIVATE AMSS_SOMMERFELD_BATCHED)
endif()
if(AMSS_BATCH_STENCIL)
  target_compile_definitions(ABE PRIVATE AMSS_BATCH_STENCIL)
endif()
if(AMSS_BATCH_STENCIL_SINGLE_REGION)
  target_compile_definitions(ABE PRIVATE AMSS_BATCH_STENCIL_SINGLE_REGION)
endif()
if(AMSS_RHS_REORDER)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_REORDER)
endif()
if(AMSS_RHS_NAN_CHECK)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_NAN_CHECK)
endif()
if(AMSS_FINE_PREDICTOR_EVOLVE)
  target_compile_definitions(ABE PRIVATE AMSS_FINE_PREDICTOR_EVOLVE)
endif()
if(AMSS_RHS_WORKSPACE_POOL)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_WORKSPACE_POOL)
endif()
if(AMSS_VACUUM_BSSN)
  target_compile_definitions(ABE PRIVATE AMSS_VACUUM_BSSN)
endif()
if(AMSS_FUSED_RHS_TAIL)
  target_compile_definitions(ABE PRIVATE AMSS_FUSED_RHS_TAIL)
endif()
if(AMSS_RHS_BULK_SIMD)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_BULK_SIMD)
endif()
if(AMSS_RHS_OMP_ASSEMBLY)
  # Requires AMSS_ENABLE_OPENMP (ON by default → OpenMP::OpenMP_Fortran linked to ABE).
  # The directive's if() clause keeps co==0/octant and OMP=1 paths fully serial, so
  # the #37 baseline is byte-identical when this stays OFF (CMake default).
  target_compile_definitions(ABE PRIVATE AMSS_RHS_OMP_ASSEMBLY)
endif()
if(AMSS_RHS_SKIP_CONSTRAINT_STORES)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_SKIP_CONSTRAINT_STORES)
endif()
if(AMSS_RHS_SKIP_CHIN1_SCAN)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_SKIP_CHIN1_SCAN)
endif()
if(AMSS_RHS_RAW_DIAG_LOPSIDED)
  target_compile_definitions(ABE PRIVATE AMSS_RHS_RAW_DIAG_LOPSIDED)
endif()
if(AMSS_FGFS_COLORED_SLAB)
  target_compile_definitions(ABE PRIVATE AMSS_FGFS_COLORED_SLAB)
endif()
if(AMSS_FGFS_HUGEPAGE_SLAB)
  target_compile_definitions(ABE PRIVATE AMSS_FGFS_HUGEPAGE_SLAB)
endif()
if(AMSS_RHS_TILE_J GREATER 0 AND AMSS_RHS_TILE_K GREATER 0)
  target_compile_definitions(ABE PRIVATE
    AMSS_RHS_TILED
    AMSS_RHS_TILE_J_VALUE=${AMSS_RHS_TILE_J}
    AMSS_RHS_TILE_K_VALUE=${AMSS_RHS_TILE_K})
endif()
if(AMSS_GREEDY_BALANCE)
  target_compile_definitions(ABE PRIVATE AMSS_GREEDY_BALANCE)
endif()

# --- TwoPunctureABE: initial-data solver -----------------------------------
add_executable(TwoPunctureABE src/TwoPunctureABE.C src/TwoPunctures.C)
target_link_libraries(TwoPunctureABE PRIVATE MPI::MPI_CXX)
target_compile_definitions(TwoPunctureABE PRIVATE AMSS_TWOP_DEFAULT_NRELAX=${AMSS_TWOP_DEFAULT_NRELAX})

# --- Optional OpenMP --------------------------------------------------------
if(AMSS_ENABLE_OPENMP OR AMSS_TWOP_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX Fortran)
endif()
if(AMSS_ENABLE_OPENMP)
  target_link_libraries(ABE PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran)
endif()
if(AMSS_ENABLE_OPENMP OR AMSS_TWOP_OPENMP)
  target_link_libraries(TwoPunctureABE PRIVATE OpenMP::OpenMP_CXX)
endif()

# --- ABEGPU: GPU evolution executable (device-resident baseline) -----------
# Same C++/Fortran sources as ABE except bssn_class.C → bssn_gpu_class.C,
# plus GPU kernel/agent files.
set(ABEGPU_CXX_SOURCES
  src/ABE.C src/Ansorg.C src/Block.C src/MPatch.C src/Parallel.C
  src/bssn_gpu_class.C src/bssn_step_gpu.C
  src/cgh.C src/checkpoint.C src/misc.C
  src/monitor.C src/perf.C src/surface_integral.C
  src/transpbh.C src/var.C
  src/helper.cpp src/Parallel_GPU.cpp
)
set(ABEGPU_CUDA_SOURCES
  src/fmisc_gpu.cu src/lopsidediff_gpu.cu src/kodiss_gpu.cu src/diff_new_gpu.cu
  src/bssn_rhs_gpu.cu src/prolongrestrict_cell_gpu.cu src/gpu_manager.cu
  src/enforce_algebra_gpu.cu src/sommerfeld_rout_gpu.cu src/rungekutta4_rout_gpu.cu
  src/fadmquantites_bssn_gpu.cu src/misc_gpu.cu src/MPatch_gpu.cu
  src/surface_integral_gpu.cu src/getnp4_gpu.cu
)

if(AMSS_ENABLE_GPU)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 80 CACHE STRING "CUDA architectures" FORCE)
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)

  add_executable(ABEGPU ${ABEGPU_CXX_SOURCES} ${ABEGPU_CUDA_SOURCES} ${AMSS_F90SOURCES})
  set_target_properties(ABEGPU PROPERTIES
    Fortran_MODULE_DIRECTORY "${AMSS_FORTRAN_MODULE_ROOT}/ABEGPU")
  set_source_files_properties(${ABEGPU_CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_target_properties(ABEGPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_compile_options(ABEGPU PRIVATE
    "$<$<COMPILE_LANGUAGE:CUDA>:-rdc=true;-lineinfo>")
  target_compile_definitions(ABEGPU PRIVATE USE_GPU fortran3 newc MPI_CUDA_AWARE=${AMSS_MPI_CUDA_AWARE})
  target_include_directories(ABEGPU PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(ABEGPU PRIVATE CUDA::cudart MPI::MPI_CXX)
  if(AMSS_ENABLE_OPENMP)
    target_link_libraries(ABEGPU PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran)
  endif()
endif()

# --- Configure summary ------------------------------------------------------
message(STATUS "==================== AMSS-NCKU Lab 4 ====================")
message(STATUS "CXX      : ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran  : ${CMAKE_Fortran_COMPILER}")
if(AMSS_ENABLE_GPU)
  message(STATUS "GPU      : ON, CUDA=${CMAKE_CUDA_COMPILER}, arch=${CMAKE_CUDA_ARCHITECTURES}")
else()
  message(STATUS "GPU      : OFF")
endif()
message(STATUS "AMSS_OPT=${AMSS_OPT} ARCH_FLAGS=${AMSS_ARCH_FLAGS} OPENMP=${AMSS_ENABLE_OPENMP} TWOP_OPENMP=${AMSS_TWOP_OPENMP} TWOP_DEFAULT_NRELAX=${AMSS_TWOP_DEFAULT_NRELAX} FAST_MATH=${AMSS_FAST_MATH} POLINT_UNIFORM6=${AMSS_POLINT_UNIFORM6} RHS_LOCALITY=${AMSS_RHS_LOCALITY} NAN_CHECK=${AMSS_RHS_NAN_CHECK} FINE_PREDICTOR_EVOLVE=${AMSS_FINE_PREDICTOR_EVOLVE} WORKSPACE_POOL=${AMSS_RHS_WORKSPACE_POOL} VACUUM_BSSN=${AMSS_VACUUM_BSSN} FUSED_RHS_TAIL=${AMSS_FUSED_RHS_TAIL} BULK_SIMD=${AMSS_RHS_BULK_SIMD} BATCH_STENCIL_SINGLE_REGION=${AMSS_BATCH_STENCIL_SINGLE_REGION} ENFORCE_GA_POINTWISE=${AMSS_ENFORCE_GA_POINTWISE} RK4_FUSE_LOWERBOUND=${AMSS_RK4_FUSE_LOWERBOUND} RK4_BATCHED=${AMSS_RK4_BATCHED} SOMMERFELD_BATCHED=${AMSS_SOMMERFELD_BATCHED} SKIP_CONSTRAINT_STORES=${AMSS_RHS_SKIP_CONSTRAINT_STORES} SKIP_CHIN1_SCAN=${AMSS_RHS_SKIP_CHIN1_SCAN} RAW_DIAG_LOPSIDED=${AMSS_RHS_RAW_DIAG_LOPSIDED} FGFS_COLORED_SLAB=${AMSS_FGFS_COLORED_SLAB} FGFS_HUGEPAGE_SLAB=${AMSS_FGFS_HUGEPAGE_SLAB} TILE_J=${AMSS_RHS_TILE_J} TILE_K=${AMSS_RHS_TILE_K} SOFTWARE_PREFETCH=${AMSS_SOFTWARE_PREFETCH} BATCH_OUTBD=${AMSS_BATCH_OUTBD} LEVEL_DIAG=${AMSS_LEVEL_DIAG} CACHE_COMM_PLANS=${AMSS_CACHE_COMM_PLANS} ASYNC_LEVEL_SYNC=${AMSS_ASYNC_LEVEL_SYNC} MPI_CUDA_AWARE=${AMSS_MPI_CUDA_AWARE}")
message(STATUS "==========================================================")
