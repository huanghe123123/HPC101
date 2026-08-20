
#ifndef BLOCK_H
#define BLOCK_H

#include <cstddef>
#include <mpi.h>
#include "macrodef.h" //need dim here; Vertex or Cell
#include "var.h"
#include "MyList.h"

#ifdef USE_GPU
#include "gpu_manager.h"   // brings in cudaStream_t / cuda_runtime.h
#endif

class Block
{

public:
   int shape[dim];
   double bbox[2 * dim];
   double *X[dim];
   int rank; // where the real data locate in
   int lev, cgpu;
   int ingfs, fngfs;
   int *(*igfs);
   double *(*fgfs); // fine grid functions
#ifdef AMSS_FGFS_COLORED_SLAB
   // Raw allocation bases follow fgfs pointers through swapList.  Each field
   // view is shifted by a small cache-line color while retaining nn
   // contiguous elements and an independently freeable allocation.
   double *(*fgfs_base);
#endif
#ifdef AMSS_FGFS_HUGEPAGE_SLAB
   // Single contiguous 2 MiB-aligned, MADV_HUGEPAGE-backed allocation holding
   // all fngfs fields.  fgfs[i] are view pointers (fgfs_slab + i*nn); swapList
   // only swaps the views, the slab owns the memory and is freed once in ~Block.
   double *fgfs_slab;
#endif

#ifdef USE_GPU
   // GPU Shadow pointers and valid flags
   double *d_X[dim];
   double *(*d_fgfs);
   bool *cpu_valid;
   bool *gpu_valid;
   cudaStream_t stream;
#endif

public:
   Block() {};
   Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfs, int levi, const int cgpui = 0);

   ~Block();

   void checkBlock();

   double getdX(int dir);
   void swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank);

#ifdef USE_GPU
   void require_on_gpu(int var_index);
   void require_on_cpu(int var_index);
   void mark_gpu_modified(int var_index);
   void mark_cpu_modified(int var_index);

   void move_to_gpu(MyList<var> *VarList);
   void move_to_cpu(MyList<var> *VarList);
#endif
};

#endif /* BLOCK_H */
