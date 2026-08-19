#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "Block.h"
#include "misc.h"

Block::Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfsi, int levi, const int cgpui) : rank(ranki), ingfs(ingfsi), fngfs(fngfsi), lev(levi), cgpu(cgpui)
{
	for (int i = 0; i < dim; i++)
		X[i] = 0;
#ifdef AMSS_FGFS_COLORED_SLAB
	fgfs_base = 0;
#endif

	if (DIM != dim)
	{
		cout << "dimension is not consistent in Block construction" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	bool flag = false;
	for (int i = 0; i < dim; i++)
	{
		shape[i] = shapei[i];
		if (shape[i] <= 0)
			flag = true;
		bbox[i] = bboxi[i];
		bbox[dim + i] = bboxi[dim + i];
	}

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (flag)
	{
		cout << "myrank: " << myrank << ", on rank: " << rank << endl;
		cout << "error shape in Block construction: (" << shape[0] << "," << shape[1] << "," << shape[2] << ")" << endl;
		cout << "box boundary: (" << bbox[0] << ":" << bbox[3] << "," << bbox[1] << ":" << bbox[4] << "," << bbox[2] << ":" << bbox[5] << ")" << endl;
		cout << "belong to level " << lev << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (myrank == rank)
	{
		for (int i = 0; i < dim; i++)
		{
			X[i] = new double[shape[i]];
			double h = (bbox[dim + i] - bbox[i]) / shape[i];
			for (int j = 0; j < shape[i]; j++)
				X[i][j] = bbox[i] + (j + 0.5) * h;
#ifdef USE_GPU
			d_X[i] = GPUManager::getInstance().allocate_device_memory(shape[i]);
			GPUManager::sync_to_gpu(X[i], d_X[i], shape[i]);
#endif
		}

		int nn = shape[0] * shape[1] * shape[2];
		fgfs = new double *[fngfs];

#ifdef USE_GPU
		// initialize GPU shadow pointers and valid flags
		d_fgfs = new double *[fngfs];
		cpu_valid = new bool[fngfs];
		gpu_valid = new bool[fngfs];
#endif

		#ifdef AMSS_FGFS_COLORED_SLAB
		// Keep allocations independent: the legacy code and Fortran runtime
		// assume each fgfs field has its own ownership.  The 32-line color
		// rotates equal point offsets through the L1D index bits.
		const size_t color_lines = 32u;
		const size_t color_doubles = color_lines * 8u;
		fgfs_base = new double *[fngfs];
		for (int i = 0; i < fngfs; i++) {
			void *raw = 0;
			if (posix_memalign(&raw, 64, sizeof(double) * (static_cast<size_t>(nn) + color_doubles)) != 0 || !raw)
			{
				cout << "on node#" << rank << ", out of memory when constructing colored Block." << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			fgfs_base[i] = static_cast<double *>(raw);
			fgfs[i] = fgfs_base[i] + (static_cast<size_t>(i) % color_lines) * 8u;
			memset(fgfs[i], 0, sizeof(double) * nn);
		}
		#else
		for (int i = 0; i < fngfs; i++) {
			fgfs[i] = (double *)malloc(sizeof(double) * nn);
			if (!(fgfs[i]))
			{
				cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			memset(fgfs[i], 0, sizeof(double) * nn);
		#endif

#ifdef USE_GPU
			d_fgfs[i] = GPUManager::getInstance().allocate_device_memory(nn);

			cpu_valid[i] = true;
			gpu_valid[i] = false;
#endif
		#ifndef AMSS_FGFS_COLORED_SLAB
		}
		#endif

		igfs = new int *[ingfs];
		for (int i = 0; i < ingfs; i++)
		{
			igfs[i] = (int *)malloc(sizeof(int) * nn);
			if (!(igfs[i]))
			{
				cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			memset(igfs[i], 0, sizeof(int) * nn);
		}
	}

#ifdef USE_GPU
	stream = GPUManager::getInstance().get_stream();
#endif
}
Block::~Block()
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == rank) {
		int nn = shape[0] * shape[1] * shape[2];
		for (int i = 0; i < dim; i ++) {
#ifdef USE_GPU
			GPUManager::getInstance().free_device_memory(d_X[i], shape[i]);
#endif
			delete[] X[i];
		}
		for (int i = 0; i < ingfs; i++)
			free(igfs[i]);
		delete[] igfs;
		#ifndef AMSS_FGFS_COLORED_SLAB
		for (int i = 0; i < fngfs; i ++) {
			free(fgfs[i]);
#ifdef USE_GPU
			GPUManager::getInstance().free_device_memory(d_fgfs[i], nn);
#endif
		}
		#else
		for (int i = 0; i < fngfs; i++)
			free(fgfs_base[i]);
		delete[] fgfs_base;
		#endif
		delete[] fgfs;
#ifdef USE_GPU
		delete[] d_fgfs;
		delete[] cpu_valid;
		delete[] gpu_valid;
#endif
		fgfs = 0;
#ifdef AMSS_FGFS_COLORED_SLAB
		fgfs_base = 0;
#endif
#ifdef USE_GPU
		d_fgfs = 0;
#endif
		X[0] = X[1] = X[2] = 0;
		igfs = 0;
	}
}
void Block::checkBlock()
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0)
	{
		cout << "belong to level " << lev << endl;
		cout << "shape: [";
		for (int i = 0; i < dim; i++)
		{
			cout << shape[i];
			if (i < dim - 1)
				cout << ",";
			else
				cout << "]";
		}
		cout << " resolution: [";
		for (int i = 0; i < dim; i++)
		{
			cout << getdX(i);
			if (i < dim - 1)
				cout << ",";
			else
				cout << "]" << endl;
		}
		cout << "locate on node " << rank << ", at (includes ghost zone):" << endl;
		cout << "(";
		for (int i = 0; i < dim; i++)
		{
			cout << bbox[i] << ":" << bbox[dim + i];
			if (i < dim - 1)
				cout << ",";
			else
				cout << ")" << endl;
		}
		cout << "has " << ingfs << " int type grids functions," << fngfs << " double type grids functions" << endl;
	}
}
double Block::getdX(int dir)
{
	if (dir < 0 || dir >= dim)
	{
		cout << "Block::getdX: error input dir = " << dir << ", this Block has direction (0," << dim - 1 << ")" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	double h;
	h = (bbox[dim + dir] - bbox[dir]) / shape[dir];
	return h;
}
void Block::swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank)
{
	if (rank == myrank)
	{
		MyList<var> *varl1 = VarList1, *varl2 = VarList2;
		while (varl1 && varl2)
		{
			misc::swap<double *>(fgfs[varl1->data->sgfn], fgfs[varl2->data->sgfn]);
#ifdef AMSS_FGFS_COLORED_SLAB
			misc::swap<double *>(fgfs_base[varl1->data->sgfn], fgfs_base[varl2->data->sgfn]);
#endif
#ifdef USE_GPU
			misc::swap<double *>(d_fgfs[varl1->data->sgfn], d_fgfs[varl2->data->sgfn]);
			misc::swap<bool>(cpu_valid[varl1->data->sgfn], cpu_valid[varl2->data->sgfn]);
			misc::swap<bool>(gpu_valid[varl1->data->sgfn], gpu_valid[varl2->data->sgfn]);
#endif
			varl1 = varl1->next;
			varl2 = varl2->next;
		}
		if (varl1 || varl2)
		{
			cout << "error in Block::swaplist, var lists does not match." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
}

#ifdef USE_GPU
void Block::require_on_gpu(int var_index) {
	if (!gpu_valid[var_index] && cpu_valid[var_index]) {
		int nn = shape[0] * shape[1] * shape[2];
		GPUManager::sync_to_gpu(fgfs[var_index], d_fgfs[var_index], nn);
		gpu_valid[var_index] = true;
	}
}

void Block::require_on_cpu(int var_index) {
	if (!cpu_valid[var_index] && gpu_valid[var_index]) {
		int nn = shape[0] * shape[1] * shape[2];
		GPUManager::sync_to_cpu(fgfs[var_index], d_fgfs[var_index], nn);
		cpu_valid[var_index] = true;
	}
}

void Block::mark_gpu_modified(int var_index) {
	gpu_valid[var_index] = true;
	cpu_valid[var_index] = false;
}

void Block::mark_cpu_modified(int var_index) {
	cpu_valid[var_index] = true;
	gpu_valid[var_index] = false;
}

void Block::move_to_gpu(MyList<var> *VarList) {
	MyList<var> *iter = VarList;
	while (iter) {
		int var_index = iter->data->sgfn;
		mark_cpu_modified(var_index);
		require_on_gpu(var_index);
		iter = iter->next;
	}
}

void Block::move_to_cpu(MyList<var> *VarList) {
	MyList<var> *iter = VarList;
	while (iter) {
		int var_index = iter->data->sgfn;
		mark_gpu_modified(var_index);
		require_on_cpu(var_index);
		iter = iter->next;
	}
}
#endif
