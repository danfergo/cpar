#include <cuda.h>
#include <device_functions.h>
#include <cuda_runtime_api.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <math.h>
#include <unistd.h>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
#define TILE_WIDTH	32
#define SYSTEMTIME struct timespec

#ifdef __INTELLISENSE__
void __syncthreads();
#endif


inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d \n", cudaGetErrorString(code), file, line);

		if (abort){
			system("pause");
			exit(code);
		}
	}
}

using namespace std;

// ------------------------------------ INIT, PRINT & FREE MATRICES; PRINT HEADERS & RESULTS; MAIN LOOP;
void initMatrices(double ** pha, double ** phb, double ** phc, int m_ar, int m_br){
	long int i, j;

	*pha = (double *)malloc((m_ar * m_br) * sizeof(double));
	*phb = (double *)malloc((m_ar * m_br) * sizeof(double));
	*phc = (double *)malloc((m_ar * m_ar) * sizeof(double));


	for (i = 0; i< m_ar; i++)
	for (j = 0; j< m_br; j++) {
		(*pha)[i*m_br + j] = i == j ? (double)i : (double)0;
	}


	for (i = 0; i< m_br; i++)
	for (j = 0; j< m_ar; j++) {
		(*phb)[i*m_ar + j] = i == j ? (double)i : (double)0;
	}

	for (i = 0; i< m_ar; i++)
	for (j = 0; j< m_ar; j++) {
		(*phc)[i*m_ar + j] = 0;
	}

}

void freeMatrices(double * pha, double * phb, double * phc){
	free(pha);
	free(phb);
	free(phc);
}


void printResultMatrix(double ** phc, int m_br) {
	int i, j;

	cout << endl << endl;

	cout << "Result matrix: " << endl;
	for (i = 0; i<10; i++)
	{

		for (j = 0; j<min(10, m_br); j++)
			cout << (*phc)[i*m_br + j] << " ";
		cout << "\n";

	}
	cout << endl << endl;
}

// ------------------------------------ CUDA HEADERS
int getSPcores(cudaDeviceProp devProp)
{
	int cores = 0;
	int mp = devProp.multiProcessorCount;
	switch (devProp.major){
	case 2: // Fermi
		if (devProp.minor == 1) cores = mp * 48;
		else cores = mp * 32;
		break;
	case 3: // Kepler
		cores = mp * 192;
		break;
	case 5: // Maxwell
		cores = mp * 128;
		break;
	default:
		printf("Unknown device type\n");
		break;
	}
	return cores;
}

void DisplayHeader()
{
	const int kb = 1024;
	const int mb = kb * kb;
	cout << "GPU" << endl << "=========" << endl;

	cout << "CUDA version:   v" << CUDART_VERSION << endl;
	//wcout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << endl << endl;

	int devCount;
	cudaGetDeviceCount(&devCount);
	cout << "CUDA Devices: " << endl << endl;

	for    (int i = 0; i < devCount; i++)
	{
		cudaDeviceProp props;
		cudaGetDeviceProperties(&props, i);
		cout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
		cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
		cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
		cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
		cout << "  Block registers: " << props.regsPerBlock << endl << endl;
		cout << "  Warp size:         " << props.warpSize << endl;
		cout << "  Number of cores: " << getSPcores(props) << endl;
		cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
		cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << endl;
		cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << endl;
		cout << endl;
	}

	cout << "=========" << endl;
}



// ------------------------------------ STAT FUNCTIONS

__global__ void MatrixMulKernelNaive(double * Md, double * Nd, double * Pd, int Width){
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int Y = by*TILE_WIDTH + ty;
	int X = bx*TILE_WIDTH + tx;
	double Pvalue = 0;
	for(int m = 0; m < Width; m++){
		Pvalue += Md[X*Width + m]* Nd[m*Width + X];
	}


	Pd[Y*Width + X] = Pvalue;
}

__global__ void MatrixMulKernel(double * Md, double * Nd, double * Pd, /*double * Td,*/ int Width){
	__shared__ float Mds[TILE_WIDTH][TILE_WIDTH];
	__shared__ float Nds[TILE_WIDTH][TILE_WIDTH];

	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int Y = by*TILE_WIDTH + ty;
	int X = bx*TILE_WIDTH + tx;
	double Pvalue = 0;

	for (int m = 0; m < Width/TILE_WIDTH; ++m){
		Mds[tx][ty] = Md[(m*TILE_WIDTH + tx)*Width + Y];
		Nds[tx][ty] = Nd[X*Width + (m*TILE_WIDTH + ty)];
		__syncthreads();

		for (int k = 0; k < TILE_WIDTH; ++k)
			Pvalue +=  Mds[tx][k] * Nds[k][ty];
		__syncthreads();

	}

	Pd[Y*Width + X] = Pvalue;
        __syncthreads();
        if(tx == 0 && ty == 0){
            
        }
}


// ------------------------------------ CUDA OPERATIONS PREPARATION
void MatrixMulOnDevice(double * M, double * N, double * P, int Width){
	long int size = Width * Width * sizeof(double);
	double * Md, *Nd, *Pd;// *Td;
	dim3 dimGrid(Width / TILE_WIDTH, Width / TILE_WIDTH);
	dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);
        
        
	// Allocate and Load M, N to device memory
	gpuErrchk(cudaMalloc(&Md, size));
	gpuErrchk(cudaMemcpy(Md, M, size, cudaMemcpyHostToDevice));

	gpuErrchk(cudaMalloc(&Nd, size));
	gpuErrchk(cudaMemcpy(Nd, N, size, cudaMemcpyHostToDevice));

	// Allocate P on the device
	gpuErrchk(cudaMalloc(&Pd, size));

//        gpuErrchk(cudaMalloc(&Td, pow((float)Width / TILE_WIDTH, 2) );

        
        
	SYSTEMTIME time1;
	if(clock_gettime(CLOCK_MONOTONIC,&time1) < 0) {
		perror("Failure getting process time (1)");
		exit(EXIT_FAILURE);
	}





	// Kernel invocation
	MatrixMulKernel << <dimGrid, dimBlock >> >(Md, Nd, Pd, /*Td,*/ Width);
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());




	SYSTEMTIME time2;

	if(clock_gettime(CLOCK_MONOTONIC,&time2) < 0) {
			perror("Failure getting process time (2)");
			exit(EXIT_FAILURE);
	}
	double delta = (double)(time2.tv_sec - time1.tv_sec) + ((double)(time2.tv_nsec - time1.tv_nsec) / pow(10.f,9));

	long long calcNrInsOp = ((long long)3*pow((float)Width,3));
	double tPerformance = calcNrInsOp/(double)(delta* pow(10.f,9));

	printf("Time  : %3.3f seconds\n", delta);
	printf("T.Prfr: %3.3f GFLOPS\n", tPerformance);

	// Read P from the device
	gpuErrchk(cudaMemcpy(P, Pd, size, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaDeviceSynchronize());

	// Free device matrices
	cudaFree(Md);
	cudaFree(Nd);
	cudaFree(Pd);
        //cudaFree(Td);
}


int main(int argc, char *argv[])
{
	DisplayHeader();

	int Width = argc > 1 ? atoi(argv[1]) : 1024;
	long int size = Width * Width * sizeof(double);
	cout << "Width: " << Width << endl;
	cout << "Global mem: " << (size * 3) / (1024 * 1024) << "mb" << endl;
	cout << "Shared mem: " << (TILE_WIDTH*TILE_WIDTH*sizeof(double) * 2) / (1024) << "kb" << endl;

	cout << "dim grid:  " << (Width / TILE_WIDTH) << ", "  << (Width / TILE_WIDTH) << endl;
	cout << "dim block: " << TILE_WIDTH << ", " << TILE_WIDTH << endl;


	double * M, *N, *P;

	initMatrices(&M, &N, &P, Width, Width);
	MatrixMulOnDevice(M, N, P, Width);
	printResultMatrix(&P, Width);

	freeMatrices(N, M, P);

	cudaDeviceReset();
	return 0;
}