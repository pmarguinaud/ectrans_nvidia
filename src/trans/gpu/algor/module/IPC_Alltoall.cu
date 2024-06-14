// AlltoAll for NVLINK-connected GPUs ithin a single server, using CUDA IPC
// All pairs of (distinct) GPUs must return 1 for canAccessPeer
// Alan Gray, NVIDIA

#include <stdio.h>
#include <assert.h>

#include <stdlib.h>
#include "cuda_profiler_api.h"
#include <cuda_runtime_api.h>

#include <unistd.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <linux/version.h>


// maximum number of devices supported
#define MAX_DEVICES          32

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

// data structure required for IPC setup 

typedef struct ipcDevices_st
{
  int count;
  int ordinals[MAX_DEVICES];
} ipcDevices_t;



// we need a seperate CUDA stream for each target GPU
static cudaStream_t streams[MAX_DEVICES];

// structure to contain pointers to remote array data, and offsets into each for destinatin data.
// we maintain 2 copies (first array dimension), corresponding to MTOL and LTOM trans comms.
// This allows us to only perform setup steps the first time, and re-use all following times. 
static double* outputptrall[2][MAX_DEVICES];
static int roff_remote[2][MAX_DEVICES];


// Initialize IP communicatins
extern "C" void initIPC(double* output_d,int* roff, int mtol_or_ltom){

}


static bool already_initialized[2]={0,0};

static bool notFullPeerAccess=0;

// main externally visible routine for performing AlltoAll comms
extern "C" int Alltoallv_CUDAIPC(double* input, int* len, int* soff, 
				  double* output, int* roff,int mtol_or_ltom){


  return 0; 
  
  
}


