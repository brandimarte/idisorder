/*  Put license and stuff here
 *
 *
 * Warning to the navigators:
 *
 * Do not dare to read this code if you are not familiar 
 * with the BLAS/LAPACK libraries and CUDA (and MAGMA).
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include "timer.h"

#ifdef MAGMA
# include <cuda.h>
# include <cuda_runtime_api.h>
# include <cublas_v2.h>
# include <magma.h>
# define doubleComplex cuDoubleComplex
#else
struct complex
{
    double Re, Im;

    complex () : Re(0.0), Im(0.0) { }
    complex (double r, double i) : Re(r), Im(i) { }

    complex operator+(const complex &C) const { complex t = *this; return t += C; }
    complex operator-(const complex &C) const { complex t = *this; return t -= C; }
    complex operator-() const { return complex( -Re, -Im); }
    complex operator*(const complex &C) const { complex t = *this; return t *= C; }
    complex& operator=(const complex &C)  { Re = C.Re;  Im = C.Im;  return *this; }
    complex& operator+=(const complex &C) { Re += C.Re; Im += C.Im; return *this; }
    complex& operator-=(const complex &C) { Re -= C.Re; Im -= C.Im; return *this; }
    complex& operator*=(const complex &C)
    {
        complex t = *this;
        Re = t.Re*C.Re - t.Im*C.Im;
        Im = t.Im*C.Re + t.Re*C.Im;
        return *this;
    }
    double mod() const { return Re*Re + Im*Im; }
    complex conj() const { return complex(Re, -Im); }
};
typedef complex doubleComplex;
#endif

#define IDX( i, j, LD ) ((i)+(j)*(LD))


// Constants to call clock* functions (which take args. by reference)
// static const int Clock_id_MatInversion     = CLOCK_ID_MatInversion;
// static const int Clock_id_MatMult          = CLOCK_ID_MatMult;
// #ifdef MAGMA
// static const int Clock_id_MatInversion_GPU = CLOCK_ID_gpuMatInversion;
// static const int Clock_id_MatMult_GPU      = CLOCK_ID_gpuMatMult;
// #endif


//----------------------------------------------------------------
// External functions used here
// We use the Fortran versions (the C ones are not stardard)
extern "C" void zgemm_(char *opA, char *opB, int *M, int *N, int *K, doubleComplex *alpha, doubleComplex *A, int *LDA, doubleComplex *B, int *LDB, doubleComplex *beta, doubleComplex *C, int *LDC);
extern "C" void zsymm_(char *side, char *uplo, int *M, int *N, doubleComplex *alpha, doubleComplex *A, int *LDA, doubleComplex *B, int *LDB, doubleComplex *beta, doubleComplex *C, int *LDC);
extern "C" void zhemm_(char *side, char *uplo, int *M, int *N, doubleComplex *alpha, doubleComplex *A, int *LDA, doubleComplex *B, int *LDB, doubleComplex *beta, doubleComplex *C, int *LDC);
extern "C" void zgetrf_(int *M, int *N, doubleComplex *A, int *LDA, int *ipiv, int *info);
extern "C" void zgetri_(int *N, doubleComplex *A, int *LDA, int *ipiv, doubleComplex *work, const int *lwork, int *info);
extern "C" void zsytrf_(char *uplo, int *N, doubleComplex *A, int *LDA, int *ipiv, doubleComplex *work, const int *lwork, int *info);
extern "C" void zsytri_(char* uplo, int *N, doubleComplex *A, int *LDA, int *ipiv, doubleComplex *work, int *info);


/* For use with fortran (mandatory, of course): */
# define HI_Init             hi_init_
# define HI_Finalize         hi_finalize_
# define HI_GetGPUcount      hi_getgpucount_
# define HI_PrintInfo        hi_printinfo_
# define HI_Info             hi_info_

# define HI_zgemm            hi_zgemm_
# define HI_zsymm            hi_zsymm_
# define HI_zhemm            hi_zhemm_
# define HI_zgeInvert        hi_zgeinvert_
# define HI_zsyInvert        hi_zsyinvert_
# define HI_Invert_prepare   hi_invert_prepare_
# define HI_Invert_expert    hi_invert_expert_
# define HI_Invert_clear     hi_invert_clear_


/* Function prototypes */
extern "C" void HI_Init(int *Proc_id, int *ProcsPerGPU);
extern "C" void HI_Finalize(void);
extern "C" void HI_GetGPUcount(int *Ngpus);
extern "C" void HI_PrintInfo(int *Proc_id, int *ProcsPerGPU);
extern "C" void HI_Info(int *Proc_id, int *ProcsPerGPU, int *fmyGPU, int *fGPUcount, int *fVirtualGPUcount, int *fIuseGPU);

#ifdef MAGMA
static inline cublasOperation_t char2cublas_op  (char *char_op);
static inline cublasSideMode_t  char2cublas_side(char *char_op);
static inline cublasFillMode_t  char2cublas_fill(char *char_op);
#endif

extern "C" void HI_zgemm(char *opA, char *opB, int *M, int *N, int *K,
                         doubleComplex *alpha, doubleComplex *A, int *LDA,
                         doubleComplex *B, int *LDB, doubleComplex *beta,
                         doubleComplex *C, int *LDC);
extern "C" void HI_zsymm(char *side, char *uplo, int *M, int *N,
                         doubleComplex *alpha, doubleComplex *A, int *LDA,
                         doubleComplex *B, int *LDB, doubleComplex *beta,
                         doubleComplex *C, int *LDC);
extern "C" void HI_zhemm(char *side, char *uplo, int *M, int *N,
                         doubleComplex *alpha, doubleComplex *A, int *LDA,
                         doubleComplex *B, int *LDB, doubleComplex *beta,
                         doubleComplex *C, int *LDC);

extern "C" void HI_zgeInvert(doubleComplex *A, int *M);
extern "C" void HI_zsyInvert(char *uplo, doubleComplex *A, int *M);

extern "C" void HI_Invert_prepare(int **ipiv, int N, doubleComplex **Work, int Nwork);
extern "C" void HI_Invert_expert(doubleComplex *SE, int N, int *ipiv, doubleComplex *Work, int Nwork, int *info);
extern "C" void HI_Invert_clear(int *ipiv, doubleComplex *Work);


template <typename T>
static inline const T max(T a, T b) { return( (a>b) ? a : b ); }
template <typename T>
static inline const T min(T a, T b) { return( (a<b) ? a : b ); }

// Globlal variables
int GPUcount = 0;
int VirtualGPUcount = 0;
int myGPU = -1;
bool IuseGPU = false;
#ifdef MAGMA
cublasHandle_t myHandle;
#endif

//-----------------------------------------
// Initializes the environment
//
// Alberto Torres, Jan 2014
//-----------------------------------------
// Inputs:
// int *Proc_id : Number of the process, rank, pid
// int *ProcsPerGPU : Number of processes to run per GPU
//-----------------------------------------
void HI_Init(int *Proc_id, int *ProcsPerGPU) //, int *FillMode)
{
#ifdef MAGMA
    cudaGetDeviceCount(&GPUcount);
	if(GPUcount == 0 && *Proc_id == 0)
	{
		printf("\nWARNING: No GPUs found. All calculations will be done in CPU.\n");
		fflush(stdout);
	}

	if(*ProcsPerGPU < 1) *ProcsPerGPU = 0;
    VirtualGPUcount = GPUcount * (*ProcsPerGPU);   // We can pin more than one process to a gpu with ProcsPerGPU > 1

    if(*Proc_id < VirtualGPUcount)
    {
        IuseGPU = true;
        magma_init();

        myGPU = (*Proc_id / *ProcsPerGPU) % GPUcount;
        cudaSetDevice(myGPU);   // Bind to GPUs in a round-robin fashion according to Proc_id
        cublasCreate( &myHandle );
        //cublasSetStream(myHandle, NULL);   // NULL is the default stream
    }
#endif
}


//-------------------------------------------------------------------
// Print information about the GPUs and wich processes will use them.
//
void HI_PrintInfo(int *Proc_id, int *ProcsPerGPU)
{
    if(*Proc_id == 0)
    {
        printf("\nUsing: (%i processes per GPU)\n", *ProcsPerGPU);
#ifdef MAGMA
        magma_print_devices();
        printf("\n");
#endif
        fflush(stdout);
    }

    if(IuseGPU)
        printf("Process nr. %i using GPU nr. %i of %i(%i)\n", *Proc_id, myGPU, GPUcount, VirtualGPUcount);
    else
        printf("Process nr. %i using CPU.\n", *Proc_id);

    fflush(stdout);
}

//-------------------------------------------------------------------
// Print and returns information about the GPUs and wich processes
// will use them.
//
void HI_Info (int *Proc_id, int *ProcsPerGPU, int *fmyGPU, int *fGPUcount, int *fVirtualGPUcount, int *fIuseGPU)
{
    if(*Proc_id == 0)
    {
        printf("procInfo: Using %i processes per GPU\n", *ProcsPerGPU);
#ifdef MAGMA
        printf("procInfo: Magma info:\n\n");
        magma_print_devices();
        printf("\n");
#endif
        fflush(stdout);
    }

    *fmyGPU = myGPU;
    *fGPUcount = GPUcount;
    *fVirtualGPUcount = VirtualGPUcount;
    if(IuseGPU)
       *fIuseGPU = 1;
    else
       *fIuseGPU = 0;

}

//-----------------------------------------
// Finalizes the environment
//
void HI_Finalize(void)
{
#ifdef MAGMA
    if(IuseGPU)
    {
        cublasDestroy(myHandle);
        magma_finalize();
    }
#endif
}


//-----------------------------------------
// Returns the number of GPUs
void HI_GetGPUcount(int *Ngpus) { *Ngpus = GPUcount; }


//-------------------------------------------------------------------
// Returns cublas operator corresponding to char_op (fortran like):
// N or n --> CUBLAS_OP_N
// T or t --> CUBLAS_OP_T
// C or c --> CUBLAS_OP_C
#ifdef MAGMA
static inline cublasOperation_t char2cublas_op(char *char_op)
{
    if      (*char_op == 'N' || *char_op == 'n') return(CUBLAS_OP_N);
    else if (*char_op == 'T' || *char_op == 't') return(CUBLAS_OP_T);
    else return(CUBLAS_OP_C);
    // else if (*char_op == 'C' || *char_op == 'c') return(CUBLAS_OP_C);
}
#endif


//-------------------------------------------------------------------
// Multiplies two double complex matrices using only 1 xPU (x=G,C).
//
// Same interface as BLAS 3 zgemm
//-------------------------------------------------------------------
// Alberto Torres, Jan 2014
//-------------------------------------------------------------------
void HI_zgemm(char *opA, char *opB, int *M, int *N, int *K, 
              doubleComplex *alpha, doubleComplex *A, int *LDA,
              doubleComplex *B, int *LDB, doubleComplex *beta,
              doubleComplex *C, int *LDC)
{
#ifdef MAGMA
    if(IuseGPU && *N>100) // Use GPU
    {
        // clock_start(&Clock_id_MatMult_GPU);

        const size_t size = sizeof(cuDoubleComplex);

        const int ldA = *LDA;
        const int ldB = *LDB;
        const int ldC = *LDC;
        const int lddA = ldA;  // This ones could be made multiples of 32
        const int lddB = ldB;
        const int lddC = ldC;

        const cuDoubleComplex za=(cuDoubleComplex)*alpha;
        const cuDoubleComplex zb=(cuDoubleComplex)*beta;

        cublasOperation_t cuOpA = char2cublas_op(opA);
        cublasOperation_t cuOpB = char2cublas_op(opB);

        int rA, cA, rB, cB;
        if(cuOpA == CUBLAS_OP_N) { rA=*M; cA=*K; }   // Rows and columns of A, *not* A^(opA).
        else                     { rA=*K; cA=*M; }

        if(cuOpB == CUBLAS_OP_N) { rB=*K; cB=*N; }   // Rows and columns of B, *not* B^(opB).
        else                     { rB=*N; cB=*K; }

        const int rC = *M;
        const int cC = *N;

        cuDoubleComplex *dA, *dB, *dC;
       
        cudaMalloc((void**)&dA, lddA*cA * size);
        cudaMalloc((void**)&dB, lddB*cB * size);
        cudaMalloc((void**)&dC, lddC*cC * size);

        cublasSetMatrix(rA, cA, size, A, ldA, dA, lddA);
        cublasSetMatrix(rB, cB, size, B, ldB, dB, lddB);

        cublasZgemm(myHandle, cuOpA, cuOpB, *M, *N, *K, &za, dA, lddA, dB, lddB, &zb, dC, lddC);
 
        cublasGetMatrix(rC, cC, size, dC, lddC, C, ldC);

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);

        // clock_stop(&Clock_id_MatMult_GPU);
    }
    else // Use CPU
#endif
    {
        // clock_start(&Clock_id_MatMult);
        zgemm_(opA, opB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
        // clock_stop(&Clock_id_MatMult);
    }
}


//-------------------------------------------------------------------
// Returns cublas operator corresponding to char_op (fortran like):
// L or l --> CUBLAS_SIDE_LEFT
// R or r --> CUBLAS_SIDE_RIGHT
#ifdef MAGMA
static inline cublasSideMode_t char2cublas_side(char *char_op)
{
   if (*char_op == 'L' || *char_op == 'l') return(CUBLAS_SIDE_LEFT);
   else return(CUBLAS_SIDE_RIGHT);
   // else if (*char_op == 'R' || *char_op == 'r') return(CUBLAS_SIDE_RIGHT);
}
#endif


//-------------------------------------------------------------------
// Returns cublas operator corresponding to char_op (fortran like):
// U or u --> CUBLAS_FILL_MODE_UPPER
// L or l --> CUBLAS_FILL_MODE_LOWER
#ifdef MAGMA
static inline cublasFillMode_t char2cublas_fill(char *char_op)
{
   if (*char_op == 'L' || *char_op == 'l') return(CUBLAS_FILL_MODE_UPPER);
   else return(CUBLAS_FILL_MODE_LOWER);
   // else if (*char_op == 'U' || *char_op == 'u') return(CUBLAS_FILL_MODE_LOWER);
}
#endif


//----------------------------------------------------------------------------
// Multiplies two double complex symmetric matrices using only 1 xPU (x=G,C).
//
// Same interface as BLAS 3 zsymm
//----------------------------------------------------------------------------
// Alberto Torres, Jan 2014
//----------------------------------------------------------------------------
void HI_zsymm(char *side, char *uplo, int *M, int *N, 
              doubleComplex *alpha, doubleComplex *A, int *LDA,
              doubleComplex *B, int *LDB, doubleComplex *beta,
              doubleComplex *C, int *LDC)
{
#ifdef MAGMA
    if(IuseGPU && *N>100) // Use GPU
    {
        // clock_start(&Clock_id_MatMult_GPU);

        const size_t size = sizeof(cuDoubleComplex);

        const int ldA = *LDA;
        const int ldB = *LDB;
        const int ldC = *LDC;
        const int lddA = ldA;  // This ones could be made multiples of 32
        const int lddB = ldB;
        const int lddC = ldC;

        const cuDoubleComplex za=(cuDoubleComplex)*alpha;
        const cuDoubleComplex zb=(cuDoubleComplex)*beta;

        cublasSideMode_t cuSide = char2cublas_side(side);
        cublasFillMode_t cuUpLo = char2cublas_fill(uplo);

        int rA;
        if(cuSide == CUBLAS_SIDE_LEFT) rA = *M;
        else rA = *N;
        const int cA = rA;
        const int rB = *M, cB = *N;
        const int rC = *M, cC = *N;

        cuDoubleComplex *dA, *dB, *dC;
       
        cudaMalloc((void**)&dA, lddA*cA * size);
        cudaMalloc((void**)&dB, lddB*cB * size);
        cudaMalloc((void**)&dC, lddC*cC * size);

        cublasSetMatrix(rA, cA, size, A, ldA, dA, lddA);
        cublasSetMatrix(rB, cB, size, B, ldB, dB, lddB);

        cublasZsymm(myHandle, cuSide, cuUpLo, rC, cC, &za, dA, lddA, dB, lddB, &zb, dC, lddC);

        cublasGetMatrix(rC, cC, size, dC, lddC, C, ldC);

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);

        // clock_stop(&Clock_id_MatMult_GPU);
    }
    else // Use CPU
#endif
    {
        // clock_start(&Clock_id_MatMult);
        zsymm_(side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
        // clock_stop(&Clock_id_MatMult);
    }
}


//----------------------------------------------------------------------------
// Multiplies two double complex hermitian matrices using only 1 xPU (x=G,C).
//
// Same interface as BLAS 3 zhemm
//----------------------------------------------------------------------------
// Alberto Torres, Jan 2014
//----------------------------------------------------------------------------
void HI_zhemm(char *side, char *uplo, int *M, int *N, 
              doubleComplex *alpha, doubleComplex *A, int *LDA,
              doubleComplex *B, int *LDB, doubleComplex *beta,
              doubleComplex *C, int *LDC)
{
#ifdef MAGMA
    if(IuseGPU  && *N>100) // Use GPU
    {
        // clock_start(&Clock_id_MatMult_GPU);

        const size_t size = sizeof(cuDoubleComplex);

        const int ldA = *LDA;
        const int ldB = *LDB;
        const int ldC = *LDC;
        const int lddA = ldA;  // This ones could be made multiples of 32
        const int lddB = ldB;
        const int lddC = ldC;

        const cuDoubleComplex za=(cuDoubleComplex)*alpha;
        const cuDoubleComplex zb=(cuDoubleComplex)*beta;

        cublasSideMode_t cuSide = char2cublas_side(side);
        cublasFillMode_t cuUpLo = char2cublas_fill(uplo);

        int rA;
        if(cuSide == CUBLAS_SIDE_LEFT) rA = *M;
        else rA = *N;
        const int cA = rA;
        const int rB = *M, cB = *N;
        const int rC = *M, cC = *N;

        cuDoubleComplex *dA, *dB, *dC;
       
        cudaMalloc((void**)&dA, lddA*cA * size);
        cudaMalloc((void**)&dB, lddB*cB * size);
        cudaMalloc((void**)&dC, lddC*cC * size);

        cublasSetMatrix(rA, cA, size, A, ldA, dA, lddA);
        cublasSetMatrix(rB, cB, size, B, ldB, dB, lddB);

        cublasZhemm(myHandle, cuSide, cuUpLo, rC, cC, &za, dA, lddA, dB, lddB, &zb, dC, lddC);

        cublasGetMatrix(rC, cC, size, dC, lddC, C, ldC);

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);

        // clock_stop(&Clock_id_MatMult_GPU);
    }
    else // Use CPU
#endif
    {
        // clock_start(&Clock_id_MatMult);
        zhemm_(side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
        // clock_stop(&Clock_id_MatMult);
    }
}


//----------------------------------------------------------------------------
// Inverts a double complex matrix using 1 xPU (x=C or G).
//----------------------------------------------------------------------------
// Alberto Torres, Jan 2014
//----------------------------------------------------------------------------
void HI_zgeInvert(doubleComplex *A, int *N)
{
    const int n=*N;
    const int lwork=64*n;
    int info;
    int *ipiv=(int*)malloc(n*sizeof(int));

#ifdef MAGMA
    if(IuseGPU  && *N>100)
    {
        // clock_start(&Clock_id_MatInversion_GPU);

        const int lddA=((n+15)/16)*16;  // Make lddA multiple of 32 (faster)
        cuDoubleComplex *dA, *dWork;
        
        cudaMalloc((void**)&dA,   n*lddA*sizeof(cuDoubleComplex));
        cudaMalloc((void**)&dWork, lwork*sizeof(cuDoubleComplex));

        cublasSetMatrix(n, n, sizeof(cuDoubleComplex), A, n, dA, lddA);

        magmablasSetKernelStream(0);

        magma_zgetrf_gpu(n, n, dA, lddA, ipiv, &info);
        magma_zgetri_gpu(n, dA, lddA, ipiv, dWork, lwork, &info);

        cublasGetMatrix(n, n, sizeof(cuDoubleComplex), dA, lddA, A, n);

        cudaFree((void*)dA);   
        cudaFree((void*)dWork);

        // clock_stop(&Clock_id_MatInversion_GPU);
    }
    else
#endif
    {
        // clock_start(&Clock_id_MatInversion);

        doubleComplex *Work = (doubleComplex*)malloc(lwork*sizeof(doubleComplex));

        zgetrf_(N, N, A, N, ipiv, &info);
        zgetri_(N, A, N, ipiv, Work, &lwork, &info);

        free(Work);

        // clock_stop(&Clock_id_MatInversion);
    }
    free(ipiv);
}


#define CACHE_L1 65536
#define BLOCK_SIZE 16
static inline void Zsymmetrize(char UpLo, doubleComplex *A, const int N)
{
//    const int Block_size = (sqrt((int)CACHE_L1))/sizeof(doubleComplex);

    if(UpLo == 'L' || UpLo == 'l')
    {
        for(int i=0; i<N; i++)
            for(int j=0; j<i; j++)
                A[IDX(j,i,N)] = A[IDX(i,j,N)];
    }
    else if(UpLo == 'U' || UpLo == 'u')
    {
        for(int i=0; i<N; i++)
            for(int j=0; j<i; j++)
                A[IDX(i,j,N)] = A[IDX(j,i,N)];
    }
}


//----------------------------------------------------------------------------
// Inverts a symmetric double complex matrix using 1 xPU (x=C or G).
//----------------------------------------------------------------------------
// Alberto Torres, Jan 2014
//----------------------------------------------------------------------------
void HI_zsyInvert(char *uplo, doubleComplex *A, int *N)
{
    const int n=*N;
    const int lwork=64*n;
    int info;
    int *ipiv=(int*)malloc(n*sizeof(int));

#ifdef MAGMA
    if(IuseGPU && *N>100)
    {
        // clock_start(&Clock_id_MatInversion_GPU);

        const int lddA=((n+31)/32)*32;  // Make lddA multiple of 32 (faster)
        cuDoubleComplex *dA, *dWork;
        
        cudaMalloc((void**)&dA,   n*lddA*sizeof(cuDoubleComplex));
        cudaMalloc((void**)&dWork, lwork*sizeof(cuDoubleComplex));

        cublasSetMatrix(n, n, sizeof(cuDoubleComplex), A, n, dA, lddA);

        magmablasSetKernelStream(0);
/*
 *      No zsytrf/i in MAGMA (yet?)
 *
        magma_zsytrf_gpu();
        magma_zsytri_gpu();
        magmablas_zsymmetrize(*uplo, n, dA, lddA );
 *
 *      So... calling zgetrf/i !
*/
        magma_zgetrf_gpu(n, n, dA, lddA, ipiv, &info);
        magma_zgetri_gpu(n, dA, lddA, ipiv, dWork, lwork, &info);

        cublasGetMatrix(n, n, sizeof(cuDoubleComplex), dA, lddA, A, n);

        cudaFree((void*)dA);   
        cudaFree((void*)dWork);

        // clock_stop(&Clock_id_MatInversion_GPU);
    }
    else
#endif
    {
        // clock_start(&Clock_id_MatInversion);

        doubleComplex *Work = (doubleComplex*)malloc(lwork*sizeof(doubleComplex));

        zsytrf_(uplo, N, A, N, ipiv, Work, &lwork, &info);
        zsytri_(uplo, N, A, N, ipiv, Work, &info);

        Zsymmetrize(*uplo, A, n);

        free(Work);

        // clock_stop(&Clock_id_MatInversion);
    }
    free(ipiv);
}


//-------------------------------------------------------------------
// Allocates buffers for MI_Invert_1_cpu_expert
//
// Alberto Torres
void HI_Invert_prepare(int **ipiv, int N, doubleComplex **Work, int Nwork)
{
    *ipiv = (int*)malloc(N*sizeof(int));
    *Work = (doubleComplex*)malloc(Nwork*sizeof(doubleComplex));
}

//-------------------------------------------------------------------
// Inverts a double complex matrix using 1 CPU.
// Requires buffers to be allocated already
//
// Alberto Torres
void HI_Invert_expert(doubleComplex *A, int N, int *ipiv, doubleComplex *Work, int Nwork, int *info)
{
    zgetrf_(&N, &N, A, &N, ipiv, info);
    zgetri_(&N, A, &N, ipiv, Work, &Nwork, info);
}

//-------------------------------------------------------------------
// Deallocates buffers for MI_Invert_1_cpu_expert
//
// Alberto Torres
void MI_Invert_clear(int *ipiv, doubleComplex *Work)
{
    free(ipiv);
    free(Work);
}
