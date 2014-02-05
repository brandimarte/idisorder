#ifdef MAGMA
#define CLOCK_N_ROUTINES 5
#else
#define CLOCK_N_ROUTINES 3
#endif

#define CLOCK_ID_idisorder         0
#define CLOCK_ID_MatInversion      1
#define CLOCK_ID_MatMult           2
#ifdef MAGMA
#define CLOCK_ID_gpuMatInversion   3
#define CLOCK_ID_gpuMatMult        4
#endif

#define CLOCK_NAME_idisorder         "i-disorder"
#define CLOCK_NAME_MatInversion      "Matrix inversion"
#define CLOCK_NAME_MatMult           "Matrix multiplication"
#ifdef MAGMA
#define CLOCK_NAME_gpuMatInversion   "Matrix inversion GPU"
#define CLOCK_NAME_gpuMatMult        "Matrix multiplication GPU"
#endif
