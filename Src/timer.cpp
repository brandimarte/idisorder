#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "timer.h"
#ifdef MACH
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#endif


struct timing
{
    int Count;
    double Last, Total, Max, Min;  // in seconds
    char *Name;
};

struct timing *Timings;



// For fortran:
#define clock_init         clock_init_
#define clock_start        clock_start_
#define clock_stop         clock_stop_
#define clock_print_last   clock_print_last_
#define clock_print_all    clock_print_all_



void clock_init()
{
    Timings = (struct timing*) malloc(CLOCK_N_ROUTINES * sizeof(struct timing));

    for(int i=0; i<CLOCK_N_ROUTINES; i++)
    {
        Timings[i].Count = 0.0;
        Timings[i].Last  = 0.0;
        Timings[i].Total = 0.0;
        Timings[i].Max   = 0.0;
        Timings[i].Min   = 1.0e10;
    }
    Timings[CLOCK_ID_idisorder].Name       = CLOCK_NAME_idisorder;
    Timings[CLOCK_ID_MatInversion].Name    = CLOCK_NAME_MatInversion;
    Timings[CLOCK_ID_MatMult].Name         = CLOCK_NAME_MatMult;
#ifdef MAGMA
    Timings[CLOCK_ID_gpuMatInversion].Name = CLOCK_NAME_gpuMatInversion;
    Timings[CLOCK_ID_gpuMatMult].Name      = CLOCK_NAME_gpuMatMult;
#endif
}



void clock_start(const int *Routine_id)
{
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);

    Timings[*Routine_id].Last = (double)time.tv_sec + (double)time.tv_nsec*1e-9;
}



void clock_stop(const int *Routine_id)
{
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    const int i = *Routine_id;
    const double t = (double)time.tv_sec + (double)time.tv_nsec*1e-9 - Timings[i].Last;

    Timings[i].Count++;
    Timings[i].Last = t;
    Timings[i].Total += t;
    if(t > Timings[i].Max) Timings[i].Max = t;
    if(t < Timings[i].Min) Timings[i].Min = t;
}



void clock_print_last(const int *Routine_id)
{
    printf("Timing: %s - %16.6f s.\n", Timings[*Routine_id].Name, Timings[*Routine_id].Last);
}



void clock_print_all(void)
{
    printf("\nTimings (in seconds):\n");
    printf("----------------------------------------------------------------------------------------------------\n");
    printf(" Routine                          Count            Max.           Min.          Total        Average\n");

    for(int i=0; i<CLOCK_N_ROUTINES; i++)
    {
        if(Timings[i].Count != 0)
            printf(" %-30s %8i %14.3g %14.3g %14.3g %14.3g\n", Timings[i].Name, Timings[i].Count,
                    Timings[i].Max, Timings[i].Min, Timings[i].Total, Timings[i].Total/Timings[i].Count);
        else
            printf(" %-30s %8i           ---            ---            ---           ---\n", Timings[i].Name, Timings[i].Count);
    }

    printf("----------------------------------------------------------------------------------------------------\n\n");
    fflush(stdout);
}



#ifdef MPI
void clock_gather(int *node, int *root, int *MyComm)
{
    // i=0 is already in the root node(=0)
    for(int i=1; i<CLOCK_N_ROUTINES; i++)
    {
        int count=0;
        MPI_Reduce(&(Timings[i].Count), &count, 1, MPI_INT, MPI_SUM, *root, MPI_Comm_f2c(*MyComm));
        if(*node == *root) Timings[i].Count = count;

        double time=0.0;
        MPI_Reduce(&(Timings[i].Total), &time, 1, MPI_DOUBLE, MPI_SUM, *root, MPI_Comm_f2c(*MyComm));
        if(*node == *root) Timings[i].Total = time;

        time=0.0;
        MPI_Reduce(&(Timings[i].Max), &time, 1, MPI_DOUBLE, MPI_MAX, *root, MPI_Comm_f2c(*MyComm));
        if(*node == *root) Timings[i].Max = time;

        time=0.0;
        MPI_Reduce(&(Timings[i].Min), &time, 1, MPI_DOUBLE, MPI_MIN, *root, MPI_Comm_f2c(*MyComm));
        if(*node == *root) Timings[i].Min = time;
    }
}
#endif
