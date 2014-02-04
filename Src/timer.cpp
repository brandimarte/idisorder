#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "timer.h"

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
    const double t = (double)time.tv_sec + (double)time.tv_nsec*1e-9;
    const int i = *Routine_id;

    Timings[i].Count++;
    Timings[i].Last = t - Timings[i].Last;
    Timings[i].Total += t;
    if(t > Timings[i].Max) Timings[i].Max = t;
    else if(t < Timings[i].Min) Timings[i].Min = t;
}



void clock_print_last(const int *Routine_id)
{
    printf("Timing: %s - %16.3f s.\n", Timings[*Routine_id].Name, Timings[*Routine_id].Last);
}



void clock_print_all(void)
{
    printf("\nTimings (in seconds):\n");
    printf("---------------------------------------------------------------------------------------------------\n");
    printf(" Routine                          Count             Max.           Min.          Total      Average\n");
    for(int i=0; i<CLOCK_N_ROUTINES; i++)
        printf(" %-30s %6i  %14.3f %14.3f %14.3f %14.3f\n", Timings[i].Name, Timings[i].Count,
                Timings[i].Max, Timings[i].Min, Timings[i].Total, Timings[i].Total/(double)Timings[i].Count);
    printf("---------------------------------------------------------------------------------------------------\n\n");
}



#ifdef MPI
#include <mpi.h>
void clock_gather(int *node, int *root, MPI_Comm *MyComm)
{
    for(int i=0; i<CLOCK_N_ROUTINES; i++)
    {
        int count=0;
        MPI_Reduce(&Timings[i].Count, &count, 1, MPI_INT, MPI_SUM, *root, MyComm);
        if(*node == *root) Timings[i].Count = count;

        double time=0.0;
        MPI_Reduce(&Timings[i].Total, &time, 1, MPI_DOUBLE, MPI_SUM, *root, MyComm);
        if(*node == *root) Timings[i].Total = time;

        time=0.0;
        MPI_Reduce(&Timings[i].Max, &time, 1, MPI_DOUBLE, MPI_MAX, *root, MyComm);
        if(*node == *root) Timings[i].Max = time;

        time=0.0;
        MPI_Reduce(&Timings[i].Min, &time, 1, MPI_DOUBLE, MPI_MIN, *root, MyComm);
        if(*node == *root) Timings[i].Min = time;
    }
}
#endif
