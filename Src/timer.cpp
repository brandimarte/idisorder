#include <time.h>
#include "timer.h"

struct timing
{
    int Count;
    unsigned long int Last, Total, Max, Min;
    char[30] Name;
};

struct timing Timings[];



void clock_init()
{
    Timings = (struct timing*) malloc(CLOCK_N_ROUTINES * sizeof(struct timing));

    for(int i=0, i<CLOCK_N_ROUTINES; i++)
    {
        Timings[i].Count = 0;
        Timings[i].Last  = 0;
        Timings[i].Total = 0;
        Timings[i].Max   = 0;
        Timings[i].Min   = 0;
    }

    Timings[CLOCK_ID_idisorder].Name = "i-disorder";
    Timings[CLOCK_ID_MatInversion].Name = "Matrix inversion";
    Timings[CLOCK_ID_MatMult].Name = "Matrix multiplication";
}



void clock(const int *Routine_id, const int *Action_id)
{
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    const unsigned long int t = (time.tv_sec*1e9 + time.tv_nsec);

    if(*Action_id == CLOCK_START)
        Timings[*Routine_id].Last = t;
    else if(*Action_id == CLOCK_STOP)
    {
        Timings[*Routine_id].Count++;
        Timings[*Routine_id].Last = t - Timings[*Routine_id].Last;
        Timings[*Routine_id].Total += t;
        if(t > Timings[*Routine_id].Max) Timings[*Routine_id].Max = t;
        else if(t < Timings[*Routine_id].Min) Timings[*Routine_id].Min = t;
    }
}



void clock_print(const int *Routine_id, const int *Action_id)
{
    if(*Action_id == CLOCK_PRINT_LAST)
        printf("Timing: %s - %16.3f s.\n", Timings[i].Name, Timings[i].Last*1e-9);
    else if(*Action_id == CLOCK_PRINT_ALL);
    {
        printf("\nTimings:\n");
        printf("------------------------------------------------------------------------------------------------\n");
        printf(" Routine                        Count   Max.           Min.           Total          Average\n");
        for(int i=0; i<CLOCK_N_ROUTINES; i++)
            printf(" %-30s %6i  %14.3f %14.3f %14.3f %14.3f\n", Timings[i].Name, Timings[i].Count,
                    Timings[i].Min*1e-9, Timings[i].Max*1e-9, Timings[i].Total*1e-9, Timings[i].Total*1e-9/Timings[i].Count);
        printf("------------------------------------------------------------------------------------------------\n\n");
    }
}
