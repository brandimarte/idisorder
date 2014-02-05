#include "timer_defs.h"
#ifdef MPI
# undef MPI
# define MPI MPI
# include "mpi.h"
//# define MPI
#endif

// For fortran:
#define clock_init         clock_init_
#define clock_start        clock_start_
#define clock_stop         clock_stop_
#define clock_print_last   clock_print_last_
#define clock_print_all    clock_print_all_
#define clock_gather       clock_gather_

extern "C"
{
    void clock_init(void);
    void clock_start(const int *Routine_id);
    void clock_stop(const int *Routine_id);
    void clock_print_last(const int *Routine_id);
    void clock_print_all(void);
#ifdef MPI
//# undef MPI
    void clock_gather(int *node, int *root, int *MyComm);
//    MPI_Comm *parallel_mp_mpi_comm_myworld_;
# define MPI_Comm_MyWorld parallel_mp_mpi_comm_myworld_
//# define MPI
#endif
}
