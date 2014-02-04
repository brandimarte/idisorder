#include "timer_defs.h"

// For fortran:
#define clock_init         clock_init_
#define clock_start        clock_start_
#define clock_stop         clock_stop_
#define clock_print_last   clock_print_last_
#define clock_print_all    clock_print_all_

extern "C" void clock_init(void);
extern "C" void clock_start(const int *Routine_id);
extern "C" void clock_stop(const int *Routine_id);
extern "C" void clock_print_last(const int *Routine_id);
extern "C" void clock_print_all(void);
