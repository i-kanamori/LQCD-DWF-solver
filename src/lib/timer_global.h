#ifndef TIMER_GLOBAL_H
#define TIMER_GLOBAL_H

#ifdef USE_TIMER_GLOBAL

#define TIMER_MPI_BARRIER_BEFORE_ALLREDUCE

#include "lib/Tools/timer.h"


void timer_global_init();
void timer_global_report();


extern unique_ptr<Timer> timer_global_comm;
extern unique_ptr<Timer> timer_global_comm_send_wait;
extern unique_ptr<Timer> timer_global_linear_algebra;
extern unique_ptr<Timer> timer_global_mult1;
extern unique_ptr<Timer> timer_global_multb;
extern unique_ptr<Timer> timer_global_mult2;
extern unique_ptr<Timer> timer_global_mult_diag;
extern unique_ptr<Timer> timer_global_copy;
extern unique_ptr<Timer> timer_global_reduction;
extern unique_ptr<Timer> timer_global_mpi_all_reduce;

extern unique_ptr<Timer> timer_global_comm_float;
extern unique_ptr<Timer> timer_global_comm_send_wait_float;
extern unique_ptr<Timer> timer_global_linear_algebra_float;
extern unique_ptr<Timer> timer_global_mult1_float;
extern unique_ptr<Timer> timer_global_multb_float;
extern unique_ptr<Timer> timer_global_mult2_float;
extern unique_ptr<Timer> timer_global_mult_diag_float;
extern unique_ptr<Timer> timer_global_copy_float;
extern unique_ptr<Timer> timer_global_reduction_float;
extern unique_ptr<Timer> timer_global_mpi_all_reduce_float;

#ifdef USE_TIMER_FLOAT
#define TIMER_GLOBAL_START(name, cond)                                       \
  if((cond)) timer_global_ ## name ## _float->start();

#define TIMER_GLOBAL_STOP(name, cond)                                        \
  if((cond)) timer_global_ ## name ## _float ->stop();
#else
#define TIMER_GLOBAL_START(name, cond)                                       \
  if((cond)) timer_global_ ## name->start();

#define TIMER_GLOBAL_STOP(name, cond)                                        \
  if((cond)) timer_global_ ## name->stop();

#endif

// with 3 arguments
#define TIMER_GLOBAL_START3(name, cond, type)                             \
  if((sizeof(type)==4) && cond) timer_global_ ## name ## _float->start(); \
  if((sizeof(type)==8) && cond) timer_global_ ## name ->start();

#define TIMER_GLOBAL_STOP3(name, cond, type)                                \
  if((sizeof(type)==4) && (cond)) timer_global_ ## name ## _float ->stop(); \
  if((sizeof(type)==8) && (cond)) timer_global_ ## name ->stop();

#else

void timer_global_init() {};
void timer_global_report() {};

#define TIMER_GLOBAL_START(name, cond)
#define TIMER_GLOBAL_STOP(name, cond)
#define TIMER_GLOBAL_START3(name, cond, size)
#define TIMER_GLOBAL_STOP3(name, cond, size)

#endif


#endif
