#ifdef USE_TIMER_GLOBAL

#include "timer_global.h"



unique_ptr<Timer> timer_global_comm;
unique_ptr<Timer> timer_global_comm_send_wait;
unique_ptr<Timer> timer_global_linear_algebra;
unique_ptr<Timer> timer_global_mult1;
unique_ptr<Timer> timer_global_multb;
unique_ptr<Timer> timer_global_mult2;
unique_ptr<Timer> timer_global_mult_diag;
unique_ptr<Timer> timer_global_copy;
unique_ptr<Timer> timer_global_reduction;
unique_ptr<Timer> timer_global_mpi_all_reduce;

unique_ptr<Timer> timer_global_comm_float;
unique_ptr<Timer> timer_global_comm_send_wait_float;
unique_ptr<Timer> timer_global_linear_algebra_float;
unique_ptr<Timer> timer_global_mult1_float;
unique_ptr<Timer> timer_global_multb_float;
unique_ptr<Timer> timer_global_mult2_float;
unique_ptr<Timer> timer_global_mult_diag_float;
unique_ptr<Timer> timer_global_copy_float;
unique_ptr<Timer> timer_global_reduction_float;
unique_ptr<Timer> timer_global_mpi_all_reduce_float;

void timer_global_init(){
  timer_global_comm.reset(new Timer("comm"));
  timer_global_comm_send_wait.reset(new Timer("comm_send_wait"));
  timer_global_linear_algebra.reset(new Timer("linear algebra"));
  timer_global_mult1.reset(new Timer("mult1"));
  timer_global_multb.reset(new Timer("multb"));
  timer_global_mult2.reset(new Timer("mult2"));
  timer_global_mult_diag.reset(new Timer("mult_diag"));
  timer_global_copy.reset(new Timer("copy"));
  timer_global_reduction.reset(new Timer("global reduction"));
  timer_global_mpi_all_reduce.reset(new Timer("global reduction (MPI_Allreduce)"));

  timer_global_comm_float.reset(new Timer("comm [float]"));
  timer_global_comm_send_wait_float.reset(new Timer("comm send wait[float]"));
  timer_global_linear_algebra_float.reset(new Timer("linear algebra [float]"));
  timer_global_mult1_float.reset(new Timer("mult1 [float]"));
  timer_global_multb_float.reset(new Timer("multb [float]"));
  timer_global_mult2_float.reset(new Timer("mult2 [float]"));
  timer_global_mult_diag_float.reset(new Timer("mult_diag [float]"));
  timer_global_copy_float.reset(new Timer("copy [float]"));
  timer_global_reduction_float.reset(new Timer("global reduction [float]"));
  timer_global_mpi_all_reduce_float.reset(new Timer("global reduction (MPI_Allreduce) [float]"));
}

void timer_global_report() {

  vout.general("\n");
  vout.general("----------- timer -------------\n");
  timer_global_comm->report();
  timer_global_comm_send_wait->report();
  timer_global_linear_algebra->report();
  timer_global_mult1->report();
  timer_global_multb->report();
  timer_global_mult2->report();
  timer_global_mult_diag->report();
  timer_global_copy->report();
  timer_global_reduction->report();
  timer_global_mpi_all_reduce->report();
  vout.general("   copy contains conversion from float to double\n");

  vout.general("\n");
  vout.general("----------- timer [float] -------------\n");
  timer_global_comm_float->report();
  timer_global_comm_send_wait_float->report();
  timer_global_linear_algebra_float->report();
  timer_global_mult1_float->report();
  timer_global_multb_float->report();
  timer_global_mult2_float->report();
  timer_global_mult_diag_float->report();
  timer_global_copy_float->report();
  timer_global_reduction_float->report();
  timer_global_mpi_all_reduce_float->report();
  vout.general("   copy [float] contains conversion from double to float\n");

  vout.general("\n");
  double mult_arithmetics_float
    = timer_global_mult1_float->elapsed_sec()
    + timer_global_mult2_float->elapsed_sec()
    + timer_global_mult_diag_float->elapsed_sec();

  double multb_float = timer_global_multb_float->elapsed_sec();
  double comm_float  = timer_global_comm_float->elapsed_sec();
  vout.general("------------------------\n");
  vout.general("[float] arithmetics in mult (mult1+mult2+mult_diag) : %6.3f sec\n", mult_arithmetics_float);
  vout.general("[float] overlap: arithmetics in mult (multb) :        %6.3f sec\n", multb_float);
  vout.general("[float] overlap: communcation in mult (comm) :        %6.3f sec\n", comm_float);
  vout.general("------------------------\n");
}

#endif
