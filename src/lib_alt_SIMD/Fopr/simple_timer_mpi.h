/*
  simple timer class with MPI_Wtime
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2019 Jan. 3

*/
#ifndef SIMPLE_TIMER_MPI_INC_H
#define SIMPLE_TIMER_MPI_INC_H
#include <mpi.h>
#include "IO/bridgeIO.h"

class simple_timer_mpi{

 public:
 simple_timer_mpi(const char* str): diff_sec(0), name(str) {}

  void reset(){
    diff_sec=0.0;
  }

  void start(){
    start_time=MPI_Wtime();
  }
  void stop(){
    end_time=MPI_Wtime();
    diff_sec+=(end_time-start_time);
  }

  double elapsed_msec() const {
    double diff_total_msec=diff_sec*1.0e+3;
    return diff_total_msec;
  }

  double elapsed_sec() const {
    return diff_sec;
  }

  void print_resolution(const Bridge::VerboseLevel vl = Bridge::GENERAL){
    double ticks=MPI_Wtick();
    double msec=1.0e3*ticks;
    Bridge::vout.general(vl, "time resolution: %g msec.  (MPI version)\n",msec);
  }

  void report(const Bridge::VerboseLevel vl = Bridge::GENERAL){
    Bridge::vout.general(vl,"Elapsed time: %s: total %12.4f sec\n", name.c_str(), elapsed_sec());
  }



 private:
  double start_time, end_time;
  double diff_sec;
  std::string name;
};
#endif
