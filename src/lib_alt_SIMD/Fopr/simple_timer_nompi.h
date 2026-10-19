/*
  simple timer class with clock_gettime in CLOCK_MONOTONIC mode
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2019 Jan. 3
*/
#ifndef SIMPLE_TIMER_NOMPI_INC_H
#define SIMPLE_TIMER_NOMPI_INC_H
#include <time.h>
#include "IO/bridgeIO.h"

class simple_timer_nompi{

 public:
 simple_timer_nompi(const char *str): diff_sec(0), diff_nsec(0), name(str){}

  void reset(){
    diff_sec=0.0;
    diff_nsec=0.0;
  }

  void start(){
    clock_gettime(CLOCK_MONOTONIC, &start_time);
  }
  void stop(){
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    diff_sec+=(end_time.tv_sec-start_time.tv_sec);
    diff_nsec+=(end_time.tv_nsec-start_time.tv_nsec);
  }

  double elapsed_msec() const {
    double diff_total_msec=diff_sec*1.0e+3+diff_nsec*1.0e-6;
    return diff_total_msec;
  }

  double elapsed_sec() const {
    double diff_total_sec=diff_sec+diff_nsec*1.0e-9;
    return diff_total_sec;
  }


  static void print_resolution(const Bridge::VerboseLevel vl = Bridge::GENERAL){
    timespec res;
    clock_getres(CLOCK_MONOTONIC, &res);
    double msec=res.tv_sec*1.0e+3+res.tv_nsec*1.0e-6;
    vout.general(vl, "time resolution: %g msec.  (clock_gettime version)\n",msec);
  }

  void report(const Bridge::VerboseLevel vl = Bridge::GENERAL){
    Bridge::vout.general(vl,"Elapsed time: %s: total %12.4f sec\n", name.c_str(), elapsed_sec());
  }

 private:
  struct timespec start_time, end_time;
  double diff_sec, diff_nsec;
  std::string name;
};
#endif
