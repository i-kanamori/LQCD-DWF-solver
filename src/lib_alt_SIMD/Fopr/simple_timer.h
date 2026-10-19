/*
  simple timer class
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2019 Jan. 3

  both w/ and w/o MPI have the same interface
*/


#ifndef SIMPLE_TIMER_INC_H
#define SIMPLE_TIMER_INC_H

#ifdef USE_MPI
#include "simple_timer_mpi.h"
typedef simple_timer_mpi SimpleTimer;
#else
#include "simple_timer_nompi.h"
typedef simple_timer_nompi SimpleTimer;
#endif

#endif // include gurad
