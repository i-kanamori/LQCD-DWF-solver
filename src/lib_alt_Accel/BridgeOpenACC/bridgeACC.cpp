/*!
      @file    bridsgeACC.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#include <cstddef>

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

#ifdef _OPENACC
#include "openacc.h"
#endif

namespace BridgeACC {

// real_t = double
//====================================================================
void enter_data_create(double* data, const size_t size, bool host_has_data)
{
#pragma acc enter data create(data[0:size])
}

//====================================================================
void exit_data_delete(double* data, const size_t size)
{
#pragma acc exit data delete(data[0:size])
}

//====================================================================
void update_device(double* data, size_t offset, size_t size)
{
#pragma acc update device(data[offset:size])
}

//====================================================================
void update_host(double* data, size_t offset, size_t size)
{
#pragma acc update host(data[offset:size])
}

// real_t = float
//====================================================================
void enter_data_create(float* data, const size_t size, bool host_has_data)
{
#pragma acc enter data create(data[0:size])
}

//====================================================================
void exit_data_delete(float* data, const size_t size)
{
#pragma acc exit data delete(data[0:size])
}

//====================================================================
void update_device(float* data, size_t offset, size_t size)
{
#pragma acc update device(data[offset:size])
}

//====================================================================
void update_host(float* data, size_t offset, size_t size)
{
#pragma acc update host(data[offset:size])
}

// char for communicator
//====================================================================
void enter_data_create(char* data, const size_t size, bool host_has_data)
{
#pragma acc enter data create(data[0:size])
  /*
#pragma acc data present(data[0:size]) copyin(size)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int i = 0; i < size; ++i){
      data[i] = 0.0;
    }
  }
  */
}

//====================================================================
void exit_data_delete(char* data, const size_t size)
{
#pragma acc exit data delete(data[0:size])
}



// void
//====================================================================
void* dev_ptr(void* phost)
{
#ifdef _OPENACC
  return acc_deviceptr(phost);
#else
  return phost;
#endif
}

} // namespace BridgeACC
//============================================================END=====

