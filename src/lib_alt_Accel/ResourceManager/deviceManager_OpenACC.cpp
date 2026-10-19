/*!
        @file    deviceManager_OpenACC.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/
#ifdef USE_ACCEL_OPENACC

#include "ResourceManager/deviceManager_OpenACC.h"

// the following code only applies to NVIDIA
#ifdef _ACCEL
#include "accel.h"
#endif

#include "lib/Communicator/communicator.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

int DeviceManager_OpenACC::m_num_devices = 1;
int DeviceManager_OpenACC::m_device_id  = 0;
bool DeviceManager_OpenACC::m_initialized = false;

Bridge::VerboseLevel DeviceManager_OpenACC::m_vl = Bridge::CRUCIAL;

const std::string DeviceManager_OpenACC::class_name
                                           = "DeviceManager_OpenACC";

namespace {

int get_local_rank_from_env(void) {
    const char *candidates[] = {
        "SLURM_LOCALID",           // SLURM (srun)
        "OMPI_COMM_WORLD_LOCAL_RANK", // Open MPI (mpirun)
        "MV2_COMM_WORLD_LOCAL_RANK",  // MVAPICH2
        "MPI_LOCALRANKID",         // Intel MPI / MS-MPI
        "PALS_LOCAL_RANKID",       // Cray PALS
        "PMI_LOCAL_RANK",          // ggeneral PMI
        "LOCAL_RANK",              // PyTorch etc.
        NULL
    };

    int val_max=-1;
    for (int i = 0; candidates[i] != NULL; i++) {
        const char *val = getenv(candidates[i]);
        //        printf("get_local_rank_from_env: i=%d, candidate=%s, val=%s\n", i, candidates[i], val);
        if (val != NULL) {
          int val_this=atoi(val);
          // some variables may be available but might be set to 0 in all processes
          if (val_this > val_max) {
            val_max=val_this;
          }
        }
    }
    return val_max; // -1 if not found
}
}

//====================================================================
void DeviceManager_OpenACC::init(int device_num_offset = 0)
{
  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: being initialized.\n", class_name.c_str());

  if ( DeviceManager_OpenACC::m_initialized ) {
    vout.general(m_vl, "%s: already done, skipping\n", class_name.c_str());
    return;
  }

  int process_id = Communicator::nodeid();
  vout.general(m_vl, "  process_id = %d\n", process_id);

#ifdef _ACCEL // the following code only applies to NVIDIA
  vout.general(m_vl, "  NVIDIA environment is used.\n");
  m_num_devices = acc_get_num_devices(acc_device_nvidia);
  m_device_id = process_id % m_num_devices + device_num_offset;
  acc_set_device_num(m_device_id, acc_device_nvidia);
#endif
  vout.general(m_vl, "  number of devices = %d\n", m_num_devices);
  vout.general(m_vl, "  device_id = %d  is used on process = %d\n",
               m_device_id, process_id);

  DeviceManager_OpenACC::m_initialized=true;
}

//====================================================================
void DeviceManager_OpenACC::init_prempi(int device_num_offset = 0)
{
  m_vl = CommonParameters::Vlevel();

  int process_id = get_local_rank_from_env();
  if ( process_id < 0 ) {
    printf("No process id is found in the enviromental variables, which is needed to initialize MPI with device aware MPI.\n");
    abort();
  }
  //  printf("init_rempi: process_id=%d\n", process_id);
  if ( process_id == 0 ) {
    vout.general(m_vl, "%s: being initialized before MPI_Init\n", class_name.c_str());
    vout.general(m_vl, "  process_id = %d\n", process_id);
  }
#ifdef _ACCEL // the following code only applies to NVIDIA
  if ( process_id == 0 ) {
    vout.general(m_vl, "  NVIDIA environment is used.\n");
  }
  m_num_devices = acc_get_num_devices(acc_device_nvidia);
  m_device_id = process_id % m_num_devices + device_num_offset;
  acc_set_device_num(m_device_id, acc_device_nvidia);
#endif
  if ( process_id == 0 ) {
    vout.general(m_vl, "  number of devices = %d\n", m_num_devices);
    vout.general(m_vl, "  device_id = %d  is used on local rank = %d\n",
                 m_device_id, process_id);
  } else {
    vout.paranoiac(m_vl, "  number of devices = %d\n", m_num_devices);
    vout.paranoiac(m_vl, "  device_id = %d  is used on local rank = %d\n",
                   m_device_id, process_id);
  }
  DeviceManager_OpenACC::m_initialized=true;
}

//====================================================================
void DeviceManager_OpenACC::finalize()
{
  vout.paranoiac(m_vl, "%s: finalization.\n", class_name.c_str());
}

#endif
//============================================================END=====
