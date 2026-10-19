
/*!
      @file    bridgeACC_char.cu
      @brief
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate:  #$
      @version $LastChangedRevision:  $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

#include "bridgeACC.h"

typedef char real_t;  // for channel with gpu_aware_mpi

namespace BridgeACC {

#include "src/bridge_cuda-inc.h"

}

//============================================================END=====
