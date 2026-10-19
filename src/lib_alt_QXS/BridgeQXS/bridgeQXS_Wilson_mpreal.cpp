/*!
      @file    bridgeQXS_Wilson_mpreal.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2492 $
*/

#ifdef USE_MPREAL

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"
#include "bridge_mpreal.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

using mpfr::mpreal;
typedef mpreal real_t;

#include "lib_alt_QXS/inline/vsimd_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_mpreal-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_mpreal-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"

#include "src/mult_common_parts_qxs-inc.h"
#include "src/mult_Wilson_parts_qxs-inc.h"
#include "src/mult_Wilson_parts_qxs2-inc.h"
#include "src/mult_Wilson_eo_parts_qxs-inc.h"

#include "src/mult_Wilson_qxs-inc.h"
#include "src/mult_Wilson_eo_qxs-inc.h"


#endif
