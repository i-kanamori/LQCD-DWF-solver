/*!

        @file    afpr_Domainwall_5din_float.cpp

        @brief   Domainwall operator (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate:  $

        @version $LastChangedRevision:  $

 */
//====================================================================
#include "afopr_Domainwall_5din.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_SIMD/Field/afield.h"
#include "site_list_utils.h"
#include "debug_wilson.h"

// inline functions
#include "lib_alt_SIMD/inline/vector_arith-inc.h"
#include "lib_alt_SIMD/inline/vector_shift_float-inc.h"   // specific to single prec.

// function template files
#include "lib_alt_SIMD/Field/afield-inc.h"
#include "lib_alt_SIMD/Field/afield_Gauge-inc.h"

// template functions
#include "afopr_Domainwall_5din-tmpl.h"


#if defined USE_GROUP_SU3
typedef AField<float,SIMD> AFIELDspec;
#include "afopr_Domainwall_5din_SU3-inc.h"
//#include "afopr_Domainwall_5din_block_SU3-inc.h"
#else
#error("valid only for SU3")
#endif


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float,SIMD> >::Factory_params::Register(
                              "Domainwall_5din", create_object_with_params);

}
#endif


// explicit instanciation.
template class AFopr_Domainwall_5din<AField<float,SIMD> >;
template<>
const std::string AFopr_Domainwall_5din<AField<float,SIMD> >::class_name = "AFopr_Domainwall_5din<AField<float,SIMD> >";

//============================================================END=====
