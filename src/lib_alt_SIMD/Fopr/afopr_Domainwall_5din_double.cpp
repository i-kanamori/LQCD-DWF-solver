/*!

        @file    afpr_Domainwall_5din_double.cpp

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
#include "lib_alt_SIMD/inline/vector_shift_double-inc.h"   // specific to double prec.

// function template files
#include "lib_alt_SIMD/Field/afield-inc.h"
#include "lib_alt_SIMD/Field/afield_Gauge-inc.h"

// template functions
#include "afopr_Domainwall_5din-tmpl.h"


#if defined USE_GROUP_SU3
typedef AField<double,SIMD> AFIELDspec;
#include "afopr_Domainwall_5din_SU3-inc.h"
//#include "afopr_Domainwall_5din_block_SU3-inc.h"
#else
#error("valid only for SU3")
#endif


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<double,SIMD> >::Factory_params::Register(
                              "Domainwall_5din", create_object_with_params);

}
#endif


// explicit instanciation.
template class AFopr_Domainwall_5din<AField<double,SIMD> >;
template<>
const std::string AFopr_Domainwall_5din<AField<double,SIMD> >::class_name = "AFopr_Domainwall_5din<AField<double,SIMD> >";

//============================================================END=====
