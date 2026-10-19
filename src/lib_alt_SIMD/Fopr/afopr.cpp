/*!
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2019-12-01 11:57:54 #$
        @version $LastChangedRevision: 2031 $
*/

#include "lib/Fopr/afopr.h"
#include "lib_alt_SIMD/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

//#include "lib_alt_SIMD/Fopr/afopr_Wilson.h"
//#include "lib_alt_SIMD/Fopr/afopr_Wilson_eo.h"
//#include "lib_alt_SIMD/Fopr/afopr_Clover.h"
//#include "lib_alt_SIMD/Fopr/afopr_Clover_eo.h"
//#include "lib_alt_SIMD/Fopr/afopr_Clover_dd.h"
//#include "lib_alt_SIMD/Fopr/afopr_Staggered.h"
//#include "lib_alt_SIMD/Fopr/afopr_Staggered_eo.h"
#include "lib_alt_SIMD/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_SIMD/Fopr/afopr_Domainwall_5din_eo.h"
//#include "lib/Fopr/afopr_Domainwall.h"
//#include "lib/Fopr/afopr_Domainwall_eo.h"
//#include "lib/Fopr/afopr_Smeared.h"
//#include "lib/Fopr/afopr_Rational.h"


template<>
bool AFopr<AField<double,SIMD> >::init_factory()
{
  typedef AField<double,SIMD> AFIELD;
  bool result = true;
  //  result &= AFopr_Wilson<AFIELD>::register_factory();
  //  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  //  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  //  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  return result;
}

template<>
bool AFopr<AField<float,SIMD> >::init_factory()
{
  typedef AField<float,SIMD> AFIELD;
  bool result = true;
  //result &= AFopr_Wilson<AFIELD>::register_factory();
  //  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  //  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  //  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


template<>
const std::string AFopr<AField<double,SIMD> >::
    class_name = "AFopr<AField<double,SIMD> >";

template<>
const std::string AFopr<AField<float,SIMD> >::
    class_name = "AFopr<AField<float,SIMD> >";


// explicit instanciation.
template class AFopr<AField<float,SIMD> >;
template class AFopr<AField<double,SIMD> >;

//============================================================END=====
