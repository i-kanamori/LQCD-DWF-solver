/*!
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib/Fopr/afopr.h"
#include "lib_alt_QXS/Field/afield.h"
#include "bridge_mpreal.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_QXS/Fopr/afopr_Wilson.h"
#include "lib_alt_QXS/Fopr/afopr_Wilson_eo.h"
//#include "lib_alt_QXS/Fopr/afopr_Clover.h"
//#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
//#include "lib_alt_QXS/Fopr/afopr_Clover_eo.h"
//#include "lib_alt_QXS/Fopr/afopr_Staggered.h"
//#include "lib_alt_QXS/Fopr/afopr_Staggered_eo.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h"
#include "lib/Fopr/afopr_Domainwall.h"
#include "lib/Fopr/afopr_Domainwall_eo.h"
//#include "lib/Fopr/afopr_Smeared.h"
//#include "lib/Fopr/afopr_Rational.h"


template<>
bool AFopr<AField<double, QXS> >::init_factory()
{
  typedef AField<double, QXS> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  //  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
  //  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  return result;
}


template<>
bool AFopr<AField<float, QXS> >::init_factory()
{
  typedef AField<float, QXS> AFIELD;
  bool result = true;
  result &= AFopr_Wilson<AFIELD>::register_factory();
  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  //  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  result &= AFopr_Domainwall<AFIELD>::register_factory();
  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
  //  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  return result;
}

#ifdef USE_MPREAL
using mpfr::mpreal;
template<>
bool AFopr<AField<mpreal, QXS> >::init_factory()
{
  typedef AField<mpreal, QXS> AFIELD;
  bool result = true;
  //  result &= AFopr_Wilson<AFIELD>::register_factory();
  //  result &= AFopr_Wilson_eo<AFIELD>::register_factory();
  //  result &= AFopr_Clover<AFIELD>::register_factory();
  //  result &= AFopr_Clover_dd<AFIELD>::register_factory();
  //  result &= AFopr_Clover_eo<AFIELD>::register_factory();
  //  result &= AFopr_Staggered<AFIELD>::register_factory();
  //  result &= AFopr_Staggered_eo<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_5din<AFIELD>::register_factory();
  result &= AFopr_Domainwall_5din_eo<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall<AFIELD>::register_factory();
  //  result &= AFopr_Domainwall_eo<AFIELD>::register_factory();
  //  result &= AFopr_Smeared<AFIELD>::register_factory();
  //  result &= AFopr_Rational<AFIELD>::register_factory();
  return result;
}
#endif

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AFopr<AField<float, QXS> >;
template class AFopr<AField<double, QXS> >;

//============================================================END=====
