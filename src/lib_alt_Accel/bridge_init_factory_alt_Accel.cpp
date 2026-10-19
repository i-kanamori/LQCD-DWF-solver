/*!
        @file    bridge_init_factory_alt_Accel.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#include "lib_alt_Accel/bridge_init_factory_alt_Accel.h"

#ifdef USE_FACTORY

// alt-code
#include "lib/Fopr/afopr.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt/Solver/asolver.h"


#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_Accel()
{
  bool result = true;

  result &= AFopr<AField<double, ACCEL> >::init_factory();
  result &= ASolver<AField<double, ACCEL> >::init_factory();

  result &= AFopr<AField<float, ACCEL> >::init_factory();
  result &= ASolver<AField<float, ACCEL> >::init_factory();

  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

void bridge_report_factory_alt_Accel()
{
  vout.general("------------------------------------------------\n");
  vout.general("Alt_Accel: Factory entries\n");
  vout.general("------------------------------------------------\n");

  typedef AFopr<AField<double, ACCEL> > AFOPR_d;
  typedef AFopr<AField<float, ACCEL> > AFOPR_f;

  AFOPR_d::Factory_noarg::print("AFopr<double, ACCEL> >(void)");
  AFOPR_d::Factory_params::print("AFopr<double, ACCEL>(Parameters&)");
  AFOPR_d::Factory_fopr::print("AFopr<double, ACCEL>(Fopr*)");
  AFOPR_d::Factory_fopr_director::print("AFopr<double, ACCEL>(Fopr*, Director*)");

  AFOPR_f::Factory_noarg::print("AFopr<float, ACCEL> >(void)");
  AFOPR_f::Factory_params::print("AFopr<float, ACCEL>(Parameters&)");
  AFOPR_f::Factory_fopr::print("AFopr<float, ACCEL>(Fopr*)");
  AFOPR_f::Factory_fopr_director::print("AFopr<float, ACCEL>(Fopr*, Director*)");


  vout.general("------------------------------------------------\n");

  return;
}

#endif /* USE_FACTORY */

//============================================================END=====
