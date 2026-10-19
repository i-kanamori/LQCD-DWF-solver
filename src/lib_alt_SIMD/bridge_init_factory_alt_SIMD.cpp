/*!
        @file    bridge_init_factory_alt_SIMD.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2019-07-30 20:13:19 #$
        @version $LastChangedRevision: 1983 $
*/

#include "bridge_init_factory_alt_SIMD.h"

#ifdef USE_FACTORY

// alt-code
#include "lib/Fopr/afopr.h"
//#include "lib/Force/Fermion/aforce_F.h"

#include "lib_alt/Solver/asolver.h"

#include "lib_alt_SIMD/Field/afield.h"

//#include "lib/Eigen/aeigensolver.h"
// adding to corelib
//#include "lib/Action/action.h"

/*
#include "Force/Gauge/force_G.h"
#include "Smear/forceSmear.h"
#include "Smear/smear.h"
#include "Smear/projection.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Measurements/Gauge/staple.h"
#include "Measurements/Fermion/source.h"
#include "Tools/randomNumbers.h"
#include "Tools/gammaMatrixSet.h"
*/

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_SIMD()
{
  bool result = true;
  typedef AField<double,SIMD> AFIELD_d;
  typedef AField<float,SIMD> AFIELD_f;

  result &= AFopr<AFIELD_d>::init_factory();
  result &= ASolver<AFIELD_d>::init_factory();
  //  result &= AEigensolver<AFIELD_d, AFopr<AFIELD_d> >::init_factory();
  // result &= AForce_F<AField<double,SIMD> >::init_factory();

  result &= AFopr<AFIELD_f>::init_factory();
  result &= ASolver<AFIELD_f>::init_factory();


  /*
  result &= Force_G::init_factory();
  result &= ForceSmear::init_factory();
  result &= Solver::init_factory();
  result &= Action::init_factory();
  result &= Projection::init_factory();
  result &= Smear::init_factory();
  result &= GaugeFixing::init_factory();
  result &= Source::init_factory();
  result &= Staple::init_factory();
  result &= RandomNumbers::init_factory();
  result &= GammaMatrixSet::init_factory();
#ifdef USE_FFTWLIB
  result &= FFT::init_factory();
#endif
  */

  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#ifdef DEBUG
void bridge_report_factory_alt_SIMD()
{
  vout.general("------------------------------------------------\n");
  vout.general("Factory entries\n");
  vout.general("------------------------------------------------\n");
  /*
  Fopr::Factory_noarg::print("Fopr(void)");
  Fopr::Factory_fopr::print("Fopr(Fopr*)");
  Fopr::Factory_fopr_director::print("Fopr(Fopr*, Director*)");
  Fopr::Factory_string::print("Fopr(string)");

  Force_G::Factory::print("Force_G");

  ForceSmear::Factory::print("ForceSmear");

  Solver::Factory::print("Solver");

  Action::Factory::print("Action");

  Projection::Factory::print("Projection");
  Smear::Factory::print("Smear");

  GaugeFixing::Factory::print("GaugeFixing");

  Source::Factory::print("Source");

  Staple::Factory::print("Staple");

  RandomNumbers::Factory_int::print("RandomNumbers(int)");
  RandomNumbers::Factory_file::print("RandomNumbers(string)");

  GammaMatrixSet::Factory::print("GammaMatrixSet");
  */
  vout.general("------------------------------------------------\n");

  return;
}

#endif

#endif /* USE_FACTORY */
