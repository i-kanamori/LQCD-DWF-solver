/*!
      @File    test_alt_QXS.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-10-25 13:27:53 #$
      @version $LastChangedRevision: 2317 $
*/

#ifdef USE_ALT_QXS

#include "lib_alt_QXS/bridge_alt_qxs.h"
#include "testlist_alt.h"
#include "spectrum_Domainwall_alt-tmpl.h"

#define IMPL QXS

namespace Test_alt_QXS{
  int test_Spectrum_Domainwall_5din();
  bool use_mpfr = false;
}

//====================================================================
int Test_alt_QXS::test_all()
{
  int result = 0;
  result += test_Spectrum_Domainwall_5din();

  vout.general("Test result for QXS = %d\n", result);
  return result;
}


/*
double calc_plaquette(const Field_G *u){
  // parepare parameter
  Parameters params_staple;
  params_staple.set_string("verbose_level", "general");
  Timer timer("calc_plqauette");
  timer.start();

  // calculate plaquette
  const unique_ptr<Staple> staple(Staple::New("Lexical"));
  staple->set_parameters(params_staple);
  const double result = staple->plaquette(*u);

  // done!
  timer.stop();
  timer.report();
  return result;
}

void apply_smearing(Field_G *u, const Parameters& params){
  const string str_smear_type = params.get_string("smear_type");
  if(str_smear_type == "None"){
    vout.general("skip smearing\n");
    return;
  }
  const string str_proj_type  = params.get_string("projection_type");
  //    const double rho_uniform    = params.get_double("rho_uniform");
  //    const int number_of_smearing= params.get_int("number_of_smearing");
  //    const string str_vlevel     = params.get_string("verbose_level");

  // parameters
  const Parameters &params_proj=params;
  const Parameters &params_smear=params;
  const Parameters &params_dr_smear=params;

  vout.general("applying smearing\n");
  Timer timer_smearing("smearig");
  timer_smearing.reset();
  timer_smearing.start();

  // projection
  unique_ptr<Projection> proj(Projection::New(str_proj_type));
  proj->set_parameters(params_proj);

  // smear
  unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
  smear->set_parameters(params_smear);

  // director
  const unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear));
  dr_smear->set_parameters(params_dr_smear);
  dr_smear->set_config(u);

  // run smearing
  const int     Nsmear  = dr_smear->get_Nsmear();
  const Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);
  copy(*u, *Usmear);

  // done!
  timer_smearing.stop();
  timer_smearing.report();
}
*/


//====================================================================
int Test_alt_QXS::test_Spectrum_Domainwall_5din()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  int result = 0;
  if(Test_alt_QXS::use_mpfr){
    result += test_dw.hadron_2ptFunction(test_file, "richardson_eo_mpreal");
    return result;
  }

  //  result += test_dw.hadron_2ptFunction(test_file, "double");
  //  result += test_dw.hadron_2ptFunction(test_file, "float");
  //  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  //  result += test_dw.hadron_2ptFunction(test_file, "float_eo");
  //  result += test_dw.hadron_2ptFunction(test_file, "mixed");
  //  result += test_dw.hadron_2ptFunction(test_file, "mixed_eo");
  result += test_dw.hadron_2ptFunction(test_file, "richardson_eo");
  // unchecked:
  //  result += test_dw.hadron_2ptFunction(test_file, "org");
  //result += test_dw.hadron_2ptFunction("double_prec");
  //test_dw.hadron_2ptFunction("float");
  //test_dw.hadron_2ptFunction("mixed");
  return result;
}


#endif
//============================================================END=====
