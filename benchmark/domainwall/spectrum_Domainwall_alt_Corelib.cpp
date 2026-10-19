/*!
      @File    spectrum_Domainwall_alt_Corelib.cpp
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2022-01-30 02:14:34 #$
      @version $LastChangedRevision: 2344 $
*/

//#include "spectrum_Domainwall_alt.h"
#include "test.h"

#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr_Domainwall.h"
#include "lib/Fopr/fopr_Domainwall_eo.h"
#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
#include "lib/Measurements/Fermion/fprop_Standard_eo.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/source.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"

#include "spectrum_Domainwall_alt.h"

// class name
const std::string Spectrum_Domainwall_alt<CORELIB>::class_name =
                                    "Spectrum_Domainwall_alt<CORELIB>";

//====================================================================
void Spectrum_Domainwall_alt<CORELIB>::init()
{
  // do nothing.
}

//====================================================================
int Spectrum_Domainwall_alt<CORELIB>::hadron_2ptFunction(
                            std::string test_file, std::string mode)
{
  const std::string test_name = class_name + ".hadron_2ptFunction";

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", test_file.c_str());
  vout.general("test mode: %s\n", mode.c_str());

  // parameter setup.
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(test_file);

  Parameters params_test = params_all.lookup("Test_Spectrum");
  Parameters params_fopr = params_all.lookup("Fopr");
  Parameters params_source = params_all.lookup("Source");

  const string str_vlevel       = params_test.get_string("verbose_level");
  Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);
  m_vl = vl;

  const bool   do_check        = params_test.is_set("expected_result");
  const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

  const string str_fopr_type   = params_fopr.get_string("fermion_type");
  const string str_source_type = params_source.get_string("source_type");

  vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
  vout.general(vl, "  source_type  = %s\n", str_source_type.c_str());

  RandomNumberManager::initialize("Mseries", 1234567UL);

  vout.general(vl, "fermion_type is %s, but is ignored\n", str_fopr_type.c_str());

  // configuration setup.
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // for correlation function
  string gmset_type  = params_fopr.get_string("gamma_matrix_type");
  vout.general(vl, "  gmset_type   = %s\n", gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));



  // full operator
  const string kernel_type = params_fopr.get_string("kernel_type");
  const double M0          = params_fopr.get_double("domain_wall_height");
  const double coeff_c     = params_fopr.get_double("coefficient_c");
  const int    Ns          = params_fopr.get_int("extent_of_5th_dimension");

  Parameters params_kernel   = params_fopr;
  double     kappa           = 1.0 / (8.0 - 2.0 * M0);
  double     one_over_2kappa = 4.0 - M0;
  params_kernel.set_double("hopping_parameter", kappa);
  unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
  unique_ptr<Fopr_Domainwall> fopr_dw(new Fopr_Domainwall(foprw.get(), params_fopr));
  fopr_dw->set_config(U);

  // even-odd operator
  unique_ptr<Fopr_Domainwall_eo> fopr_dw_eo;

  // solver
  Parameters params_solver = params_all.lookup("Solver");
  string str_solver_type = params_solver.get_string("solver_type");
  unique_ptr<Solver> solver;  // itrative solver
  unique_ptr<Fprop> fprop;    // contains pre/post process

  if(mode == "double"){
    solver.reset(Solver::New(str_solver_type, fopr_dw.get(), params_solver));
    fprop.reset(new Fprop_Standard_lex(solver.get()));
    fprop->set_mode("D");
  } else if(mode == "double_eo"){
    fopr_dw_eo.reset(new Fopr_Domainwall_eo(params_fopr));
    solver.reset(Solver::New(str_solver_type, fopr_dw_eo.get(), params_solver));
    fprop.reset(new Fprop_Standard_eo(solver.get()));
    fprop->set_mode("D");
    fprop->set_config(U.get());
  }else{
    vout.crucial(vl, "%s: irrelevant mode =%s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  //  fprop->mult_performance("DdagD", 100);

  // source
  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);


  // ####  Execution main part  ####
  Timer timer(test_name);
  timer.start();

  std::vector<Field_F> sq(Nc * Nd);
#pragma omp parallel
  {
    for (int i_cd = 0; i_cd < Nc * Nd; ++i_cd) {
      sq[i_cd].set(0.0);
    }
  }

  Field_F b;
  Field_F vtmp_p, vtmp_m, vtmp;

  Field b5(b.nin(), b.nvol(), Ns);
  Field y5(b.nin(), b.nvol(), Ns);
  Field x5(b.nin(), b.nvol(), Ns);

  vout.general(vl, "\n");
  vout.general(vl, "Solving quark propagator:\n");
  vout.general(vl, "  color spin   Nconv      diff           diff2\n");

  Timer timer_prop("calculating porpagator");
  timer_prop.reset();
  timer_prop.start();

  for (int ispin = 0; ispin < Nd; ++ispin) {
    for (int icolor = 0; icolor < Nc; ++icolor) {
      int i_cd = icolor + Nc * ispin;

      source->set(b, i_cd);

      // set 5d source as
      //  - D_- ( P_+ 0 ....0 P_-)^T b4
#pragma omp parallel
      {
        // build 5dim vector
        b5.set(0.0);

        // s5=0
        fopr_dw->mult_chproj_4d(vtmp_p, b, +1);
        fopr_dw->mult_chproj_4d(vtmp_m, b, -1);

        foprw->mult(vtmp, vtmp_p);
        axpy(vtmp_p, -one_over_2kappa * coeff_c, vtmp);  // (-cD+1) (1+gm5)b/2

        foprw->mult(vtmp, vtmp_m);
        axpy(vtmp_m, -one_over_2kappa * coeff_c, vtmp);  // (-cD+1) (1-gm5)b/2

        copy(b5, 0, vtmp_p, 0);
        copy(b5, Ns - 1, vtmp_m, 0);
      } // omp parallel

      int    Nconv;
      double diff;
      fprop->invert_D(x5, b5, Nconv, diff);
#pragma omp parallel
      {
        fopr_dw->set_mode("D");
        fopr_dw->mult(y5, x5);
        axpy(y5, -1.0, b5);
        double diff2 = y5.norm2()/b5.norm2();

        vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                     icolor, ispin, Nconv, diff, diff2);

        // convert to 4D propagator
        copy(vtmp, 0, x5, 0);
        fopr_dw->mult_chproj_4d(vtmp_m, vtmp, -1);
        copy(vtmp, 0, x5, Ns - 1);
        fopr_dw->mult_chproj_4d(vtmp_p, vtmp, +1);

        copy(sq[i_cd], vtmp_m);
        axpy(sq[i_cd], 1.0, vtmp_p);
      } // omp parallel
    }
  }

  timer_prop.stop();
  timer_prop.report();


  //- meson correlators
  Corr2pt_4spinor corr(gmset.get());
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  vout.general(vl, "\n");
  vout.general(vl, "2-point correlator:\n");

  double result = corr.meson_all(sq, sq);

  timer.stop();
  timer.report();

  RandomNumberManager::finalize();

  vout.general(vl, " corelib code.\n");


  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }

}


//============================================================END=====
