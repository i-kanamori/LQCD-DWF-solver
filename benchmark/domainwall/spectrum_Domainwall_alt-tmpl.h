/*!
      @File    spectrum_Domainwall_alt.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2022-01-30 02:14:34 #$
      @version $LastChangedRevision: 2344 $
*/

#include "spectrum_Domainwall_alt.h"
#include "test.h"

#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/fopr_Domainwall.h"
#include "lib/Fopr/fopr_Domainwall_eo.h"
//#include "lib/Measurements/Fermion/fprop_Standard_lex.h"
//#include "lib/Measurements/Fermion/fprop_Standard_eo.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Solver/solver.h"
#include "lib/Measurements/Fermion/source.h"
//#include "lib/Measurements/Fermion/fprop_Standard_Precond.h"
#include "lib/Measurements/Fermion/corr2pt_4spinor.h"


//#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"
//#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex_Mixedprec.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
//#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Mixedprec.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"

#include "lib/fapp_macros.h"
#ifdef USE_MPREAL
#include "bridge_mpreal.h"
#endif

// class name
template<Impl IMPL>
const std::string Spectrum_Domainwall_alt<IMPL>::class_name =
                                    "Spectrum_Domainwall_alt";

#ifdef USE_ALT_SIMD
#define USE_CORELIB_WILSON
#endif

//====================================================================
template<Impl IMPL>
void Spectrum_Domainwall_alt<IMPL>::init()
{
  // do nothing.
}

//====================================================================
template<Impl IMPL>
int Spectrum_Domainwall_alt<IMPL>::hadron_2ptFunction(
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
  //Parameters params_fopr_ref = params_all.lookup("Fopr_ref");
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

  // configuration setup.
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // for correlation function
  string gmset_type  = params_fopr.get_string("gamma_matrix_type");
  vout.general(vl, "  gmset_type   = %s\n", gmset_type.c_str());
  unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(gmset_type));

  unique_ptr<AFopr<AField<double, IMPL> > >
    afopr(AFopr<AField<double, IMPL> >::New(str_fopr_type, params_fopr));
  afopr->set_config(U.get());

  unique_ptr<Source> source(Source::New(str_source_type));
  source->set_parameters(params_source);

  unique_ptr<Timer> timer(new Timer(test_name));

  // setup fermion propgator computer.
  Parameters params_solver = params_all.lookup("Solver");

  //  unique_ptr<Solver> solver; // only used for original Fprop.

  //  check(fopr, params_fopr);
  vout.general(vl, "  check finished\n\n");
  //return 0;

  //  unique_ptr<Fprop> fprop;
  unique_ptr<Fprop> fprop;

  if(mode == "org"){
    //    fprop.reset(new Fprop_Standard_Precond(fopr.get()));
    //    fprop->set_mode("D");
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    abort();
  }else if(mode == "double"){
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    abort();
    //    fprop.reset(new Fprop_alt_Standard_lex<AField<double,IMPL> >(
    //                                     params_fopr, params_solver) );
    //    fprop->set_mode("D");

  }else if(mode == "double_prec"){
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    exit(EXIT_FAILURE);
    //    fprop.reset(new Fprop_alt_Standard_lex<AField<double,IMPL> >(
    //                                     params_fopr, params_solver) );
    //    fprop->set_mode("D_prec");

  }else if(mode == "double_eo"){
    fprop.reset(new Fprop_alt_Standard_eo<AField<double,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");

  }else if(mode == "float"){
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    exit(EXIT_FAILURE);
    //    params_solver.set_double("convergence_criterion_squared", 1.0e-14);
    //    fprop.reset(new Fprop_alt_Standard_lex<AField<float,IMPL> >(
    //                                     params_fopr, params_solver) );
    //    fprop->set_mode("D");

  }else if(mode == "float_eo"){
    params_solver.set_double("convergence_criterion_squared", 1.0e-14);
    fprop.reset(new Fprop_alt_Standard_eo<AField<float,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");
  }else if(mode == "richardson_eo"){
    fprop.reset(new Fprop_alt_Standard_eo_Richardson<AField<double,IMPL>,
                                                    AField<float,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");
#ifdef USE_MPREAL
  }else if(mode == "richardson_eo_mpreal"){
    int prec=params_solver.get_int("mpreal_prec");
    mpreal::set_default_prec(prec);
    fprop.reset(new Fprop_alt_Standard_eo_Richardson<AField<double,IMPL>,
                                                    AField<mpreal,IMPL> >(
                                     params_fopr, params_solver) );
    fprop->set_mode("D");
#endif
  }else if(mode == "mixed_eo"){
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    exit(EXIT_FAILURE);
    //    fprop.reset(new Fprop_alt_Standard_eo_Mixedprec<AField<double,IMPL>,
    //                                                    AField<float,IMPL> >(
    //                                     params_fopr, params_solver) );
    //    fprop->set_mode("D");
  }else if(mode == "mixed"){
    vout.crucial("mode = %s is temprary unsupported.\n", mode.c_str());
    exit(EXIT_FAILURE);
    //fprop.reset(new Fprop_alt_Standard_lex_Mixedprec<AField<double,IMPL>,
    //                                                    AField<float,IMPL> >(
    //                                     params_fopr, params_solver) );
    //    fprop->set_mode("D");
  }else{
    vout.crucial(vl, "%s: irrelevant mode =%s\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  fprop->set_config(U.get());

  if(mode=="double_prec" || mode == "float_prec"){
    //    fprop->mult_performance("DdagD_prec", 100);
    fprop->mult_performance("DdagD_prec", 10);
  } else {
    //    fprop->mult_performance("DdagD", 100);
    //    fprop->mult_performance("D", 100);
    fprop->mult_performance("D", 10);
    fprop->mult_performance("DdagD", 10);
  }

  //return 0;
  // setup source and progator.
  int Ns = params_fopr.get_int("extent_of_5th_dimension");
  params_fopr.fetch_int("extent_of_5th_dimension", Ns);

  Field_F b, vt1, vt2;
  Field_F b5(Nvol, Ns), x5(Nvol, Ns), y5(Nvol, Ns);
  std::vector<Field_F> sq(Nc * Nd);

  // kernel of Domainwall operator
  typedef AField<double,IMPL> AFIELD;
  std::string kernel_type;
  double M0;
  double coeff_c;
  int err = params_fopr.fetch_string("kernel_type", kernel_type);
  if(err > 0){
    vout.crucial(vl, "%s: Error: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("domain_wall_height", M0);
  if(err > 0){
    vout.crucial(vl, "Error at %s: domain_wall_height is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  err += params_fopr.fetch_double("coefficient_c", coeff_c);
  if(err > 0){
    vout.crucial(vl, "Error at %s: coefficient_c is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel = params_fopr;
  double kappa = 1.0/(8.0 - 2.0 * M0);
  double one_over_2kappa=4.0-M0;
  params_kernel.set_double("hopping_parameter", kappa);

  // Factory is assumed to work
#ifdef USE_CORELIB_WILSON
  unique_ptr<Fopr> foprw(Fopr::New(kernel_type, params_kernel));
#else
  unique_ptr<AFopr<AFIELD> > foprw(AFopr<AFIELD>::New(kernel_type, params_kernel));
#endif
  foprw->set_mode("D");
  foprw->set_config(U.get());

  // afields
  AIndex_lex<double,IMPL> index_alt;

  int nin=b.nin();

  AFIELD ab5, ax5, ay5;
  //  AFIELD ab5(nin, Nvol, Ns);
  //  AFIELD ax5(nin, Nvol, Ns);
  //  AFIELD ay5(nin, Nvol, Ns);
  if(str_fopr_type == "Domainwall_5din"){
    ab5.reset(nin*Ns, Nvol, 1);
    ax5.reset(nin*Ns, Nvol, 1);
    ay5.reset(nin*Ns, Nvol, 1);
  } else {
    ab5.reset(nin, Nvol, Ns);
    ax5.reset(nin, Nvol, Ns);
    ay5.reset(nin, Nvol, Ns);

  }
#ifdef USE_CORELIB_WILSON
  //  Field_F vt1;
  //  Field_F vt2;
#else
  AFIELD ab(nin, Nvol, 1);
  AFIELD avt1(nin, Nvol, 1);
  AFIELD avt2(nin, Nvol, 1);
#endif
  
#pragma omp parallel
  {
    b.set(0.0);
  }
  source->set(b, 0);

  // Solver main part.
  timer->start();

  int nconv;
  double diff;

  Timer timer_prop("calculating porpagator");
  timer_prop.reset();
  timer_prop.start();

  vout.general(vl, "\n");
  vout.general(vl, "Solving quark propagator:\n");
  vout.general(vl, "  color spin   Nconv      diff           diff2\n");

  for (int ispin = 0; ispin < Nd; ++ispin) {
    for (int icolor = 0; icolor < Nc; ++icolor) {
      //  for (int ispin = 0; ispin < 1; ++ispin) {
      //    for (int icolor = 0; icolor < 1; ++icolor) {
      int idx = icolor + Nc * ispin;
      source->set(b, idx);

      // set 5d source as
      //   D_- ( P_+ 0 ....0 P_-)^T b4
#pragma omp parallel
      {
#ifdef USE_CORELIB_WILSON
        b5.set(0.0);
        foprw->mult_gm5(vt1, b);
        axpy(vt1, +1.0, b);    //   (1+gm5)b
        foprw->mult(vt2, vt1);
        aypx(-one_over_2kappa*coeff_c, vt2, vt1); // (-cD+1) (1+gm5)bb
        scal(vt2, -0.5);          // 0.5*(cD-1) (1+gm5)b
        copy(b5, 0, vt2, 0);

        foprw->mult_gm5(vt1,b);
        axpy(vt1, -1.0, b);    //   (-1+gm5)b  [ = -(1-gm5)b ]
        foprw->mult(vt2, vt1);
        aypx(-coeff_c*one_over_2kappa, vt2, vt1); // (cD-1) (1-gm5)b
        scal(vt2, 0.5);          // 0.5*(cD-1) (1-gm5)b
        copy(b5, Ns-1, vt2, 0);
#else
        if(foprw->needs_convert()){
          foprw->convert(ab, b);
        }else{
          convert_spinor(index_alt, ab, b);
        }
#pragma omp barrier

        // s5=0
        foprw->mult_gm5(avt1,ab);
        axpy(avt1, +1.0, ab);    //   (1+gm5)ab
        foprw->mult(avt2, avt1);
        aypx(-one_over_2kappa*coeff_c, avt2, avt1); // (-cD+1) (1+gm5)ab
        scal(avt2, -0.5);          // 0.5*(cD-1) (1+gm5)ab
        if(foprw->needs_convert()){
          foprw->reverse(vt1, avt2);
        }else{
          reverse_spinor(index_alt, vt1, avt2);
        }

        foprw->mult_gm5(avt1,ab);
        axpy(avt1, -1.0, ab);    //   (-1+gm5)ab  [ = -(1-gm5)ab ]
        foprw->mult(avt2, avt1);
        aypx(-coeff_c*one_over_2kappa, avt2, avt1); // (cD-1) (1-gm5)ab
        scal(avt2, 0.5);          // 0.5*(cD-1) (1-gm5)ab
        if(foprw->needs_convert()){
          foprw->reverse(vt2, avt2);
        }else{
          reverse_spinor(index_alt, vt2, avt2);
        }
#pragma omp barrier
        // build 5dim vector
        b5.set(0.0);
#pragma omp barrier
        copy(b5, 0, vt1, 0);
        copy(b5, Ns-1, vt2, 0);
#endif
      }

      START_FAPP("fprop_invert", 1, 0);
      fprop->invert(x5, b5, nconv, diff);
      STOP_FAPP("fprop_invert", 1, 0);

#pragma omp parallel
      {
        {// check the solution
          if(afopr->needs_convert()){
            afopr->convert(ab5, b5);
            afopr->convert(ax5, x5);
          }else{
            convert_spinor(index_alt, ab5, b5);
            convert_spinor(index_alt, ax5, x5);
          }
          copy(ay5, ab5);
          afopr->set_mode("D");
          afopr->mult(ay5, ax5);
          axpy(ay5, -1.0, ab5);
          double diff2= ay5.norm2() / ab5.norm2();
          vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                       icolor, ispin, nconv, diff, diff2);

        }

        // convert to 4D propagator
#ifdef USE_CORELIB_WILSON
        copy(sq[idx], 0, x5, 0);
#pragma omp barrier
        foprw->mult_gm5(vt1, sq[idx]);
        axpy(sq[idx], -1.0, vt1);      // sq[idx] = (1-gm5)x5[0]
        copy(vt2, 0, x5, Ns-1);
#pragma omp barrier
        foprw->mult_gm5(vt1, vt2);
        axpy(sq[idx], 1.0, vt1);
        axpy(sq[idx], 1.0, vt2);       // sq[idx] = (1-gm5)x5[0] + (1+gm5)x5[Ns-1]
        scal(sq[idx], 0.5);
#else
        copy(vt1, 0, x5, 0);
        copy(vt2, 0, x5, Ns-1);

        if(foprw->needs_convert()){
          foprw->convert(ab, vt1);
          foprw->convert(avt2, vt2);
        }else{
          convert_spinor(index_alt, ab, vt1);
          convert_spinor(index_alt, avt2, vt2);
        }
#pragma omp barrier
        foprw->mult_gm5(avt1, ab);
        axpy(ab, -1.0, avt1);      // ab = (1-gm5)vt1 converted
#pragma omp barrier
        foprw->mult_gm5(avt1, avt2);
        axpy(ab, 1.0, avt1);
        axpy(ab, 1.0, avt2);       // ab = (1-gm5)vt1 + (1+gm5)vt2 converted
        scal(ab, 0.5);
#pragma omp barrier
        if(foprw->needs_convert()){
          foprw->reverse(sq[idx], ab);
        }else{
          reverse_spinor(index_alt, sq[idx], ab);
        }
#endif
      }
    }
  }

  fprop->report_performance();
  timer_prop.stop();
  timer_prop.report();

  //- meson correlators
  vout.general(vl, "skippling meson correlators\n");

  Corr2pt_4spinor corr(gmset.get());
  corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));

  vout.general(vl, "\n");
  vout.general(vl, "2-point correlator:\n");

  double result = corr.meson_all(sq, sq);

  //  double result = 0.0;
  timer->report();

  RandomNumberManager::finalize();

  vout.general(vl, " alt code.\n");


  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }

}


//============================================================END=====
