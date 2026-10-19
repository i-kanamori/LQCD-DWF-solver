/*!
         @file    main.cpp
         @brief
         @author  Hideo Matsufuru (matufuru)
                  $LastChangedBy: kanamori $
         @date    $LastChangedDate:: 2021-07-25 17:30:02 #$
         @version $LastChangedRevision: 2294 $
*/

#include <vector>
#include <string>
using namespace std;


//#include "mpi.h"
//#include "bridge.h"

// corelib
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Parameters/parameterManager_YAML.h"
#include "lib/Tools/timer.h"
#include "lib/bridge_init_factory.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/timer_global.h"

// alt-code
#ifdef USE_ALT_CODE
#include "lib_alt/bridge_alt_init.h"
#endif

// local
#include "testlist_alt.h"
#include "show_version.h"

const std::string filename_main_input = "main.yaml";
// const std::string filename_main_input = "stdin";

//- prototype declarations
//int run_test();

//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  Communicator::init(&argc, &argv);

  // ####  banner  ####
  vout.general(vl, "Bridge++ %s\n\n", BRIDGE_VERSION);
  Bridge::show_version();
  string target("all");

  if(argc>1){
    target=string(argv[1]);
    vout.general("target: %s\n", target.c_str());
  } else {
    vout.general("no target is given, assumes all.\n");
    vout.general("  available target:\n    corelib\n");
#ifdef USE_ALT_CODE
#ifdef USE_ALT_QXS
    vout.general("    alt_qxs\n");
#endif
#ifdef USE_ALT_ACCEL
    vout.general("    alt_accel\n");
#endif
#ifdef USE_ALT_SIMD
    vout.general("    alt_simd\n");
#endif

#endif
  }

  Parameters params_all  = ParameterManager::read(filename_main_input);
  Parameters params_main = params_all.lookup("Main");

  const vector<int> lattice_size = params_main.get_int_vector("lattice_size");
  const vector<int> grid_size    = params_main.get_int_vector("grid_size");
  const int number_of_thread     = params_main.get_int("number_of_thread");
  const int number_of_color      = params_main.get_int("number_of_color");
  const string str_logfile       = params_main.get_string("log_filename");
  const string str_ildg_logfile  = params_main.get_string("ildg_log_filename");
  const string str_vlevel        = params_main.get_string("verbose_level");


  //- initializations
  vl = vout.set_verbose_level(str_vlevel);
  CommonParameters::init_Vlevel(vl);

  if (str_logfile != "stdout") {
    vout.init(str_logfile);
  }

  if (str_ildg_logfile != "stdout") {
    vout.ildg_init(str_ildg_logfile);
  }


  // CommonParameters::init(lattice_size, grid_size);
  CommonParameters::init(lattice_size, grid_size, number_of_color);
  Communicator::setup();

  ThreadManager_OpenMP::init(number_of_thread);

  timer_global_init();

  //- print input parameters
  vout.general(vl, "Main: input parameters\n");
  vout.general(vl, "  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
  vout.general(vl, "  grid_size        = %s\n", Parameters::to_string(grid_size).c_str());

  vout.general(vl, "  number of thread = %d\n", number_of_thread);
  vout.general(vl, "  number of color  = %d\n", number_of_color);
  vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
  vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
  vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
  vout.general(vl, "\n");

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_logfile);
  err += ParameterCheck::non_NULL(str_ildg_logfile);

  if (err) {
    vout.crucial(vl, "Error at main: input parameters have not been set.\n");
    exit(EXIT_FAILURE);
  }

#ifdef USE_FACTORY
#ifdef USE_FACTORY_AUTOREGISTER
#else
  bridge_init_factory();
#endif
#endif


  // initialization of alt-code
#ifdef USE_ALT_CODE
  bridge_alt_init(params_main);
#endif

  //- timestamp (starting time)
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->timestamp();
  timer->start();

  //####  here function is called explicitly  ####

  vout.general("target: %s\n", target.c_str());
  if(target=="all" || target=="corelib"){
    Test_alt_Corelib::test_all();
  }
#ifdef USE_ALT_QXS
  else if(target=="all" || target=="alt_qxs"){
    Test_alt_QXS::test_all();
  }
  else if(target=="all" || target=="alt_qxs_mpfr"){
    Test_alt_QXS::use_mpfr = true;
    Test_alt_QXS::test_all();
  }
#endif
#ifdef USE_ALT_ACCEL
  else if(target=="all" || target=="alt_accel"){
    Test_alt_Accel::test_all();
  }
#endif
#ifdef USE_ALT_SIMD
  else if(target=="all" || target=="alt_simd"){
    Test_alt_SIMD::test_all();
  }
#endif
  else {
    vout.general("unknown target: %s\n", target.c_str());
  }


  timer_global_report();

  //- timestamp (end time)
  timer->report();
  timer->timestamp();

  // finalization of alt-code
#ifdef USE_ALT_CODE
  Communicator::sync();
  vout.general("calling bridge_alt_fin()\n");
  bridge_alt_fin();
  Communicator::sync();
#endif


  vout.general("calling ThreadManager_OpenMP::finalize()\n");
  ThreadManager_OpenMP::finalize();
  vout.general("calling ThreadManager_OpenMP::finalize(), done\n");
  fflush(stdout);
  Communicator::sync();
  Communicator::finalize();

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====
