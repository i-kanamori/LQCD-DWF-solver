/*!
      @file    spectrum_alt.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "spectrum_alt.h"
#include "test.h"

//#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Tools/randomNumberManager.h"

#include "job_Utils.h"

// class name
const std::string Spectrum_alt::class_name = "Spectrum_alt";
//====================================================================
void Spectrum_alt::init()
{
  m_vl = CommonParameters::Vlevel();

  //  vout.general(m_vl, "hello, spectrum.\n");

  //  RandomNumberManager::initialize("Mseries", 1234567UL);

}

//====================================================================
void Spectrum_alt::tidyup()
{

  //  RandomNumberManager::finalize();

}


//====================================================================
void Spectrum_alt::setup_config(unique_ptr<Field_G>& U,
                                Parameters& params_test)
{
  string gconf_status = params_test.get_string("gauge_config_status");
  string gconf_read   = params_test.get_string("gauge_config_type_input");
  string readfile     = params_test.get_string("config_filename_input");

  int err = 0;
  err += ParameterCheck::non_NULL(gconf_status);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameters have not been set\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Configuration setup:\n");
  vout.general(m_vl, "  gconf_status = %s\n", gconf_status.c_str());
  vout.general(m_vl, "  gconf_read   = %s\n", gconf_read.c_str());
  vout.general(m_vl, "  readfile     = %s\n", readfile.c_str());

#if defined USE_GROUP_SU3

  if (gconf_status == "Continue") {
    GaugeConfig(gconf_read).read(*U, readfile);
  } else if (gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(*U);
  } else if (gconf_status == "Hot_start") {
    GaugeConfig("Random").read(*U);
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported gconf status \"%s\"\n",
                 class_name.c_str(), gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

#else

  vout.general(m_vl, "  Nun-SU3 gauge group.\n");

  unique_ptr<FieldIO> gconf;
  if(gconf_read == "Binary"){
    gconf.reset(new FieldIO_Binary(IO_Format::Trivial));
  }else if(gconf_read == "Text"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else if(gconf_read == "Text_4x4x4x8"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else if(gconf_read == "NO_OUTPUT"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else{
    vout.crucial(m_vl, "  read config type not supported.\n");
    exit(EXIT_FAILURE);
  }

  if (gconf_status == "Continue") {
    gconf->read_file(U.get(), readfile);
  } else if (gconf_status == "Cold_start") {
    U->set_unit();
  } else {
    vout.crucial(m_vl, "Error at Test_Spectrum: unsupported gconf"
                       " status \"%s\"\n", gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

#endif

  //  Staple_lex staple;
  //  double plaq = staple.plaquette((Field_G)*U);
  //  vout.general(m_vl, "  plaq = %f\n", plaq);

  return;

}


//============================================================END=====
