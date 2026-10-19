/*!
      @File    test_alt_Accel.cpp
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
      @version $LastChangedRevision: 2513 $
*/

#ifdef USE_ALT_ACCEL

#include "lib_alt_Accel/bridge_alt_accel.h"
#include "testlist_alt.h"
#include "spectrum_Domainwall_alt-tmpl.h"

#define IMPL ACCEL


namespace Test_alt_Accel{
  int test_Spectrum_Domainwall();
}

//====================================================================
int Test_alt_Accel::test_all()
{
  int result = 0;
  result += test_Spectrum_Domainwall();

  vout.general("Test result for Accel = %d\n", result);
  return result;
}


//====================================================================
int Test_alt_Accel::test_Spectrum_Domainwall()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file
     = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  int result = 0;
  //  result += test_dw.hadron_2ptFunction(test_file, "double");
  //  result += test_dw.hadron_2ptFunction(test_file, "float");
  //  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  //  result += test_dw.hadron_2ptFunction(test_file, "float_eo");

  result += test_dw.hadron_2ptFunction(test_file, "richardson_eo");

  return result;
}

//====================================================================

#endif
//============================================================END=====
