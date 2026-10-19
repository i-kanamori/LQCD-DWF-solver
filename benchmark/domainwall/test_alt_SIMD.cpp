/*!
      @File    test_alt_SIMD.cpp
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: $
      @date    $LastChangedDate:  $
      @version $LastChangedRevision: 2317 $
*/

#ifdef USE_ALT_SIMD

#include "lib_alt_SIMD/bridge_alt_simd.h"
#include "testlist_alt.h"
#include "spectrum_Domainwall_alt-tmpl.h"

#define IMPL SIMD

namespace Test_alt_SIMD{
  int test_Spectrum_Domainwall_5din();
}

//====================================================================
int Test_alt_SIMD::test_all()
{
  int result = 0;
  result += test_Spectrum_Domainwall_5din();

  vout.general("Test result for SIMD = %d\n", result);
  return result;
}



//====================================================================
int Test_alt_SIMD::test_Spectrum_Domainwall_5din()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  int result = 0;
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
