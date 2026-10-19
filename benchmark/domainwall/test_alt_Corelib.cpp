/*!
      @File    test_Corelib.cpp
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-10-25 13:27:53 #$
      @version $LastChangedRevision: 2317 $
*/



#include "testlist_alt.h"
#include "spectrum_Domainwall_alt.h"

#define IMPL CORELIB

namespace Test_alt_Corelib{
  int test_Spectrum_Domainwall();
}

//====================================================================
int Test_alt_Corelib::test_all()
{
  int result = 0;
  result += test_Spectrum_Domainwall();

  vout.general("Test result for Corelib = %d\n", result);
  return result;
}

//====================================================================
int Test_alt_Corelib::test_Spectrum_Domainwall()
{
  Spectrum_Domainwall_alt<IMPL> test_dw;
  string test_file = "test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml";

  int result = 0;
  //result += test_dw.hadron_2ptFunction(test_file, "double");
  result += test_dw.hadron_2ptFunction(test_file, "double_eo");
  return result;
}


//============================================================END=====
