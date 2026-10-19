/*
        @file    testlist_accel.h
        @brief
        @author  Hideo Matsufuru
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2013-04-08 18:00:27 #$
        @version $LastChangedRevision: 2294 $
*/

#ifndef TESTLIST_ALT_INCLUDED
#define TESTLIST_ALT_INCLUDED


namespace Test_alt_Corelib {
  int test_all();
}

#ifdef USE_ALT_SIMD
namespace Test_alt_SIMD {
  int test_all();
}
#endif

#ifdef USE_ALT_ACCEL
namespace Test_alt_Accel {
  int test_all();
}
#endif

#ifdef USE_ALT_QXS
namespace Test_alt_QXS {
  int test_all();
  extern bool use_mpfr;
}
#endif


#endif
