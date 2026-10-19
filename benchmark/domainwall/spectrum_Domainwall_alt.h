/*!
        @file    $Id: spectrum_Domainwall_alt.h #$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$

        @version $LastChangedRevision: 2160 $
*/

#ifndef SPECTRUM_DOMAINWALL_ALT_INCLUDED
#define SPECTRUM_DOMAINWALL_ALT_INCLUDED

#include <vector>
#include <string>

#include  "spectrum_alt.h"

#ifdef USE_ALT_CODE
#include  "lib_alt/alt_impl.h"
#else
enum Impl
{
  CORELIB
};
#endif

#include  "lib/Parameters/commonParameters.h"
#include  "lib/Parameters/parameters.h"
#include  "lib/Field/field_G.h"
#include  "lib/Fopr/fopr.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;


template<Impl IMPL>
class Spectrum_Domainwall_alt : public Spectrum_alt
{

public:
  static const std::string class_name;

private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

public:

  //! constructor
  Spectrum_Domainwall_alt()
    : Spectrum_alt(), m_vl(CommonParameters::Vlevel())
    { init(); }

  //! destructor
  ~Spectrum_Domainwall_alt(){}

  //! calculation of hadron two-point correlators.
  int hadron_2ptFunction(std::string test_file, std::string mode);

 private:

  //! initial setup
  void init();

  //! check operator
  void check(unique_ptr<Fopr>&, Parameters&);

};





template<>
class Spectrum_Domainwall_alt<CORELIB> : public Spectrum_alt
{

public:
  static const std::string class_name;

private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

public:

  //! constructor
  Spectrum_Domainwall_alt()
    : Spectrum_alt(), m_vl(CommonParameters::Vlevel())
    { init(); }

  //! destructor
  ~Spectrum_Domainwall_alt(){}

  //! calculation of hadron two-point correlators.
  int hadron_2ptFunction(std::string test_file, std::string mode);

 private:

  //! initial setup
  void init();

};

#endif // SPECTRUM_DOMAINWALL_ALT_INCLUDED
