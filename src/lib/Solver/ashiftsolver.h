/*!
        @file    ashiftsolver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef ASHIFTSOLVER_INCLUDED
#define ASHIFTSOLVER_INCLUDED

#include "bridge_defs.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

//! Shiftsolver class as an abstract base class for multi-shift solvers.
template<typename FIELD>
class AShiftsolver
{
 public:
  AShiftsolver() {}

  virtual ~AShiftsolver() {}

 private:
  // non-copyable
  AShiftsolver(const AShiftsolver&);
  AShiftsolver& operator=(const AShiftsolver&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void solve(std::vector<FIELD>& solution,
                     const std::vector<double>& shift,
                     const FIELD& source,
                     int& Nconv,
                     double& diff) = 0;

  virtual double flop_count() = 0;
};

#endif
