// @(#)root/lvmini:$Id$
// Author: A. Burgmeier

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013  DESY                                           *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class LVMiniMinimizer

#include "LVMini/LVMiniMinimizer.h"

#include "Math/IFunction.h"
#include "Math/IOptions.h"

// Declare the fortran symbols
extern "C" {
   extern void lvmeps_(float* eps, float* wlf1, float* wlf2);
   extern int lvmdim_(int* npar, int* mvec);
   extern void lvmini_(int* npar, int* mvec, int* nfcn, double* aux);
   extern void lvmfun_(double* x, double* f, int* iret, double* aux);
   extern int lvmind_(int* iarg);
}

namespace ROOT { 

namespace LVMini {

LVMiniMinimizer::LVMiniMinimizer():
   fCalcErrors(true), fFunc(NULL), fAux(NULL), fMin(0.), fIterations(0)
{
   common_init();
}

LVMiniMinimizer::LVMiniMinimizer(const char* type):
   fCalcErrors(true), fFunc(NULL), fAux(NULL), fMin(0.), fIterations(0)
{
   common_init();
}

void LVMiniMinimizer::common_init()
{
   // Default parameters:
   // TODO: Make these configurable
   float eps = 1e-4;
   float wlf1 = 1e-4;
   float wlf2 = 0.9;
   lvmeps_(&eps, &wlf1, &wlf2);
}

LVMiniMinimizer::~LVMiniMinimizer()
{
   delete[] fAux;
}

void LVMiniMinimizer::SetFunction(const ROOT::Math::IMultiGenFunction& func)
{
   std::cerr << "LVMiniMinimizer::SetFunction: LVMini needs gradient, use LVMiniMinimizer::SetFunction(const IMultiGradFunction&) instead!" << std::endl;
   std::cerr << "If you are using TH1::Fit, use the fitting option \"G\", and, ideally, provide an analytic gradient calculation" << std::endl;
}

void LVMiniMinimizer::SetFunction(const ROOT::Math::IMultiGradFunction& func)
{
   fFunc = &func;
}

bool LVMiniMinimizer::SetVariable(unsigned int var, const std::string& varname, double start, double step)
{
   if(var != fVariables.size())
   {
      std::cerr << "LVMiniMinimizer::SetVariable: Wrong index! Variables need to be added consecutively, use index=" << fVariables.size() << std::endl;
      return false;
   }
   else
   {
      fVariables.push_back(start);
      fVariableNames.push_back(varname);
      return true;
   }
}

bool LVMiniMinimizer::Minimize()
{
   if(fVariables.empty())
   {
      std::cerr << "LVMiniMinimizer::Minimize: No free variables!" << std::endl;
      return false;
   }

   if(!fFunc)
   {
      std::cerr << "LVMiniMinimizer::Minimize: No function to minimize provided!" << std::endl;
      return false;
   }

   int npar = fVariables.size();
   if(fDebug) npar = -npar;

   // I haven't seen a clear guide on how many vector pairs
   // should be used. This is a rough estimate from table 1 of
   // Blobel's manual and the recommendation to use values between
   // 6 and 29.
   int mvec = std::min(std::max(6, npar / 5), 29);
   if(fCalcErrors) mvec = -mvec;

   int mdim = lvmdim_(&npar, &mvec);

   delete[] fAux;
   fAux = new double[mdim];

   int nfcn = fMaxCalls;
   lvmini_(&npar, &mvec, &nfcn, fAux);

   unsigned int fMaxIterations = fMaxIter;
   if(fMaxIterations == 0) fMaxIterations = 10000;

   int iret = -1;
   fIterations = 0;

   do
   {
      fFunc->FdF(&fVariables[0], fMin, &fAux[0]);
      lvmfun_(&fVariables[0], &fMin, &iret, &fAux[0]);
   } while(++fIterations < fMaxIterations && iret < 0);
   fStatus = iret;

   if(iret == -1)
   {
      // Max number of iterations reached
      fValidError = false;
      std::cerr << "LVMiniMinimizer::Minimize: Max iterations reached" << std::endl;
      return false;
   }
   else if(iret > 0)
   {
      fValidError = false;
      std::cerr << "LVMiniMinimizer::Minimize: Minimization failed, error code=" << iret << std::endl;
      return false;
   }
   else
   {
      if(fCalcErrors)
         fValidError = true;

      // Minimization complete
      return true;
   }
}

double LVMiniMinimizer::MinValue() const
{
   int iarg = 0;
   int ind = lvmind_(&iarg);
   return fAux[ind - 1];
}

double LVMiniMinimizer::Edm() const
{
   return 0.0; // Can we obtain it? Maybe remember the difference between the last two evalutions?
}

const double* LVMiniMinimizer::X() const
{
   int iarg = 1;
   int ind = lvmind_(&iarg);
   return &fAux[ind];
}

const double* LVMiniMinimizer::MinGradient() const
{
   // TODO: Evaluate the gradient function at the minimum, and store the results somewhere...
   return NULL;
}

unsigned int LVMiniMinimizer::NCalls() const
{
   int iarg = -2;
   return lvmind_(&iarg);
}

unsigned int LVMiniMinimizer::NDim() const
{
   return fVariables.size();
}

unsigned int LVMiniMinimizer::NFree() const
{
   return fVariables.size();
}

bool LVMiniMinimizer::ProvidesError() const
{
   return fCalcErrors;
}

const double* LVMiniMinimizer::Errors() const
{
   int iarg = 2;
   if(fCalcErrors) iarg = 3;

   int ind = lvmind_(&iarg);
   if(ind == 0) return NULL; // should not happen

   return &fAux[ind];
}

double LVMiniMinimizer::CovMatrix(unsigned int i, unsigned int j) const
{
   int iarg = 5;
   int ind = lvmind_(&iarg);
   if(ind == 0) return -1.;

   // Matrix is symmetric
   if(j > i) std::swap(i, j);

   const unsigned int ijInd = ind - 1 + (j + 1) + ((i + 1) * (i + 1) - (i + 1)) / 2;
   return fAux[ijInd];
}

} // end namespace LVMini

} // end namespace ROOT

