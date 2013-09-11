// @(#)root/minuit2:$Id$
// Author: L. Moneta    10/2006  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FCNGradAdapter
#define ROOT_Minuit2_FCNGradAdapter

#ifndef ROOT_Minuit2_FCNGradientBase
#include "Minuit2/FCNGradientBase.h"
#endif

//#define DEBUG
#ifdef DEBUG
#include <iostream> 
#endif

namespace ROOT {

   namespace Minuit2 {

/** 


template wrapped class for adapting to FCNBase signature a IGradFunction

@author Lorenzo Moneta

@ingroup Minuit

*/

template< class Function> 
class FCNGradAdapter : public FCNGradientBase {

public:

   FCNGradAdapter(const Function & f, double up = 1.) : 
      fFunc(f) , 
      fUp (up)  
      //fGrad(std::vector<double>(fFunc.NDim() ) )

   {}

   ~FCNGradAdapter() {}

  
   double operator()(const std::vector<double>& v) const { 
      return fFunc.operator()(&v[0]); 
   }
   double operator()(const double *  v) const { 
      return fFunc.operator()(v); 
   }

   double Up() const {return fUp;}
   
   std::vector<double> Gradient(const std::vector<double>& v) const{
      std::vector<double> grad(v.size() );
      fFunc.Gradient(&v[0], &grad[0]);

#ifdef DEBUG
      std::cout << " gradient in FCNAdapter = { " ;
      for (unsigned int i = 0; i < grad.size(); ++i)
         std::cout << grad[i] << "\t";
      std::cout << "}" << std::endl;
#endif

      return grad;
   }
  
    void Gradient(const std::vector<double>& v, std::vector<double>& grad, std::vector<double>& g2) const {
       if (g2.size() > 0) {
          fFunc.Gradient(&v[0], &grad[0], &g2[0] );
          // case g2 is not implemented 
          if (g2[0] == std::numeric_limits<double>::quiet_NaN() ) g2.clear(); 
       }   
       else {
          fFunc.Gradient(&v[0], &grad[0]);
       }

#ifdef DEBUG
      std::cout << " gradient in FCNAdapter = { " ;
      for (unsigned int i = 0; i < grad.size(); ++i)
         std::cout << grad[i] << "\t";
      std::cout << "}" << std::endl;
#endif

   }
   // forward interface
   //virtual double operator()(int npar, double* params,int iflag = 4) const;
   bool CheckGradient() const { return false; } 

private:
   const Function & fFunc; 
   double fUp; 
   //mutable std::vector<double> fGrad;
};

   } // end namespace Minuit2

} // end namespace ROOT



#endif //ROOT_Minuit2_FCNGradAdapter
