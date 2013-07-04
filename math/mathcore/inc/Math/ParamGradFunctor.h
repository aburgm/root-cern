// @(#)root/mathcore:$Id$
// Author: A. Burgmeier Thu Jul 04 9:59:04 2013

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013  DESY                                           *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for Functor classes. 
// design is inspired by the Loki Functor

#ifndef ROOT_Math_ParamGradFunctor
#define ROOT_Math_ParamGradFunctor

// #ifndef ROOT_Math_IFunction
// #include "Math/IFunction.h"
// #endif

// #ifndef Root_Math_StaticCheck
// #include "Math/StaticCheck.h"
// #endif

//#ifndef __CINT__
//#include <memory> 

#include <vector>
#include <iostream>

namespace ROOT { 

namespace Math { 


/** class defining the signature for multi-dim parametric gradient functions 

   @ingroup  ParamGradFunctor_int
 */

class ParamGradFunctionBase { 
  public: 
   virtual ~ParamGradFunctionBase() {}
//   virtual double operator() (const double * x, const double *p) const = 0; 
   virtual double operator() (double * x, double *p, double *u) = 0; 
   virtual ParamGradFunctionBase * Clone() const = 0; 
};



/** 
   ParamGradFunctor Handler class is responsible for wrapping any other functor and pointer to 
   free C functions.
   It can be created from any function implementing the correct signature 
   corresponding to the requested type

   @ingroup  ParamGradFunctor_int

*/ 
#ifndef __CINT__

template<class ParentFunctor, class Func >
class ParamGradFunctorHandler : public ParentFunctor::Impl { 

   typedef typename ParentFunctor::Impl Base; 

public: 

   // constructor 
   ParamGradFunctorHandler(const Func & fun) : fFunc(fun) {}


   virtual ~ParamGradFunctorHandler() {}


   // for 1D functions
   inline double operator() (double x, double *p, double *u)  { 
      return fFunc(x,p,u); 
   }  
//    inline double operator() (double x, const double *p, double *u) const { 
//       return fFunc(x,p,u); 
//    }  
   // for multi-dimensional functions
//    inline double operator() (const double * x, const double *p, double *u) const { 
//       return fFunc(x,p,u); 
//    }  
   inline double operator() (double * x, double *p, double *u)  { 
      return FuncEvaluator<Func>::Eval(fFunc,x,p,u); 
   }  

   // clone (use same pointer)
   ParamGradFunctorHandler  * Clone() const { 
      return new ParamGradFunctorHandler(fFunc); 
   } 


private :
   
   Func fFunc; 

   // structure to distinguish pointer types
   template <typename F> struct FuncEvaluator { 
      inline static double Eval( F & f, double *x, double * p, double *u) { 
         return f(x,p,u);
      }
   };
   template <typename F> struct FuncEvaluator<F*> { 
      inline static double Eval( F * f, double *x, double * p, double *u) { 
         return (*f)(x,p,u);
      }
   };
   template <typename F> struct FuncEvaluator<F* const> { 
      inline static double Eval( const F * f, double *x, double * p, double *u) { 
         return (*f)(x,p,u);
      }
   };
   // need maybe also volatile ? 
};


#if defined(__MAKECINT__) || defined(G__DICTIONARY) 
// needed since CINT initialize it with TRootIOCtor
//class TRootIOCtor; 
template<class ParentFunctor> 
class ParamGradFunctorHandler<ParentFunctor,TRootIOCtor *> : public ParentFunctor::Impl 
{
public:

   ParamGradFunctorHandler(TRootIOCtor  *) {}

   double operator() (double *, double * )  { return 0; } 
   // clone (use same pointer)
   ParamGradFunctorHandler  * Clone() const { 
      return 0; 
   } 

}; 
#endif   


/**
   ParamGradFunctor Handler to Wrap pointers to member functions 

   @ingroup  ParamGradFunctor_int
*/
template <class ParentFunctor, typename PointerToObj,
          typename PointerToMemFn>
class ParamGradMemFunHandler : public ParentFunctor::Impl
{
   typedef typename ParentFunctor::Impl Base;

   
public:
   
   /// constructor from a pointer to the class and a pointer to the function
   ParamGradMemFunHandler(const PointerToObj& pObj, PointerToMemFn pMemFn) 
      : fObj(pObj), fMemFn(pMemFn)
   {}

   virtual ~ParamGradMemFunHandler() {}
        
//    inline double operator() (double x, const double * p, double *u) const { 
//       return ((*fObj).*fMemFn)(x,p);  
//    }  

   inline double operator() (double x, double * p, double *u)  { 
      return ((*fObj).*fMemFn)(x,p,u);  
   }  
       
//    inline double operator() (const double * x, const double * p, double *u) const { 
//       return ((*fObj).*fMemFn)(x,p,u);  
//    }

   inline double operator() (double * x, double * p, double *u)  { 
      return ((*fObj).*fMemFn)(x,p,u);  
   }  

   // clone (use same pointer)
   ParamGradMemFunHandler  * Clone() const { 
      return new ParamGradMemFunHandler(fObj, fMemFn); 
   } 


private :
   ParamGradMemFunHandler(const ParamGradMemFunHandler&); // Not implemented
   ParamGradMemFunHandler& operator=(const ParamGradMemFunHandler&); // Not implemented
       
   PointerToObj fObj;
   PointerToMemFn fMemFn; 

};

#endif  



/**
   Param Functor class for Multidimensional parametric gradient functions. 
   It is used to wrap in a very simple and convenient way 
   any other C++ callable object (implemention double operator( const double *, const double *, double *) ) 
   or a member function with the correct signature, 
   like Foo::EvalPar(const double *, const double *, double *)

   @ingroup  ParamFunc

 */


class ParamGradFunctor   { 


public: 

   typedef  ParamGradFunctionBase Impl;   


   /** 
      Default constructor
   */ 
   ParamGradFunctor ()  : fImpl(0) {}  


   /** 
       construct from a pointer to member function (multi-dim type)
    */ 
   template <class PtrObj, typename MemFn>
   ParamGradFunctor(const PtrObj& p, MemFn memFn)
      : fImpl(new ParamMemFunHandler<ParamGradFunctor, PtrObj, MemFn>(p, memFn))
   {}



   /**
      construct from another generic Functor of multi-dimension 
    */
   template <typename Func> 
   explicit ParamGradFunctor( const Func & f) : 
      fImpl(new ParamGradFunctorHandler<ParamGradFunctor,Func>(f) )
   {}



   // specialization used in TF1
   typedef double (* FreeFunc ) (double * , double *, double *);
   ParamGradFunctor(FreeFunc f) : 
      fImpl(new ParamGradFunctorHandler<ParamGradFunctor,FreeFunc>(f) )
   {
   }


   /** 
      Destructor (no operations)
   */ 
   virtual ~ParamGradFunctor ()  {
      if (fImpl) delete fImpl;
   }  

   /** 
      Copy constructor
   */ 
   ParamGradFunctor(const ParamGradFunctor & rhs) : 
      fImpl(0)
   {
//       if (rhs.fImpl.get() != 0) 
//          fImpl = std::auto_ptr<Impl>( (rhs.fImpl)->Clone() ); 
      if (rhs.fImpl != 0)  fImpl = rhs.fImpl->Clone(); 
   } 

   /** 
      Assignment operator
   */ 
   ParamGradFunctor & operator = (const ParamGradFunctor & rhs)  {
//      ParamGradFunctor copy(rhs); 
      // swap auto_ptr by hand
//       Impl * p = fImpl.release(); 
//       fImpl.reset(copy.fImpl.release());
//       copy.fImpl.reset(p);

      
      if (fImpl) delete fImpl;
      fImpl = 0; 
      if (rhs.fImpl != 0) 
         fImpl = rhs.fImpl->Clone();

      return *this;
   }

   void * GetImpl() { return (void *) fImpl; }


   double operator() (double * x, double * p, double *u)  { 
      return (*fImpl)(x,p,u); 
   }  



   bool Empty() { return fImpl == 0; }


   void SetFunction(Impl * f) { 
      fImpl = f;
   }

private :


   //std::auto_ptr<Impl> fImpl; 
   Impl * fImpl; 


}; 



   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_ParamGradFunctor */
