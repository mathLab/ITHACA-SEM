///////////////////////////////////////////////////////////////////////////////
//
// File Arpack.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: wrapper of functions around standard Arpack routines.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_ARPACK_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_ARPACK_HPP

#include <LibUtilities/LinearAlgebra/TransF77.hpp>

// Translations for using Fortran version of arpack
namespace Arpack
{
    extern "C"
    {
		//ARPACK ROUTINES
		
		//Computation of eigenvalues

		void F77NAME(dsaupd) (int& ido,           const char* bmat,
							  const int& n,       const char* which,
							  const int& nev,     const double& tol,
							  double* resid,      const int& ncv,
							  double* v,          const int& ldv,
							  int* iparam,        int* ipntr,         
							  double* workd,      double* workl,      
							  const int& lworkl,  int& info );
		
		void F77NAME(dnaupd) (int& ido,           const char* bmat,
							  const int& n,       const char* which,
							  const int& nev,     const double& tol,
							  double* resid,      const int& ncv,
							  double* v,          const int& ldv,
							  int* iparam,        int* ipntr,         
							  double* workd,      double* workl,      
							  const int& lworkl,  int& info);
		

		
		//Computation of eigenvectors
		void F77NAME(dseupd) (const int& rvec,   const char* howmny,
							  const int* select, double* d,
							  double* z ,         const int& ldz,
							  const double& sigma,const char* bmat,
							  const int& n,       const char* which,
							  const int& nev,     const double& tol,
							  double* resid,      const int& ncv,
							  double* v,          const int& ldv,
							  const int* iparam,  int* ipntr,         
							  const double* workd,double* workl,      
							  const int& lworkl,  int& info);
		
		void F77NAME(dneupd) (const int& rvec,    const char* howmny,
							  const int* select,  double* dr,
							  double* di,          double* z ,         
							  const int& ldz,      const double& sigmar,      
							  const double& sigmai,double* workev,
							  const char* bmat,    const int& n,
							  const char* which,   const int& nev,      
							  const double& tol,   double* resid,
							  const int& ncv,      double* v,
							  const int& ldv,      int* iparam,   
							  int* ipntr,          double* workd,       
							  double* workl,       const int& lworkl,   int& info);
	
	}

//#ifdef NEKTAR_USING_ARPACK


    ///Top level reverse communication interface to solve real double-precision symmetric problems
	
	static inline void Dsaupd ( int& ido,           const char* bmat,
							    const int& n,       const char* which,
							    const int& nev,     const double& tol,
								double* resid,      const int& ncv,
	                            double* v,          const int& ldv,
							    int* iparam,  int* ipntr,         
							    double* workd,      double* workl,      
							    const int& lworkl,  int& info)
    {
        F77NAME(dsaupd) (ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
    }
	
	

	/// Post-processing routine to computed eigenvector of computed eigenvalues in Dsaupd
	static inline void Dseupd (const int& rvec,   const char* howmny,
						  const int* select, double* d,
						  double* z ,         const int& ldz,
						  const double& sigma,      const char* bmat,
						  const int& n,       const char* which,
						  const int& nev,     const double& tol,
						  double* resid,      const int& ncv,
						  double* v,          const int& ldv,
						  int* iparam,        int* ipntr,         
						  double* workd,      double* workl,      
						  const int& lworkl,  int& info)
	{
        F77NAME(dseupd) (rvec,howmny, select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
    }
	
	
	///Top level reverse communication interface to solve real double-precision non-symmetric problems
	static inline void Dnaupd (    int& ido,           const char* bmat,
								   const int& n,       const char* which,
								   const int& nev,     const double& tol,
								   double* resid,      const int& ncv,
								   double* v,          const int& ldv,
								   int* iparam,        int* ipntr,         
								   double* workd,      double* workl,      
								   const int& lworkl,  int& info)
		{
			F77NAME(dnaupd) (ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
		}

		/// Post-processing routine to computed eigenvector of computed eigenvalues in Dnaupd
	
 static inline void Dneupd(const int& rvec,    const char* howmny,
						   const int* select,  double* dr,
						   double* di,          double* z ,         
                           const int& ldz,      const double& sigmar,      
                           const double& sigmai,double* workev,
                           const char* bmat,    const int& n,
                           const char* which,   const int& nev,      
                           const double& tol,   double* resid,       
						   const int& ncv,      double* v,
						   const int& ldv,      int* iparam,   
						   int* ipntr,          double* workd,       
                           double* workl,       const int& lworkl,   int& info)
		{
			F77NAME(dneupd) (rvec,howmny,select,dr,di,z,ldz,sigmar,sigmai,workev,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
		}
		
		
			
//#endif //NEKTAR_USING_ARPACK
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_ARPACK_HPP

