///////////////////////////////////////////////////////////////////////////////
//
// File ContExpList1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Continusou Expansion list definition in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CONTEXPLIST1D_H
#define CONTEXPLIST1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>

#include <StdRegions/StdMatrix.h>

namespace Nektar
{
  namespace MultiRegions
  {
    
    class ContExpList1D: 
      public ExpList1D 
    {
    public:
      ContExpList1D();
      ContExpList1D(const StdRegions::BasisKey &Ba, 
		    SpatialDomains::MeshGraph1D &graph1D);
      ~ContExpList1D();
      
      inline int getContNcoeffs()
      {
	return m_contNcoeffs;
      }

      inline double *get_cont_coeffs()
      {
	return m_contCoeffs;
      }
    
      inline void ContToLocal()
      {
	ContToLocal(m_contCoeffs,m_coeffs);
      }
    
      inline void ContToLocal(const double *cont,double *loc)
      {
	Vmath::Gathr(m_ncoeffs,cont,m_locToContMap,loc);
      }

      inline void LocalToCont()
      {
	LocalToCont(m_coeffs,m_contCoeffs);
      }

      inline void LocalToCont(const double *loc, double *cont)
      {
	Vmath::Scatr(m_ncoeffs,loc,m_locToContMap,cont);
      }

      inline void Assemble()
      {
	Assemble(m_coeffs,m_contCoeffs);
      }

      inline void Assemble(const double *loc, double *cont)
      {
	Vmath::Zero(m_contNcoeffs,cont,1);
	Vmath::Assmb(m_ncoeffs,loc,m_locToContMap,cont);
      }
		       
      void IProductWRTBase(const double *inarray, double *outarray);
      
      void FwdTrans(const double *inarray);

      void BwdTrans(double *outarray);
      
      void GenMassMatrix(void);

    protected:
    

    private:
      int    m_contNcoeffs;
      int    *m_locToContMap;
      
      double *m_contCoeffs;
      
      StdRegions::StdMatContainer *m_mass;
      
      virtual void v_BwdTrans(double *outarray)
      {
	BwdTrans(outarray);
      }
      
      
    };
  } //end of namespace
} //end of namespace

#endif // end of define
