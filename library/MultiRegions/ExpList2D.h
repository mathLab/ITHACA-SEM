///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2D.h
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
// Description: Expansion list 2D header definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST2D_H
#define EXPLIST2D_H

#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <StdRegions/StdBasis.h>

#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <SpatialDomains/MeshGraph2D.h>

namespace Nektar
{
  namespace MultiRegions
  {


    class ExpList2D:
      public ExpList
    {
    public:
      ExpList2D();
      ~ExpList2D();

      ExpList2D(const StdRegions::BasisKey &TriBa, 
		const StdRegions::BasisKey &TriBb, 
		const StdRegions::BasisKey &QuadBa, 
		const StdRegions::BasisKey &QuadBb, 
		SpatialDomains::MeshGraph2D &graph2D);
      

      double Integral(const double *inarray);
      void   IProductWRTBase(const double *inarray, double *outarray);
      void   IProductWRTBase(ExpList2D &S1, ExpList2D &S2);
      void   IProductWRTBase(ExpList2D &S1, double * outarray);
      void   Deriv      (const int n, double **outarray);
      void   Deriv      (const int n, const double *inarray, double ** outarray);
      void   FwdTrans   (const double *inarray);
      void   BwdTrans   (double *outarray); 
      
      void   GetCoords  (double **coords);
      void   WriteToFile(std::ofstream &out);

      inline int GetCoordim(int eid)
      {
	  if(m_tri.size())
	  {
	      ASSERTL2(eid <= m_tri.size(),
		       "eid is larger than number of elements");
	      return m_tri[eid]->GetCoordim();
	  }
	  else
	  {
	      ASSERTL2(eid <= m_quad.size(),
		       "eid is larger than number of elements");
	      return m_quad[eid]->GetCoordim();
	  }
      }

      double Linf (const double *sol);
      double L2   (const double *sol);

    protected:
      LocalRegions::TriExpVector  m_tri;
      LocalRegions::QuadExpVector m_quad;

    private:
    
      virtual void v_BwdTrans(double *outarray)
      {
	BwdTrans(outarray);
      }
    };

  } //end of namespace
} //end of namespace

#endif
