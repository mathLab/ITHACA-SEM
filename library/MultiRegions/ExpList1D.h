///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.h
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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST1D_H
#define EXPLIST1D_H

#include <vector>
#include <fstream>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <StdRegions/StdBasis.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph1D.h>

namespace Nektar
{
  namespace MultiRegions
  {

    typedef boost::shared_ptr<SegExp> SharedSegExpPtr;
    typedef std::vector< SharedSegExpPtr > SegExpVector;
    typedef std::vector< SharedSegExpPtr >::iterator SegExpVectorIter;


    class ExpList1D: 
      public ExpList
    {

    public:

      ExpList1D();

      ExpList1D(const StdRegions::BasisKey &Ba, 
		SpatialDomains::MeshGraph1D &graph1D);

      ~ExpList1D();

      double Integral(const double *inarray);
      void   IProductWRTBase(const double *inarray, double *outarray);
      void   IProductWRTBase(ExpList1D &S1, ExpList1D &S2);
      void   IProductWRTBase(ExpList1D &S1, double * outarray);
      void   Deriv    (const int n, double **outarray);
      void   Deriv    (const int n, const double *inarray, double ** outarray);
      void   FwdTrans (const double *inarray);
      void   BwdTrans (double *outarray); 
      
      void   GetCoords(double **coords);
      void   WriteToFile(ofstream &out);
    
      inline int GetCoordim(int eid)
      {
	ASSERTL2(eid <= m_seg.size(),"ExpList1D:get_coordim()",
		 "eid is larger than number of elements");
	return m_seg[eid]->GetCoordim();
      }
      
      double Linf (const double *sol);
      double L2   (const double *sol);
      
    protected:
      SegExpVector  m_seg;

    private:
      
      virtual void v_BwdTrans(double *outarray)
      {
	BwdTrans(outarray);
      }
      
    };
    
  } //end of namespace
} //end of namespace

#endif//EXPLIST1D_H
