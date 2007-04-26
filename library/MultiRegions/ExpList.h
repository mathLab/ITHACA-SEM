///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.h
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
// Description: Expansion list top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H
#define NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H

#include <MultiRegions/MultiRegions.hpp>
#include <StdRegions/StdExpansion.h>

namespace Nektar
{
  namespace MultiRegions
  {

    class ExpList
    {
    public:
      ExpList();
      virtual ~ExpList();
      
      inline int GetNcoeffs(void) const
      {
	return m_ncoeffs;
      }
      
      inline int GetPointsTot(void) const
      {
	return m_npoints;
      }
      
      inline void SetTransState(const TransState transState)
      {
	m_transState = transState;
      }
      
      inline void SetPhysState(const bool physState)
      {
	m_physState = physState;
      }

      inline bool GetPhysState(void) const
      {
	return m_physState;
      }
      
      NekDouble Integral(ConstNekDoubleSharedArray inarray);
      void   IProductWRTBase(ConstNekDoubleSharedArray inarray, 
			     NekDoubleSharedArray &outarray);
      void   IProductWRTBase(const ExpList &S1);
      void   IProductWRTBase(const ExpList &S1, 
			     NekDoubleSharedArray &outarray);
      void   PhysDeriv(ConstNekDoubleSharedArray inarray,
		       NekDoubleSharedArray &out_d0, 
		       NekDoubleSharedArray &out_d1 = NullNekDoubleSharedArray,
		       NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray);
      void   FwdTrans (ConstNekDoubleSharedArray inarray,
		       NekDoubleSharedArray &outarray);
      void   FwdTrans (ConstNekDoubleSharedArray inarray, 
                       ExpList &Sout);
      void   FwdTrans (const ExpList &Sin);
      void   BwdTrans (ConstNekDoubleSharedArray inarray, 
		       NekDoubleSharedArray &outarray); 
      void   BwdTrans (const ExpList &Sin); 
      
      void   GetCoords(NekDoubleSharedArray &coord_0,
		       NekDoubleSharedArray &coord_1 = NullNekDoubleSharedArray,
		       NekDoubleSharedArray &coord_2 = NullNekDoubleSharedArray);
      void   GetCoords(NekDoubleArrayVector &coords);

      void   WriteToFile(std::ofstream &out);
    
      inline int GetCoordim(int eid)
      {
	  ASSERTL2(eid <= m_exp.size(),"eid is larger than number of elements");
	
	  return m_exp[eid]->GetCoordim();
      }
      

      inline ConstNekDoubleSharedArray GetCoeffs() const
      {
          return m_coeffs;
      }

      inline NekDoubleSharedArray &UpdateCoeffs()
      {
          return m_coeffs;
      }

      NekDouble Linf (ConstNekDoubleSharedArray sol);
      NekDouble L2   (ConstNekDoubleSharedArray sol);


      inline int GetExpSize(void)
      {
	  return m_exp.size();
      }
      

      inline StdRegions::StdExpansionSharedPtr& GetExp(int n)
      {
	  return m_exp[n];
      }
      
    protected:
      int m_ncoeffs; 
      int m_npoints;
      NekDoubleSharedArray m_coeffs;

      TransState m_transState;
      bool       m_physState;
     
      StdRegions::StdExpansionVector m_exp;
      
    private:
      
    };

    
  } //end of namespace
} //end of namespace
  
#endif // EXPLIST_H
