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

#ifndef H_EXPLIST
#define H_EXPLIST

#include <MultiRegions/MultiRegions.hpp>
#include <LibUtilities/NekMemoryManager.hpp>

namespace Nektar
{
  namespace MultiRegions
  {

    class ExpList
    {
    public:
      ExpList();
      ~ExpList();
      
      inline int GetNcoeffs(void)
      {
	return m_ncoeffs;
      }
      
      inline int GetPointsTot(void)
      {
	return m_npoints;
      }
      
      inline double *GetCoeffs(void)
      {
	return (m_coeffs);
      }
      
      inline double *GetPhys(void)
      {
	return (m_phys);
      }
      
      inline void SetTransState(TransState transState)
      {
	m_transState = transState;
      }
      
      inline void SetPhysState(bool physState)
      {
	m_physState = physState;
      }
  
      void  BwdTrans (double *outarray)
      {
	v_BwdTrans (outarray);
      }
    
    protected:
      int m_ncoeffs; 
      int m_npoints;
      
      double *m_coeffs;
      double *m_phys;
      
      TransState m_transState;
      bool       m_physState;
      
    private:
      //virtuals
      virtual void   v_BwdTrans (double *outarray)  = 0;
    };
    
  } //end of namespace
} //end of namespace
  
#endif // EXPLIST_H
