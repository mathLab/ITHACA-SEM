///////////////////////////////////////////////////////////////////////////////
//
// File BasisManager.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Basis definition and manager 
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/BasisManager.h>
#include <StdRegions/StdBasis.h>

namespace Nektar
{
namespace StdRegions
{

  const Basis* BasisManager::GetBasis (const BasisKey &bkey)
  {

    std::vector<Basis*>::iterator def;
    int  id_btype = (int)bkey.GetBasisType();
    
    def = find(m_basis[id_btype].begin(),
	       m_basis[id_btype].end(),bkey);
    
    if(def != m_basis[id_btype].end())
    {

      // check that this definition has sufficiently high order
      if(bkey.GetBasisOrder() > def[0]->GetBasisOrder())
      {
	// reset bdata in this defintiion.
	def[0]->ResetBasisOrder(bkey.GetBasisOrder());
      }
      
      return (def[0]); 
    }
    else
    { // set up new defintiion
      int id = m_basis[id_btype].size();
      
      Basis * B = new Basis(bkey);
      
      m_basis[id_btype].push_back(B);
      
      return (B);
    }
  }

  const Basis* BasisManager::GetBasis (const BasisType btype, const int order,
				       const PointsType ptype, const int nq,
				       const double alpha, const double beta)
  {
    
    BasisKey local(btype,order,ptype,nq,alpha,beta);
    return(BasisManager::GetBasis(local));
  }
  
  const double* BasisManager::GetBasisArray (const BasisKey & bkey)
  {

    const Basis * B = GetBasis(bkey);
    return(B->GetBdata());

  }

  const double* BasisManager::GetBasisArray (const BasisType btype, 
			      const int order, const PointsType ptype, 
			      const int nq, const double alpha, 
			      const double beta)
  {
    const Basis * B = GetBasis(btype,order,ptype,nq,alpha,beta);
    return(B->GetBdata());
  }

}//end namespace
}//end namespace


/** 
 * $Log: BasisManager.cpp,v $
 * Revision 1.20  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 * Revision 1.19  2006/02/26 23:37:28  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 * Revision 1.18  2006/02/19 13:26:12  sherwin
 *
 * Coding standard revisions so that libraries compile
 *
 * Revision 1.17  2006/02/15 08:06:35  sherwin
 *
 * Put files into coding standard (although they do not compile)
 *
 *
 **/ 
