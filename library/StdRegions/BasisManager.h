///////////////////////////////////////////////////////////////////////////////
//
// File BasisManager.h
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
// Description: Header file of Basis manager
//
///////////////////////////////////////////////////////////////////////////////

#ifndef BASISMANAGER_H
#define BASISMANAGER_H

#include <StdRegions/PolyManager.h>
#include <StdRegions/StdBasis.h>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
namespace StdRegions
{

  class BasisManager: public PolyManager
  {
  public:
    BasisManager()
    {
      m_basis = new std::vector<Basis*>[SIZE_BasisType];
    }

    ~BasisManager()
    {
      std::vector<Basis*>::iterator def;

      for(int i=0;i<SIZE_BasisType;i++){
    for(def = m_basis[i].begin(); def != m_basis[i].end(); ++def)
      delete def[0];
      }

      delete[] m_basis;
    }

    const Basis* GetBasis (const BasisKey & bkey);

    const Basis* GetBasis (const BasisType btype, const int order,
               const PointsType ptype, const int nq,
               const double alpha, const double beta);

    const double* GetBasisArray (const BasisKey & bkey);

    const double* GetBasisArray (const BasisType btype, const int order,
                 const PointsType ptype, const int nq,
                 const double alpha, const double beta);

  protected:

  private:
    std::vector<Basis*> *m_basis;

  };

} //end of namespace
} //end of namespace

#endif //BASISMANAGER_H


/**
 * $Log: BasisManager.h,v $
 * Revision 1.19  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.18  2006/02/26 23:37:28  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 * Revision 1.17  2006/02/19 13:26:12  sherwin
 *
 * Coding standard revisions so that libraries compile
 *
 * Revision 1.16  2006/02/15 08:06:35  sherwin
 *
 * Put files into coding standard (although they do not compile)
 *
 *
 **/
