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
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/NodalTriExp.h>
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
            
            ExpList2D(const ExpList2D &In);   
            
            ExpList2D(const LibUtilities::BasisKey &TriBa, 
                      const LibUtilities::BasisKey &TriBb, 
                      const LibUtilities::BasisKey &QuadBa, 
                      const LibUtilities::BasisKey &QuadBb, 
                      const SpatialDomains::MeshGraph2D &graph2D,
                      const LibUtilities::PointsType 
                      TriNb = LibUtilities::SIZE_PointsType);

            ExpList2D(SpatialDomains::MeshGraph2D &graph2D);
            
            ~ExpList2D();
            
            void   PhysDeriv  (ExpList &S0,
                               ExpList &S1, 
                               ExpList &S2 = NullExpList)
            {
                ExpList::PhysDeriv(S0,S1,S2);
            }
            
        protected:
            
        private:
            
        };
        
        typedef boost::shared_ptr<ExpList2D>      ExpList2DSharedPtr;
        typedef std::vector< ExpList2DSharedPtr > ExpList2DVector;
        typedef std::vector< ExpList2DSharedPtr >::iterator ExpList2DVectorIter;
    } //end of namespace
} //end of namespace

#endif//EXPLIST2D_H

/**
* $Log: ExpList2D.h,v $
* Revision 1.11  2007/07/22 23:04:21  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.10  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

