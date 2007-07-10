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
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

	class ExpList1D: 
	    public ExpList
	{

	public:
	    
	    ExpList1D();

	    ExpList1D(const ExpList1D &In);
	    
            //	    ExpList1D(const LibUtilities::BasisKey &Ba, 
	    //      const SpatialDomains::MeshGraph1D &graph1D);

            ExpList1D(const LibUtilities::BasisKey &Ba, 
                      const SpatialDomains::Composite &cmps
);
	    
	    ~ExpList1D();


            void   PhysDeriv  (ExpList &S0,
                               ExpList &S1 = NullExpList, 
                               ExpList &S2 = NullExpList)
            {
                ExpList::PhysDeriv(S0,S1,S2);
            }

	protected:
	    
	private:
            
	};
	
        typedef boost::shared_ptr<ExpList1D>      ExpList1DSharedPtr;
        typedef std::vector<ExpList1DSharedPtr>   ExpList1DVector;
        typedef std::vector<ExpList1DSharedPtr>::iterator ExpList1DVectorIter;
	
    } //end of namespace
} //end of namespace

#endif//EXPLIST1D_H

/**
* $Log: ExpList1D.h,v $
* Revision 1.12  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.11  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
