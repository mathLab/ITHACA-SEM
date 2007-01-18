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
	    ~ExpList2D();
	    
	    
	    ExpList2D(const StdRegions::BasisKey &TriBa, 
		      const StdRegions::BasisKey &TriBb, 
		      const StdRegions::BasisKey &QuadBa, 
		      const StdRegions::BasisKey &QuadBb, 
		      SpatialDomains::MeshGraph2D &graph2D,
		      SpatialDomains::Domain &domain = NULL,
		      const StdRegions::NodalBasisType 
		              TriNb = StdRegions::SIZE_NodalBasisType);
	    
	protected:
	    
	private:
	    
	    virtual void v_BwdTrans(double *outarray)
	    {
		BwdTrans(outarray);
	    }
	};
	
        typedef boost::shared_ptr<ExpList2D>      ExpList2DSharedPtr;
        typedef std::vector< ExpList2DSharedPtr > ExpList2DVector;
        typedef std::vector< ExpList2DSharedPtr >::iterator ExpList2DVectorIter;
    } //end of namespace
} //end of namespace

#endif
