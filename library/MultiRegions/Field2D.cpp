///////////////////////////////////////////////////////////////////////////////
//
// File Field2D.cpp
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
// Description: Field definition for 2D domains
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/Field2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	Field2D::Field2D(void)
	{
	}

	Field2D::Field2D(const StdRegions::BasisKey &TriBa, 
			 const StdRegions::BasisKey &TriBb, 
			 const StdRegions::NodalBasisType  TriNb,
			 const StdRegions::BasisKey &QuadBa, 
			 const StdRegions::BasisKey &QuadBb, 
			 SpatialDomains::Domain &Domain2D)
	{
	    int i,nbnd;
	    ExpList1DSharedPtr  bndry;
	    BoundaryVectorIter def;
	    BoundaryVector Bndry = Domain2D.getBoundaries();

	    m_field.reset(new ContExpList2D(TriBa,TriBb,TriNb,QuadBa,QuadBa,
					    Domain2D.GetGraph2D()));

	    nbnd = Domain2D.GetBoundaries().size();

	    for(def = Bndry.begin(); def < Bndry.end(); ++def)
	    {
		bndry.reset(new ExpList1D(Ba,graph1D));
		m_bndConstraint.push_back(Bndry);
		m_bndType(Def of Bndry Type);/////////////////////////
	    } 

 -- Get element and fac id's from graph - ask joe

Set up numbering to local and global definition. 

	}


	Field2D::~Field2D()
	{
	}

    } // end of namespace
} //end of namespace
