///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.cpp
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

#include <MultiRegions/ExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ExpList1D::ExpList1D()
	{
	}
	
	ExpList1D::~ExpList1D()
	{
	}
	
        ExpList1D::ExpList1D(const ExpList1D &In):
            ExpList(In)
        {
            m_coeffs = MemoryManager::AllocateSharedArray<NekDouble>(m_ncoeffs);
            m_phys   = MemoryManager::AllocateSharedArray<NekDouble>(m_npoints);
        }

	ExpList1D::ExpList1D(const LibUtilities::BasisKey &Ba, 
			     const SpatialDomains::MeshGraph1D &graph1D)
	{
	    LocalRegions::SegExpSharedPtr seg;
	    SpatialDomains::SegGeomVector SegGeoms = graph1D.GetSeggeoms();
	    SpatialDomains::SegGeomVectorIter def;
	    
	    m_ncoeffs = SegGeoms.size()*Ba.GetNumModes();
	    m_npoints = SegGeoms.size()*Ba.GetNumPoints();
	    
	    m_transState = eNotSet; 
	    m_physState  = false;
	    
	    for(def = SegGeoms.begin(); def != SegGeoms.end(); ++def)
	    {
		seg = MemoryManager::AllocateSharedPtr<LocalRegions::SegExp>(Ba, *def);
		(*m_exp).push_back(seg);
	    }

            m_coeffs = MemoryManager::AllocateSharedArray<NekDouble>(m_ncoeffs);
            m_phys   = MemoryManager::AllocateSharedArray<NekDouble>(m_npoints);
	}
    } //end of namespace
} //end of namespace
