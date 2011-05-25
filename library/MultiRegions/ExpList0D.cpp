///////////////////////////////////////////////////////////////////////////////
//
// File ExpList0D.cpp
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
// Description: Expansion list 0D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList0D.h>
#include <LibUtilities/Polylib/Polylib.h>

namespace Nektar
{
    namespace MultiRegions
    {

		/**
         * Default constructor ExpList0D object.
         */
        ExpList0D::ExpList0D():
		ExpList()
        {
        }
		
        /**
         * Creates an identical copy of another ExpList0D object.
         */
        ExpList0D::ExpList0D(const ExpList0D &In, bool DeclareCoeffPhysArrays):
		ExpList(In,DeclareCoeffPhysArrays)
        {
        }


        /**

         */
        ExpList0D::ExpList0D(const SpatialDomains::VertexComponentSharedPtr &m_geom):
		ExpList()
        {
			m_point = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(m_geom);
			
			m_ncoeffs = 1;
			m_npoints = 1;
			
			// Set up m_coeffs, m_phys.
			m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
			m_phys   = Array<OneD, NekDouble>(m_npoints);
			
			m_coeffs[0] = m_point->GetCoeff(0);
			m_phys[0]  = m_point->GetPhys(0);
			
			m_coeffs = m_point->UpdateCoeffs();
			m_phys   = m_point->UpdatePhys();
        }

        /**
         *
         */
        ExpList0D::~ExpList0D()
        {
        }
		
		void ExpList0D::v_GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
        {
			m_point->GetCoords(x,y,z);
		}
		
		void ExpList0D::v_GetCoord(Array<OneD,NekDouble> &coords)
        {
			m_point->GetCoords(coords);
		}
		
		void ExpList0D::v_SetCoeff(NekDouble val)
        {
			m_point->SetCoeff(val);
		}
		
		void ExpList0D::v_SetPhys(NekDouble val)
        {
			m_point->SetPhys(val);
		}
		
		const SpatialDomains::VertexComponentSharedPtr &ExpList0D::v_GetGeom(void) const
		{
			return m_point->GetGeom();
		}
		
		const SpatialDomains::VertexComponentSharedPtr &ExpList0D::v_GetVertex(void) const
		{
			return m_point->GetVertex();
		}
		
    } //end of namespace
} //end of namespace

/**
 * $Log: v $
 *
 **/