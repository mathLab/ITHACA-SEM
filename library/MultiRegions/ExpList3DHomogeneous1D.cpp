///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3DHomogeneous1D.cpp
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
// Description: A 2D field which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList3Dhomogeneous1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs

        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D():
            ExpList(),
            m_nzplanes(0)
        {
        }

        // Constructor for ExpList3DHomogeneous1D to act as a Explist2D field
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const int nzplanes, 
                                     const NekDouble lz, 
                                     SpatialDomains::MeshGraph2D &graph2D):
            ExpList(),
            m_nzplanes(nzplanes),
            m_planes(nzplanes),
            m_lz(lz)
        {
            int n,j,nel;
            bool value = false; 

            m_planes[0] = MemoryManager<ExpList2D>::AllocateSharedPtr(graph2D,
                                                                      value); 
            m_globalOptParam = m_planes[0]->GetGlobalOptParam();

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_planes[0]->GetExpSize();

            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(j));
            }

            for(n = 1; n < m_nzplanes; ++n)
            {
                m_planes[n] = MemoryManager<ExpList2D>::AllocateSharedPtr(*m_planes[0],value);
                for(j = 0; j < nel; ++j)
                {
                    (*m_exp).push_back((*m_exp)[j]);
                }
            }            

            SetCoeffPhys();

        }


        /**
         * @param   In          ExpList3DHomogeneous1D object to copy.
         */
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const ExpList3DHomogeneous1D &In):
            ExpList(In,false),
            m_nzplanes(In.m_nzplanes),
            m_planes(In.m_planes),    // soft copy
            m_lz(m_lz)
        {
            SetCoeffPhys();
        }

        /**
         * Destructor
         */
        ExpList3DHomogeneous1D::~ExpList3DHomogeneous1D()
        {
        }

        void ExpList3DHomogeneous1D::SetCoeffPhys(void)
        {
            int i,n,cnt;
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            // Set total coefficients and points
            m_ncoeffs = ncoeffs_per_plane*m_nzplanes;
            m_npoints = npoints_per_plane*m_nzplanes;
            
            m_coeffs = Array<OneD, NekDouble> (m_ncoeffs);
            m_phys   = Array<OneD, NekDouble> (m_npoints);

            int nel = m_planes[0]->GetExpSize();
            m_coeff_offset   = Array<OneD,int>(nel*m_nzplanes);
            m_phys_offset    = Array<OneD,int>(nel*m_nzplanes);
            m_offset_elmt_id = Array<OneD,int>(nel*m_nzplanes);
            Array<OneD, NekDouble> tmparray;

            for(cnt  = n = 0; n < m_nzplanes; ++n)
            {
                m_planes[0]->SetCoeffsArray(tmparray= m_coeffs + ncoeffs_per_plane*n);
                m_planes[0]->SetPhysArray(tmparray = m_phys + npoints_per_plane*n);

                for(i = 0; i < nel; ++i)
                {
                    m_coeff_offset[cnt] = m_planes[n]->GetCoeff_Offset(i) + n*ncoeffs_per_plane;
                    m_phys_offset[cnt] =  m_planes[n]->GetPhys_Offset(i) + n*npoints_per_plane;
                    m_offset_elmt_id[cnt++] = m_planes[n]->GetOffset_Elmt_Id(i) + n*nel;
                }
            }
        }
        

        /**
         * The operation calls the 2D plane coordinates through the
         * function ExpList#GetCoords and then evaluated the third
         * coordinate using the member \a m_lz
         *
         * @param coord_0 After calculation, the \f$x_1\f$ coordinate
         *                          will be stored in this array.
         *
         * @param coord_1 After calculation, the \f$x_2\f$ coordinate
         *                          will be stored in this array.
         *
         * @param coord_2 After calculation, the \f$x_3\f$ coordinate
         *                          will be stored in this array. This
         *                          coordinate is evaluated using the
         *                          predefined value \a m_lz
         */
        void ExpList3DHomogeneous1D::v_GetCoords(Array<OneD, NekDouble> &xc0,
                                                 Array<OneD, NekDouble> &xc1,
                                                 Array<OneD, NekDouble> &xc2)
        {
            int n;
            Array<OneD, NekDouble> tmp_xc_0;
            NekDouble val,dz;
            int npoints = m_planes[0]->GetTotPoints();

            ExpList::v_GetCoords(xc0,xc1);

            dz = m_lz/m_nzplanes;
            
            for(n = 0; n < m_nzplanes; ++n)
            {
                val = n*dz; 
                Vmath::Fill(npoints,val,tmp_xc_0 = xc2 + npoints*n,1);
            }
        }
        
        
        // Possibly could just leave this to the default method which
        // would set up a block array over all planes. 
        void ExpList3DHomogeneous1D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_nzplanes; ++n)
            {
                m_planes[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1, UseContCoeffs);
                cnt   += m_planes[n]->GetTotPoints();
                cnt1  += m_planes[n]->GetNcoeffs();
            }
        }

    } //end of namespace
} //end of namespace


/**
* $Log: v $
*
**/

