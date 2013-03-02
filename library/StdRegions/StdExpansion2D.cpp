///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion2D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion2D::StdExpansion2D()
        {
        }

        StdExpansion2D::StdExpansion2D(int numcoeffs,
                       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb):
        StdExpansion(numcoeffs,2, Ba, Bb)
        {
        }

        StdExpansion2D::StdExpansion2D(const StdExpansion2D &T):
                StdExpansion(T)
        {
        }

        StdExpansion2D::~StdExpansion2D()
        {
        }

        //----------------------------
        // Differentiation Methods
        //----------------------------
        void StdExpansion2D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, NekDouble> &outarray_d0,
                         Array<OneD, NekDouble> &outarray_d1)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            if (outarray_d0.num_elements() > 0) // calculate du/dx_0
            {
                DNekMatSharedPtr D0 = m_base[0]->GetD();
                if(inarray.data() == outarray_d0.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &wsp[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
                else
                {
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &inarray[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
            }

            if (outarray_d1.num_elements() > 0) // calculate du/dx_1
            {
                DNekMatSharedPtr D1 = m_base[1]->GetD();
                if(inarray.data() == outarray_d1.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &wsp[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
                else
                {
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
            }
        }

        NekDouble StdExpansion2D::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
        {
            return PhysEvaluate(coords,m_phys);
        }


        NekDouble StdExpansion2D::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble> & physvals)
        {
            NekDouble val;
            int i;
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp1(nq1);

            DNekMatSharedPtr I = m_base[0]->GetI(coords);;

            ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
            ASSERTL2(coords[0] <  1 + NekConstants::kNekZeroTol, "coord[0] >  1");
            ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
            ASSERTL2(coords[1] <  1 + NekConstants::kNekZeroTol, "coord[1] >  1");

            // interpolate first coordinate direction
            for (i = 0; i < nq1;++i)
            {
                wsp1[i] = Blas::Ddot(nq0, &(I->GetPtr())[0], 1,
                                     &physvals[i * nq0], 1);
            }

            // interpolate in second coordinate direction
            I = m_base[1]->GetI(coords+1);
            val = Blas::Ddot(nq1, I->GetPtr(), 1, wsp1, 1);

            return val;
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        NekDouble StdExpansion2D::Integral(const Array<OneD, const NekDouble>& inarray,
                       const Array<OneD, const NekDouble>& w0,
                       const Array<OneD, const NekDouble>& w1)
        {
            int i;
            NekDouble Int = 0.0;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad0 * nquad1);

            // multiply by integration constants
            for (i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, &inarray[0] + i*nquad0, 1, w0.get(),
                            1, &tmp[0] + i*nquad0, 1);
            }

            for (i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1, &tmp[0]+ i, nquad0, w1.get(), 1,
                            &tmp[0] + i, nquad0);
            }
            Int = Vmath::Vsum(nquad0 * nquad1, tmp, 1);

            return Int;
        }

        void StdExpansion2D::v_SetUpPhysNormals(const int edge)
        {
           ComputeEdgeNormal(edge);
        }

        const NormalVector & StdExpansion2D::v_GetEdgeNormal(const int edge) const
        {
            std::map<int, NormalVector>::const_iterator x;
            x = m_edgeNormals.find(edge);
            ASSERTL0 (x != m_edgeNormals.end(),
                        "Edge normal not computed.");
            return x->second;
        }
        
        void StdExpansion2D::v_NegateEdgeNormal(const int edge)
        {
            m_negatedNormals[edge] = true;
            for (int i = 0; i < GetCoordim(); ++i)
            {
                Vmath::Neg(m_edgeNormals[edge][i].num_elements(), 
                           m_edgeNormals[edge][i], 1);
            }
        }

        bool StdExpansion2D::v_EdgeNormalNegated(const int edge)
        {
            return m_negatedNormals[edge];
        }
    } //end namespace
} //end namespace
