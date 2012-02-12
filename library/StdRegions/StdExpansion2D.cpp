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
#include <SpatialDomains/Geometry2D.h>

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

    } //end namespace
} //end namespace


/**
* $Log: StdExpansion2D.cpp,v $
* Revision 1.30  2008/08/28 15:03:54  pvos
* small efficiency updates
*
* Revision 1.29  2008/08/20 09:14:57  sherwin
* In TensorDeriv replaced comparison of arrays with comparison of stored pointer
*
* Revision 1.28  2008/08/14 22:09:50  sherwin
* Modifications to remove HDG routines from StdRegions and removed StdExpMap
*
* Revision 1.27  2008/07/05 16:15:44  sherwin
* Corrected issue when not using Nektar_using_Blac for physderiv if intput and output are the same
*
* Revision 1.26  2008/07/04 10:18:40  pvos
* Some updates
*
* Revision 1.25  2008/06/05 15:06:06  pvos
* Added documentation
*
* Revision 1.24  2008/05/07 16:04:57  pvos
* Mapping + Manager updates
*
* Revision 1.23  2008/04/22 05:22:15  bnelson
* Speed enhancements.
*
* Revision 1.22  2008/04/06 06:04:15  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.21  2008/04/03 16:12:11  pvos
* updates for NEKTAR_USING_DIRECT_BLAS_CALLS
*
* Revision 1.20  2008/03/18 14:15:45  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.19  2008/03/12 15:25:09  pvos
* Clean up of the code
*
* Revision 1.18  2007/11/08 16:55:14  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.17  2007/10/15 20:37:14  ehan
* Make changes of column major matrix
*
* Revision 1.16  2007/09/27 12:55:57  pvos
* Column major Blas calls corrections
*
* Revision 1.15  2007/07/20 02:16:53  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.14  2007/05/30 23:56:55  sherwin
* Silly errors
*
* Revision 1.13  2007/05/22 02:01:41  bnelson
* Changed Array::size to Array::num_elements.
*
* Fixed some compiler errors in assertions.
*
* Revision 1.12  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.11  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.10  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.9  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.8  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.7  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.6  2007/03/05 19:28:50  bcarmo
* StdExpansion2D.cpp modified according to StdExpansion1D. Compiles.
*
* Revision 1.6  2007/03/05 22:35:21  bcarmo
* Version with StdExpansion2D compiling
*
* Revision 1.5  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.4  2007/01/18 18:44:45  bnelson
* Updates to compile on Visual Studio 2005.
*
* Revision 1.3  2007/01/17 16:05:35  pvos
* updated doxygen documentation
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.16  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.15  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.14  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.13  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.12  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/


