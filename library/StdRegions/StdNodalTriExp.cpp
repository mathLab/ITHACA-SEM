///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.cpp
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
// Description: Nodal triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/NekMemoryManager.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        StdMatrix StdNodalTriExp::s_elmtmats;

        StdNodalTriExp::StdNodalTriExp(const BasisKey &Ba, const BasisKey &Bb, 
            NodalBasisType Ntype):
        StdTriExp(Ba,Bb)
        {
            m_nbtype = Ntype;

            ASSERTL0(m_base[0]->GetBasisOrder() == m_base[1]->GetBasisOrder(),
                "Nodal basis initiated with different orders in the a and b directions");
        }

        // Destructor
        StdNodalTriExp::~StdNodalTriExp()
        { 
        }

        void StdNodalTriExp::GenNBasisTransMatrix(double * outarray)
        {
            int             i,j;
            int             tot_order = GetNcoeffs();
            int             nquad;
            const double*   r; 
            const double*   s; 
            const double*   t;
            double          c[2];

            NBasisManagerSingleton::Instance().GetNodePoints(m_nbtype,
                m_base[0]->GetBasisOrder(),r,s,t);

            for(i = 0; i < tot_order; ++i)
            {
                // fill physical space with mode i
                StdTriExp::FillMode(i,m_phys);

                // interpolate mode i to the Nodal points 'j' and store in outarray
                for(j = 0; j < tot_order; ++j)
                {
                    c[0] = r[j];
                    c[1] = s[j];
                    // define matrix in row major format to have rows of 
                    // all the different expansion bases defined at the nodal point 
                    outarray[j*tot_order+i] = StdTriExp::Evaluate(c);
                }
            }
        }

        void StdNodalTriExp::NodalToModal()
        {
            StdMatContainer *M;

            M = GetNBasisTransMatrix();
            M->Solve(m_coeffs,1);
        }

        StdMatContainer * StdNodalTriExp::GetNBasisTransMatrix() 
        {
            StdMatContainer * mat;
            mat = s_elmtmats.GetNBasisTrans(this);
            return mat;
        }

        void StdNodalTriExp::ModalToNodal()
        {
            StdMatContainer *M;

            M = GetNBasisTransMatrix();
            M->Mxv(m_coeffs,m_coeffs);
        }


        //////////////////////////////
        /// Integration Methods
        //////////////////////////////


        void StdNodalTriExp::IProductWRTBase(const double * inarray, 
            double * outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                inarray,outarray);
        }

        /** \brief Calculate the inner product of inarray with respect to
        the basis B=base0[p]*base1[pq] and put into outarray:

        This function uses the StdTriExp routine and then 
        calls ModalToNodal to transform to Nodal basis
        **/

        void StdNodalTriExp:: IProductWRTBase(const double *base0, 
            const double *base1, 
            const double *inarray,
            double *outarray)
        {
            // Take inner product with respect to Orthgonal basis using
            // StdTri routine
            StdTriExp::IProductWRTBase(base0,base1,inarray,outarray);

            // Transform to Nodal Expansion
            ModalToNodal();
        }


        /** Fill outarray with nodal mode 'mode' of expansion and put in m_phys
        */

        void StdNodalTriExp::FillMode(const int mode, double *outarray)
        {

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs,m_coeffs,1);

            m_coeffs[mode] = 1.0;
            NodalToModal();

            StdTriExp::BwdTrans(outarray);
        }

        void StdNodalTriExp::GenMassMatrix(double * outarray)
        {
            int      i;
            BstShrDArray tmp = GetDoubleTmpSpace(GetPointsOrder(0)*
                GetPointsOrder(1));

            for(i = 0; i < m_ncoeffs; ++i)
            {
                FillMode(i,tmp.get());
                IProductWRTBase(tmp.get(),outarray+i*m_ncoeffs);
            }
        }

        void StdNodalTriExp::GenLapMatrix(double * outarray)
        {
            ASSERTL0(false, "Not implemented");
        }

        StdMatContainer * StdNodalTriExp::GetMassMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalMass(this);
            return tmp;
        }

        StdMatContainer * StdNodalTriExp::GetLapMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalLap(this);
            return tmp;
        }

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        void StdNodalTriExp::Deriv(const double *inarray, double *outarray_d0, 
            double *outarray_d1)
        {
            StdTriExp::Deriv(inarray, outarray_d0, outarray_d1);
        }

        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////

        // Currently convert nodal values into tranformed values and
        // backward transform

        void StdNodalTriExp::BwdTrans(double * outarray)
        {
            double *tmp = new double[m_ncoeffs];

            // save nodal values
            Blas::Dcopy(m_ncoeffs,m_coeffs,1,tmp,1);
            NodalToModal();
            StdTriExp::BwdTrans(outarray);
            Blas::Dcopy(m_ncoeffs,tmp,1,m_coeffs,1);
        }

        void StdNodalTriExp::FwdTrans(const double * inarray)
        {
            StdTriExp::FwdTrans(inarray);
            ModalToNodal();
        }

        double StdNodalTriExp::Evaluate(const double * coords)
        {
            return StdTriExp::Evaluate(coords);
        }
    } // end of namespace
} // end of namespace

/** 
* $Log: StdNodalTriExp.cpp,v $
* Revision 1.19  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.18  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.17  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.16  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
**/ 

