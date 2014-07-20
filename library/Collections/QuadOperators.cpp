///////////////////////////////////////////////////////////////////////////////
//
// File: QuadOperators.cpp
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
// Description: Operators specific to the quadrilateral
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Operator.h>

namespace Nektar 
{
    namespace Collections 
    {

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */       
        class PhysDeriv_SumFac_Quad : public Operator
        {
        public:
            PhysDeriv_SumFac_Quad(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom)
                : Operator (pExp, pGeom),
                  m_nquad0 (pExp->GetNumPoints(0)),
                  m_nquad1 (pExp->GetNumPoints(1)),
                  m_df     (2*pExp->GetCoordim(), m_nquad0*m_nquad1*m_numElmt)
            {
                const int coordim = pExp->GetCoordim();
                const int nqtot   = m_nquad0 * m_nquad1;
                int cnt = 0;
                
                for (int i = 0; i < pGeom.size(); ++i)
                {
                    const Array<TwoD, const NekDouble>& df =
                        pGeom[i]->GetGeomFactors()->GetDerivFactors(
                                                                    pExp->GetPointsKeys());
                    
                    for (int j = 0; j < 2*coordim; ++j)
                    {
                        if (pGeom[i]->GetGeomFactors()->GetGtype() ==
                            SpatialDomains::eDeformed)
                        {
                            Vmath::Vcopy(nqtot, &df[j][0], 1, &m_df[j][cnt], 1);
                        }
                        else
                        {
                            Vmath::Fill(nqtot, df[j][0], &m_df[j][cnt], 1);
                        }
                    }
                    
                    cnt += nqtot;
                }
                
                m_wspSize = 2 * m_df[0].num_elements();
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &wsp)
            {
                const int coordim = m_stdExp->GetCoordim();
                const int nqcol   = m_df[0].num_elements();
                const int nqtot   = m_nquad0 * m_nquad1;
                
                ASSERTL1(wsp.num_elements() == m_wspSize,
                         "Incorrect workspace size");
                ASSERTL1(input.num_elements() == nqcol,
                         "Incorrect input size");
                ASSERTL1(output.num_elements() == coordim*nqcol,
                         "Incorrect output size");
                
                Array<OneD, NekDouble> diff0(nqcol, wsp             );
                Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);
                Array<OneD, NekDouble> out0 (nqcol, output          );
                Array<OneD, NekDouble> out1 (nqcol, output +   nqcol);
                
                DNekMatSharedPtr D0 = m_stdExp->GetBasis(0)->GetD();
                DNekMatSharedPtr D1 = m_stdExp->GetBasis(1)->GetD();
                
                Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt, m_nquad0, 1.0,
                            D0->GetRawPtr(), m_nquad0, input.get(), m_nquad0, 0.0,
                            diff0.get(), m_nquad0);
                
                int cnt = 0;
                for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                                input.get() + cnt, m_nquad0, D1->GetRawPtr(), m_nquad1, 0.0,
                                diff1.get() + cnt, m_nquad0);
                }
                
                Vmath::Vmul  (nqcol, m_df[0], 1, diff0, 1, out0, 1);
                Vmath::Vvtvp (nqcol, m_df[1], 1, diff1, 1, out0, 1, out0, 1);
                Vmath::Vmul  (nqcol, m_df[2], 1, diff0, 1, out1, 1);
                Vmath::Vvtvp (nqcol, m_df[3], 1, diff1, 1, out1, 1, out1, 1);
                
                if (coordim == 3)
                {
                    Array<OneD, NekDouble> out2 (nqcol, output + 2*nqcol);
                    Vmath::Vmul  (nqcol, m_df[4], 1, diff0, 1, out2, 1);
                    Vmath::Vvtvp (nqcol, m_df[5], 1, diff1, 1, out2, 1, out2, 1);
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Quad)
            
            protected:
            const int m_nquad0;
            const int m_nquad1;
            Array<TwoD, NekDouble> m_df;
        };
        
        OperatorKey PhysDeriv_SumFac_Quad::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eQuadrilateral, ePhysDeriv, 
                                                eSumFac),
                                    PhysDeriv_SumFac_Quad::create, "PhysDeriv_SumFac_Quad");


    }
}
