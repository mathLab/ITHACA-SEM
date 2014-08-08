///////////////////////////////////////////////////////////////////////////////
//
// File: TriOperators.cpp
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
#include <Collections/Collection.h>

namespace Nektar 
{
    namespace Collections 
    {

        /*
         * ----------------------------------------------------------
         * BwdTrans operators
         * ----------------------------------------------------------
         */       
        class BwdTrans_SumFac_Tri : public Operator
        {
        public:
            BwdTrans_SumFac_Tri(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1))
            {
                m_wspSize = m_nquad0*m_nmodes1*m_numElmt;
                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopVertex = true;
                }
                else
                {
                    m_sortTopVertex = false;
                }
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                Array<OneD, const NekDouble> base0  = m_stdExp->GetBasis(0)->GetBdata();
                Array<OneD, const NekDouble> base1  = m_stdExp->GetBasis(1)->GetBdata();
                

                int ncoeffs = m_stdExp->GetNcoeffs(); 
                int i,mode;
                for (i = mode = 0; i < m_nmodes0; ++i)
                {
                    
                    Blas::Dgemm('N','N', m_nquad1,m_numElmt,m_nmodes1-i,1.0, base1.get()+mode*m_nquad1,
                                m_nquad1, &input[0]+mode, ncoeffs,0.0,&wsp[i*m_nquad1*m_numElmt], m_nquad1);
                    mode += m_nmodes1-i;
                }
                
                // fix for modified basis by splitting top vertex mode
                if(m_sortTopVertex)
                {
                    for(i = 0; i < m_numElmt; ++i)
                    {
                        Blas::Daxpy(m_nquad1,input[1+i*ncoeffs],base1.get()+m_nquad1,1,
                                    &wsp[m_nquad1*m_numElmt]+i*m_nquad1,1);
                    }

                }
                
                Blas::Dgemm('N','T', m_nquad0,m_nquad1*m_numElmt,m_nmodes0,1.0, base0.get(),m_nquad0,
                            &wsp[0], m_nquad1*m_numElmt,0.0, &output[0], m_nquad0);
            } 
            
            OPERATOR_CREATE(BwdTrans_SumFac_Tri)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            bool m_sortTopVertex;
        };
        
        OperatorKey BwdTrans_SumFac_Tri::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTriangle, eBwdTrans, eSumFac),
                                    BwdTrans_SumFac_Tri::create, "BwdTrans_SumFac_Tri");


        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Tri : public Operator
        {
        public:
            IProductWRTBase_SumFac_Tri(StdRegions::StdExpansionSharedPtr pExp,
                                       vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1))
            {
                m_jac     = GeomData->GetJacWithStdWeights(pExp,pGeom);
                m_base0   = GeomData->GetBase(0,pExp);
                m_base1   = GeomData->GetBase(1,pExp);
                m_wspSize = 2*m_numElmt*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1));
                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopVertex = true;
                }
                else
                {
                    m_sortTopVertex = false;
                }
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                int totmodes  = m_stdExp->GetNcoeffs(); 
                int totpoints = m_nquad0 *m_nquad1;
                
                Vmath::Vmul(m_numElmt*totpoints,m_jac,1,input,1,wsp,1);

                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                Array<OneD, NekDouble> wsp1 = wsp  + max(totpoints,totmodes)*m_numElmt; 
                
                Blas::Dgemm('T','N', m_nquad1*m_numElmt,m_nmodes0,m_nquad0,1.0,
                            &wsp[0],m_nquad0, m_base0.get(), m_nquad0, 
                            0.0,&wsp1[0], m_nquad1*m_numElmt);
                
                int i, mode;
                // Inner product with respect to 'b' direction 
                for (mode=i=0; i < m_nmodes0; ++i)
                {
                    Blas::Dgemm('T','N',m_nmodes1-i,m_numElmt,m_nquad1,1.0,m_base1.get()+mode*m_nquad1,
                                m_nquad1,wsp1.get() + i*m_nquad1*m_numElmt,m_nquad1, 
                                0.0, &output[mode],totmodes);
                    
                    mode += m_nmodes1 - i;
                }

                // fix for modified basis by splitting top vertex mode
                if (m_sortTopVertex)
                {
                    Blas::Dgemv('T', m_nquad1,m_numElmt,1.0,wsp1.get()+m_nquad1*m_numElmt,m_nquad1,
                                m_base1.get()+m_nquad1,1,1.0, &output[1],totmodes);
                }
            }
            
            OPERATOR_CREATE(IProductWRTBase_SumFac_Tri)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            bool m_sortTopVertex;
        };
        
        OperatorKey IProductWRTBase_SumFac_Tri::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eSumFac),IProductWRTBase_SumFac_Tri::create, "IProductWRTBase_SumFac_Tri");

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */       
        class PhysDeriv_SumFac_Tri : public Operator
        {
        public:
            PhysDeriv_SumFac_Tri(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator (pExp, pGeom, GeomData),
                  m_nquad0 (pExp->GetNumPoints(0)),
                  m_nquad1 (pExp->GetNumPoints(1))
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_coordim = pExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);

                const Array<OneD, const NekDouble>& z0 = pExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z1 = pExp->GetBasis(1)->GetZ();
                
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
                // set up geometric factor: 0.5*(1+z0)
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        m_fac0[i+j*m_nquad0] = 0.5*(1+z0[i]);
                    }
                }

                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
                // set up geometric factor: 2/(1-z1)
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        m_fac1[i+j*m_nquad0] = 2.0/(1-z1[j]);
                    }
                }

                
                m_Deriv0 = &((pExp->GetBasis(0)->GetD())->GetPtr())[0];
                m_Deriv1 = &((pExp->GetBasis(1)->GetD())->GetPtr())[0];
                m_wspSize = 2 * m_nquad0*m_nquad1*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                const int nqtot   = m_nquad0 * m_nquad1;
                const int nqcol   = nqtot*m_numElmt;
                
                ASSERTL1(wsp.num_elements() == m_wspSize,
                         "Incorrect workspace size");
                ASSERTL1(input.num_elements() >= nqcol,
                         "Incorrect input size");
                
                Array<OneD, NekDouble> diff0(nqcol, wsp             );
                Array<OneD, NekDouble> diff1(nqcol, wsp    +   nqcol);
                
                // Tensor Product Derivative 
                Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt, 
                            m_nquad0, 1.0, m_Deriv0, m_nquad0, 
                            input.get(), m_nquad0, 0.0,
                            diff0.get(), m_nquad0);
                
                int cnt = 0;
                for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
                {
                    // scale diff0 by geometric factor: 2/(1-z1)
                    Vmath::Vmul(nqtot,&m_fac1[0],1,diff0.get()+cnt,1,
                                diff0.get()+cnt,1);
                    
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                                input.get() + cnt, m_nquad0, 
                                m_Deriv1, m_nquad1, 0.0,
                                diff1.get() + cnt, m_nquad0);

                    // add to diff1 by diff0 scaled by: (1_z0)/(1-z1)
                    Vmath::Vvtvp(nqtot,m_fac0.get(),1,diff0.get()+cnt,1,
                                 diff1.get()+cnt,1,diff1.get()+cnt,1);
                }

                
                Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[1], 1, diff1, 1, 
                              output0, 1, output0, 1);
                Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output1, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[3], 1, diff1, 1, 
                              output1, 1, output1, 1);
                
                if (m_coordim == 3)
                {
                    Vmath::Vmul  (nqcol, m_derivFac[4], 1, diff0, 1, 
                                  output2, 1);
                    Vmath::Vvtvp (nqcol, m_derivFac[5], 1, diff1, 1, 
                                  output2, 1, output2, 1);
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Tri)
            
            protected:
            int m_coordim;
            const int m_nquad0;
            const int m_nquad1;
            Array<TwoD, const NekDouble> m_derivFac;
            NekDouble *m_Deriv0;
            NekDouble *m_Deriv1;
            Array<OneD, NekDouble> m_fac0;
            Array<OneD, NekDouble> m_fac1;
        };
        
        OperatorKey PhysDeriv_SumFac_Tri::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTriangle, 
                                                ePhysDeriv, eSumFac),
                    PhysDeriv_SumFac_Tri::create, "PhysDeriv_SumFac_Tri");



    }
}
