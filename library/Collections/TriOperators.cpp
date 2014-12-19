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
#include <Collections/IProduct.h>
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
            BwdTrans_SumFac_Tri(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata())
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
                

                int ncoeffs = m_stdExp->GetNcoeffs(); 
                int i,mode;
                for (i = mode = 0; i < m_nmodes0; ++i)
                {
                    
                    Blas::Dgemm('N','N', m_nquad1,m_numElmt,m_nmodes1-i,1.0, m_base1.get()+mode*m_nquad1,
                                m_nquad1, &input[0]+mode, ncoeffs,0.0,&wsp[i*m_nquad1*m_numElmt], m_nquad1);
                    mode += m_nmodes1-i;
                }
                
                // fix for modified basis by splitting top vertex mode
                if(m_sortTopVertex)
                {
                    for(i = 0; i < m_numElmt; ++i)
                    {
                        Blas::Daxpy(m_nquad1,input[1+i*ncoeffs],m_base1.get()+m_nquad1,1,
                                    &wsp[m_nquad1*m_numElmt]+i*m_nquad1,1);
                    }

                }
                
                Blas::Dgemm('N','T', m_nquad0,m_nquad1*m_numElmt,m_nmodes0,1.0, m_base0.get(),m_nquad0,
                            &wsp[0], m_nquad1*m_numElmt,0.0, &output[0], m_nquad0);
            } 
            
            OPERATOR_CREATE(BwdTrans_SumFac_Tri)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            bool m_sortTopVertex;
        };
        
        OperatorKey BwdTrans_SumFac_Tri::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTriangle, eBwdTrans, eSumFac,false),
                                    BwdTrans_SumFac_Tri::create, "BwdTrans_SumFac_Tri");


        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Tri : public Operator
        {
        public:
            IProductWRTBase_SumFac_Tri(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata())
            {
                m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
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
                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");

                TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1, 
                            m_nmodes0, m_nmodes1,m_base0,m_base1,m_jac, input,
                            output,wsp);
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
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eSumFac,false),IProductWRTBase_SumFac_Tri::create, "IProductWRTBase_SumFac_Tri");

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */       
        class PhysDeriv_SumFac_Tri : public Operator
        {
        public:
            PhysDeriv_SumFac_Tri(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator (pCollExp, GeomData),
                  m_nquad0 (m_stdExp->GetNumPoints(0)),
                  m_nquad1 (m_stdExp->GetNumPoints(1))
            {
                LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
                m_coordim = m_stdExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pCollExp);

                const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();
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

                
                m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
                m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
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
        
        OperatorKey PhysDeriv_SumFac_Tri::m_typeArr[] = 
            {
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                                ePhysDeriv, eSumFac,false),
                      PhysDeriv_SumFac_Tri::create, "PhysDeriv_SumFac_Tri"),

                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                   ePhysDeriv, eSumFac,true),
                      PhysDeriv_SumFac_Tri::create, "PhysDeriv_SumFac_NodalTri")
            };



        /*
         * ----------------------------------------------------------
         * IProductWRTDerivBase operator
         * ----------------------------------------------------------
         */       
        class IProductWRTDerivBase_SumFac_Tri : public Operator
        {
        public:
            IProductWRTDerivBase_SumFac_Tri(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                   CoalescedGeomDataSharedPtr GeomData)
                : Operator(pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
                  m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
                  m_derbase0   (m_stdExp->GetBasis(0)->GetDbdata()),
                  m_derbase1   (m_stdExp->GetBasis(1)->GetDbdata())
            {
                LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
                m_coordim = m_stdExp->GetCoordim();
                
                m_derivFac = GeomData->GetDerivFactors(pCollExp);
                m_jac      = GeomData->GetJacWithStdWeights(pCollExp);
                m_wspSize  = 4*m_numElmt*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1));

                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopVertex = true;
                }
                else
                {
                    m_sortTopVertex = false;
                }

                const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();
                
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
                // set up geometric factor: 2/(1-z1)
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        m_fac0[i+j*m_nquad0] = 2.0/(1-z1[j]);
                    }
                }

                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1);
                // set up geometric factor: (1+z0)/(1-z1)
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        m_fac1[i+j*m_nquad0] = (1+z0[i])/(1-z1[j]);
                    }
                }
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                                    Array<OneD, NekDouble> &entry1,
                                    Array<OneD, NekDouble> &entry2,
                                    Array<OneD, NekDouble> &entry3,
                                    Array<OneD, NekDouble> &wsp)
            {
                unsigned int nPhys  = m_stdExp->GetTotPoints();
                unsigned int ntot   = m_numElmt*nPhys;
                unsigned int nmodes = m_stdExp->GetNcoeffs();
                unsigned int nmax   = max(ntot,m_numElmt*nmodes);
                Array<OneD, Array<OneD, const NekDouble> > in(3);
                Array<OneD, NekDouble> output, wsp1;
                Array<OneD, Array<OneD, NekDouble> > tmp(2);
                
                in[0] = entry0; in[1] = entry1; in[2] = entry2; 
                
                output = (m_coordim == 2)? entry2: entry3;
                
                tmp[0] = wsp; tmp[1] = wsp + nmax;
                wsp1   = wsp + 2*nmax;
                
                
                // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx),in[0])
                //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy),in[1]) 
                //
                // Note dphi/dxi_0  = 
                //             dphi/deta_0 deta_0/dxi_0 = dphi/deta_0 2/(1-eta_1) 
                // 
                //      dphi/dxi_1  = 
                //             dphi/deta_1 deta_1/dxi_1 + dphi/deta_1 deta_1/dxi_1 = 
                //             dphi/deta_0 (1+eta_0)/(1-eta_1) + dphi/deta_1
                // 
                // and so the full inner products are 
                //
                // (dphi/dx,in[0]) + (dphi/dy,in[1])
                //    = (dphi/deta_0, ((2/(1-eta_1) (dxi_0/dx in[0] + dxi_0/dy in[1])
                //            + (1_eta_0)/(1-eta_1) (dxi_1/dx in[0] + dxi_1/dy in[1]))
                //    + (dphi/deta_1, (dxi_1/dx in[0] + dxi_1/dy in[1]))
                
                for(int i = 0; i < 2; ++i)
                {
                    Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1, tmp[i],1);

                    for(int j = 1; j < m_coordim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i +j*2],1,
                                      in[j],1, tmp[i], 1, tmp[i],1);
                    }
                }
            
                // Multiply by factor: 2/(1-z1)
                for (int i = 0; i < m_numElmt; ++i)
                {
                    // scale tmp[0] by geometric factor: 2/(1-z1)
                    Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                                tmp[0].get()+i*nPhys,1);
                    
                    // scale tmp[1] by geometric factor (1+z0)/(1-z1)
                    Vmath::Vvtvp(nPhys,&m_fac1[0],1,tmp[1].get()+i*nPhys,1,
                                 tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);
                }
                
                // Iproduct wrt derivative of base 0 
                TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1,
                            m_nmodes0,  m_nmodes1, m_derbase0, m_base1,
                            m_jac, tmp[0], output, wsp1);
                
                // Iproduct wrt derivative of base 1 
                TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1,
                            m_nmodes0,  m_nmodes1, m_base0, m_derbase1,
                            m_jac, tmp[1], tmp[0], wsp1);
                
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Tri)
            
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const bool m_colldir0;
            const bool m_colldir1;
            int m_coordim;
            Array<TwoD, const NekDouble> m_derivFac;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_derbase0;
            Array<OneD, const NekDouble> m_derbase1;
            Array<OneD, NekDouble> m_fac0;
            Array<OneD, NekDouble> m_fac1;
            bool m_sortTopVertex;
        };
        
        OperatorKey IProductWRTDerivBase_SumFac_Tri::m_type = 
            GetOperatorFactory().RegisterCreatorFunction(
                OperatorKey(LibUtilities::eTriangle, 
                            eIProductWRTDerivBase, eSumFac,false),
                IProductWRTDerivBase_SumFac_Tri::create, 
                "IProductWRTDerivBase_IterPerExp_Tri");
    }
}
