///////////////////////////////////////////////////////////////////////////////
//
// File: PrismOperators.cpp
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
        class BwdTrans_SumFac_Prism : public Operator
        {
        public:
            BwdTrans_SumFac_Prism(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_nmodes2 (pExp->GetBasisNumModes(2)),
                  m_base0   (pExp->GetBasis(0)->GetBdata()),
                  m_base1   (pExp->GetBasis(1)->GetBdata()),
                  m_base2   (pExp->GetBasis(2)->GetBdata())
            {
                m_wspSize = m_numElmt*m_nmodes0*(m_nmodes1*m_nquad2 + m_nquad1*m_nquad2);

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
                
                // Assign second half of workspace for 2nd DGEMM operation.
                int totmodes  = m_stdExp->GetNcoeffs();
                
                Array<OneD, NekDouble> wsp2 = wsp + m_nmodes0*m_nmodes1*m_nquad2*m_numElmt;
                
                Vmath::Zero(m_nmodes0*m_nmodes1*m_nquad2*m_numElmt,wsp,1);
                int i,j,mode,mode1,cnt;
                for (i = mode = mode1 = 0; i < m_nmodes0; ++i)
                {
                    cnt = i*m_nquad2*m_numElmt; 
                    for (j = 0; j < m_nmodes1; ++j)
                    {
                        Blas::Dgemm('N','N',m_nquad2,m_numElmt,m_nmodes2-i,
                                    1.0,  m_base2.get()+mode*m_nquad2,
                                    m_nquad2, &input[0]+mode1, totmodes,0.0,
                                    &wsp[j*m_nquad2*m_numElmt*m_nmodes0+ cnt],
                                    m_nquad2);
                        mode1 += m_nmodes2-i;
                    }
                    mode += m_nmodes2-i;
                }

                // fix for modified basis by splitting top vertex mode
                if(m_sortTopVertex)
                {
                    for(j = 0; j < m_nmodes1; ++j)
                    {
                        for(i = 0; i < m_numElmt; ++i)
                        {
                            Blas::Daxpy(m_nquad2,input[1+i*totmodes+j*m_nmodes2],
                                        m_base2.get()+m_nquad2,1,
                                        &wsp[j*m_nquad2*m_numElmt*m_nmodes0 +
                                             m_nquad2*m_numElmt]+i*m_nquad2,1);
                        }
                    }
                    // Believe this could be made into a m_nmodes1
                    // dgemv if we made an array of m_numElmt copies
                    // of m_base2[m_quad2] (which are of size
                    // m_nquad2.
                }

                // Perform summation over '1' direction
                Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt*m_nmodes0, 
                            m_nmodes1,1.0, m_base1.get(),  m_nquad1,
                            wsp.get(), m_nquad2*m_numElmt*m_nmodes0,
                            0.0, wsp2.get(), m_nquad1);


                // Perform summation over '0' direction
                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt, 
                            m_nmodes0, 1.0, m_base0.get(),  m_nquad0,
                            wsp2.get(), m_nquad1*m_nquad2*m_numElmt,
                            0.0, output.get(), m_nquad0);
                
            }
            
            OPERATOR_CREATE(BwdTrans_SumFac_Prism)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const int  m_nmodes2;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_base2;
            bool m_sortTopVertex;
        };
        
        OperatorKey BwdTrans_SumFac_Prism::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::ePrism, 
                                                eBwdTrans, eSumFac),
                                    BwdTrans_SumFac_Prism::create, "BwdTrans_SumFac_Prism");
        
        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Prism : public Operator
        {
        public:
            IProductWRTBase_SumFac_Prism(StdRegions::StdExpansionSharedPtr pExp,
                                       vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_nmodes2 (pExp->GetBasisNumModes(2)),
                  m_base0    (pExp->GetBasis(0)->GetBdata()),
                  m_base1    (pExp->GetBasis(1)->GetBdata()),
                  m_base2    (pExp->GetBasis(2)->GetBdata())
                
            {
                m_jac = GeomData->GetJacWithStdWeights(pExp,pGeom);

                m_wspSize = m_numElmt*m_nquad2*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1)) + 
                    m_nquad1*m_nquad2*m_numElmt*m_nmodes0;

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
                                    Array<OneD, NekDouble> &output,
                                    Array<OneD, NekDouble> &output1,
                                    Array<OneD, NekDouble> &output2,
                                    Array<OneD, NekDouble> &wsp)
            {

                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                PrismIProduct(m_sortTopVertex, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_base2,
                            m_jac,input,output,wsp);
            }
            
            OPERATOR_CREATE(IProductWRTBase_SumFac_Prism)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const int  m_nmodes2;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_base2;
            bool m_sortTopVertex;
        };
        
        OperatorKey IProductWRTBase_SumFac_Prism::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eSumFac),
                                    IProductWRTBase_SumFac_Prism::create, "IProductWRTBase_SumFac_Prism");

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */
        
        class PhysDeriv_SumFac_Prism : public Operator
        {
        public:
            PhysDeriv_SumFac_Prism(StdRegions::StdExpansionSharedPtr pExp,
                                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2))
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();

                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);

                const Array<OneD, const NekDouble>& z0 = pExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z2 = pExp->GetBasis(2)->GetZ();
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        for(int k = 0; k < m_nquad2; ++k)
                        {
                            m_fac0[i+j*m_nquad0 + k*m_nquad0*m_nquad1] = 
                                2.0/(1-z2[k]);
                            m_fac1[i+j*m_nquad0 + k*m_nquad0*m_nquad1] =
                                0.5*(1+z0[i]);
                        }
                    }
                }



                m_Deriv0 = &((pExp->GetBasis(0)->GetD())->GetPtr())[0];
                m_Deriv1 = &((pExp->GetBasis(1)->GetD())->GetPtr())[0];
                m_Deriv2 = &((pExp->GetBasis(2)->GetD())->GetPtr())[0];

                m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                int nPhys = m_stdExp->GetTotPoints();
                int ntot = m_numElmt*nPhys;
                Array<OneD, NekDouble> tmp0,tmp1,tmp2;
                Array<OneD, Array<OneD, NekDouble> > Diff(3);
                Array<OneD, Array<OneD, NekDouble> > out(3);
                out[0] = output0; out[1] = output1; out[2] = output2;

                for(int i = 0; i < m_dim; ++i)
                {
                    Diff[i] = wsp + i*ntot;
                }

                // dEta0
                Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                            m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                            m_nquad0,0.0,&Diff[0][0],m_nquad0);
                
                int cnt = 0; 
                for(int  i = 0; i < m_numElmt; ++i)
                {
                    
                    // dEta 1
                    for (int j = 0; j < m_nquad2; ++j)
                    {
                        Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                    1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0, m_Deriv1, m_nquad1, 0.0, 
                                    &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0);
                    }

                    // dEta 2
                    Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                                1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                                m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                                m_nquad0*m_nquad1);

                    // dxi0 = 2/(1-eta_2) d Eta_0
                    Vmath::Vmul(nPhys,&m_fac0[0],1,Diff[0].get()+cnt,1,
                                Diff[0].get()+cnt,1);

                    // dxi2 = (1+eta0)/(1-eta_2) d Eta_0 + d/dEta2;
                    Vmath::Vvtvp(nPhys,&m_fac1[0],1,Diff[0].get()+cnt,1,
                                 Diff[2].get()+cnt,1,Diff[2].get()+cnt,1);
                    cnt += nPhys; 
                }

                // calculate full derivative 
                for(int i = 0; i < m_coordim; ++i)
                {
                    Vmath::Vmul(ntot,m_derivFac[i*m_dim],1,Diff[0],1,out[i],1);
                    for(int j = 1; j < m_dim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i*m_dim+j],1,Diff[j],1, out[i], 1, out[i],1);
                    }
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Prism)

            Array<TwoD, const NekDouble> m_derivFac;
            int m_dim;
            int m_coordim;
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            NekDouble *m_Deriv0;
            NekDouble *m_Deriv1;
            NekDouble *m_Deriv2;
            Array<OneD, NekDouble> m_fac0;
            Array<OneD, NekDouble> m_fac1;
        };

        OperatorKey PhysDeriv_SumFac_Prism::m_typeArr[] =
        {
            GetOperatorFactory().RegisterCreatorFunction(
                OperatorKey(LibUtilities::ePrism, ePhysDeriv, eSumFac),
                PhysDeriv_SumFac_Prism::create, "PhysDeriv_SumFac_Prism")
        };

        /*
         * ----------------------------------------------------------
         * IProductWRTDerivBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTDerivBase_SumFac_Prism : public Operator
        {
        public:
            IProductWRTDerivBase_SumFac_Prism(StdRegions::StdExpansionSharedPtr pExp,
                                       vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_nmodes2 (pExp->GetBasisNumModes(2)),
                  m_base0    (pExp->GetBasis(0)->GetBdata()),
                  m_base1    (pExp->GetBasis(1)->GetBdata()),
                  m_base2    (pExp->GetBasis(2)->GetBdata()),
                  m_derbase0    (pExp->GetBasis(0)->GetDbdata()),
                  m_derbase1    (pExp->GetBasis(1)->GetDbdata()),
                  m_derbase2    (pExp->GetBasis(2)->GetDbdata())
                
            {
                m_jac      = GeomData->GetJacWithStdWeights(pExp,pGeom);
                m_wspSize  = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);

                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopVertex = true;
                }
                else
                {
                    m_sortTopVertex = false;
                }

                const Array<OneD, const NekDouble>& z0 = pExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z2 = pExp->GetBasis(2)->GetZ();
                
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);

                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        for(int k = 0; k < m_nquad2; ++k)
                        {
                            // set up geometric factor: 2/(1-z1)
                            m_fac0[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0/(1-z2[k]);
                            // set up geometric factor: (1+z0)/(1-z1)
                            m_fac1[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z0[i])/(1-z2[k]);

                        }
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
                unsigned int nmax  = max(ntot,m_numElmt*nmodes);
                Array<OneD, Array<OneD, const NekDouble> > in(3);
                Array<OneD, NekDouble> output, wsp1;
                Array<OneD, Array<OneD, NekDouble> > tmp(3);

                in[0] = entry0; in[1] = entry1; 
                in[2] = entry2; 
                
                output =  entry3;

                for(int i = 0; i < 3; ++i)
                {
                    tmp[i] = wsp + i*nmax; 
                }
                
                // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx),in[0])
                //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy),in[1]) 
                //     +     (dphi/dz,in[2]) = ((dphi/dxi_0 dxi_0/dz + dphi/dxi_1 dxi_1/dz),in[2]) 
                //
                // Note dphi/dxi_0  = 
                //             dphi/deta_0 deta_0/dxi_0 = dphi/deta_0 2/(1-eta_2) 
                // 
                //      dphi/dxi_2  = 
                //             dphi/deta_0 deta_0/dxi_2 + dphi/deta_2 deta_2/dxi_2 = 
                //             dphi/deta_0 (1+eta_0)/(1-eta_2) + dphi/deta_2
                // 
                // and so the full inner products are 
                //
                // (dphi/dx,in[0]) + (dphi/dy,in[1]) + (dphi/dz,in[2])
                //    = (dphi/deta_0, ((2/(1-eta_2) (dxi_0/dx in[0] + dxi_0/dy in[1] + dxi_0/dz in[2]   )
                //            + (1_eta_0)/(1-eta_2) (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2] ))
                //    + (dphi/deta_1, (dxi_1/dx in[0] + dxi_1/dy in[1] + dxi_1/dz in[2]))
                //    + (dphi/deta_2, (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2]))



                for(int i = 0; i < 3; ++i)
                {
                    Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1, 
                                 tmp[i],1);
                    for(int j = 1; j < 3; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i+3*j],1,
                                      in[j],1, tmp[i], 1, tmp[i],1);
                    }
                }
                wsp1   = wsp + 3*nmax;

                // Sort into eta factors 
                for (int i = 0; i < m_numElmt; ++i)
                {
                    // scale tmp[0] by fac0 
                    Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                                tmp[0].get()+i*nPhys,1);
                    
                    // scale tmp[2] by fac1 and add to tmp0
                    Vmath::Vvtvp(nPhys,&m_fac1[0],1,tmp[2].get()+i*nPhys,1,
                                 tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);
                }
                
                // calculate Iproduct WRT Std Deriv
                PrismIProduct(m_sortTopVertex, m_numElmt,
                            m_nquad0,   m_nquad1,  m_nquad2,
                            m_nmodes0,  m_nmodes1, m_nmodes2,
                            m_derbase0, m_base1,   m_base2,
                            m_jac,tmp[0],output,wsp1);

                PrismIProduct(m_sortTopVertex, m_numElmt,
                            m_nquad0,  m_nquad1,   m_nquad2,
                            m_nmodes0, m_nmodes1,  m_nmodes2,
                            m_base0,   m_derbase1, m_base2,
                            m_jac,tmp[1],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
                
                PrismIProduct(m_sortTopVertex, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_derbase2,
                            m_jac,tmp[2],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Prism)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const int  m_nmodes2;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_base2;
            Array<OneD, const NekDouble> m_derbase0;
            Array<OneD, const NekDouble> m_derbase1;
            Array<OneD, const NekDouble> m_derbase2;
            Array<TwoD, const NekDouble> m_derivFac;
            Array<OneD, NekDouble> m_fac0;
            Array<OneD, NekDouble> m_fac1;
            bool m_sortTopVertex;
        };
        
        OperatorKey IProductWRTDerivBase_SumFac_Prism::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::ePrism,
                                                eIProductWRTDerivBase, eSumFac),
                                    IProductWRTDerivBase_SumFac_Prism::create, 
                                    "IProductWRTDerivBase_SumFac_Prism");
        
    }
}
