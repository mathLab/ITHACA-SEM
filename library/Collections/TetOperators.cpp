///////////////////////////////////////////////////////////////////////////////
//
// File: TetOperators.cpp
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
        class BwdTrans_SumFac_Tet : public Operator
        {
        public:
            BwdTrans_SumFac_Tet(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
                  m_base2   (m_stdExp->GetBasis(2)->GetBdata())
            {
                m_wspSize = m_numElmt*(m_nquad2*m_nmodes0*(2*m_nmodes1-m_nmodes0+1)/2+
                                       m_nquad2*m_nquad1*m_nmodes0);

                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopEdge = true;
                }
                else
                {
                    m_sortTopEdge = false;
                }
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                Array<OneD, NekDouble > tmp  = wsp;
                Array<OneD, NekDouble > tmp1 = tmp + m_numElmt*m_nquad2*m_nmodes0*(2*m_nmodes1-m_nmodes0+1)/2;
                
                int mode, mode1, cnt;
                int ncoeffs = m_stdExp->GetNcoeffs();
                
                // Perform summation over '2' direction
                mode = mode1 = cnt = 0;
                
                for(int i = 0; i < m_nmodes0; ++i)
                {
                    for(int j = 0; j < m_nmodes1-i; ++j, ++cnt)
                    {
                        Blas::Dgemm('N', 'N', m_nquad2, m_numElmt, m_nmodes2-i-j,
                                    1.0, m_base2.get()+mode*m_nquad2, m_nquad2,
                                    input.get()+ mode1,  ncoeffs, 0.0, tmp.get() + 
                                    cnt*m_nquad2*m_numElmt, m_nquad2);
                        mode  += m_nmodes2-i-j;
                        mode1 += m_nmodes2-i-j;
                    }
                    
                    //increment mode in case m_nmodes1!=m_nmodes2
                    mode +=  (m_nmodes2-m_nmodes1)*(m_nmodes2-m_nmodes1+1)/2;
                }

                // vertex mode - currently (1+c)/2 x (1-b)/2 x (1-a)/2
                // component is evaluated
                if(m_sortTopEdge)
                {
                    for(int n = 0; n < m_numElmt; ++n)
                    {
                        // top singular vertex - (1+c)/2 x (1+b)/2 x (1-a)/2 component
                        Blas::Daxpy(m_nquad2,input[1+n*ncoeffs],m_base2.get()+m_nquad2,1,
                                    &tmp[m_nquad2*m_numElmt]+n*m_nquad2,1);
                        
                        // top singular vertex - (1+c)/2 x (1-b)/2 x (1+a)/2 component
                        Blas::Daxpy(m_nquad2,input[1+n*ncoeffs],m_base2.get()+m_nquad2,1,
                                    &tmp[m_nmodes1*m_nquad2*m_numElmt]+n*m_nquad2,1);
                    }
                }

                // Perform summation over '1' direction
                mode = 0;
                for(int i = 0; i < m_nmodes0; ++i)
                {
                    Blas::Dgemm('N', 'T', m_nquad1, m_nquad2*m_numElmt, m_nmodes1-i,
                                1.0, m_base1.get()+mode*m_nquad1,  m_nquad1,
                                tmp.get()+mode*m_nquad2*m_numElmt, m_nquad2*m_numElmt,
                                0.0, tmp1.get()+i*m_nquad1*m_nquad2*m_numElmt, m_nquad1);
                    mode  += m_nmodes1-i;
                }

                // fix for modified basis by adding additional split of
                // top and base singular vertex modes as well as singular
                // edge
                if(m_sortTopEdge)
                {
                    // this could probably be a dgemv or higher if we
                    // made a specialised m_base1[m_nuqad1] array
                    // containing multiply copies
                    for(int n = 0; n < m_numElmt; ++n)
                    {
                        // sort out singular vertices and singular
                        // edge components with (1+b)/2 (1+a)/2 form
                        for(int i = 0; i < m_nquad2; ++i)
                        {
                            Blas::Daxpy(m_nquad1,tmp[m_nquad2*m_numElmt+n*m_nquad2+i], 
                                        m_base1.get()+m_nquad1,1,
                                        &tmp1[m_nquad1*m_nquad2*m_numElmt]+n*m_nquad1*m_nquad2+i*m_nquad1,1);
                        }
                    }
                }
                    
                // Perform summation over '0' direction
                Blas::Dgemm('N', 'T', m_nquad0, m_nquad1*m_nquad2*m_numElmt, m_nmodes0,
                            1.0, m_base0.get(),    m_nquad0,  tmp1.get(),     
                            m_nquad1*m_nquad2*m_numElmt, 0.0, output.get(), m_nquad0);
                
            }
            OPERATOR_CREATE(BwdTrans_SumFac_Tet)
            
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
            bool m_sortTopEdge;
        };
        
        OperatorKey BwdTrans_SumFac_Tet::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eSumFac,false),
                                    BwdTrans_SumFac_Tet::create, "BwdTrans_SumFac_Tet");


        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Tet : public Operator
        {
        public:
            IProductWRTBase_SumFac_Tet(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
                  m_base2   (m_stdExp->GetBasis(2)->GetBdata())
            {
                m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
                m_wspSize = m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,
                            m_nquad2*m_nmodes0*(2*m_nmodes1-m_nmodes0+1)/2)+
                                       m_nquad2*m_nquad1*m_nmodes0);

                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopEdge = true;
                }
                else
                {
                    m_sortTopEdge = false;
                }
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                TetIProduct(m_sortTopEdge, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_base2,
                            m_jac,input,output,wsp);
                
            }
            OPERATOR_CREATE(IProductWRTBase_SumFac_Tet)
            
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
            bool m_sortTopEdge;
        };
        
        OperatorKey IProductWRTBase_SumFac_Tet::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eSumFac,false), IProductWRTBase_SumFac_Tet::create, "IProductWRTBase_SumFac_Tet");


        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */
        
        class PhysDeriv_SumFac_Tet : public Operator
        {
        public:
            PhysDeriv_SumFac_Tet(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2))
            {
                LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();

                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pCollExp);

                m_Deriv0 = &((m_stdExp->GetBasis(0)->GetD())->GetPtr())[0];
                m_Deriv1 = &((m_stdExp->GetBasis(1)->GetD())->GetPtr())[0];
                m_Deriv2 = &((m_stdExp->GetBasis(2)->GetD())->GetPtr())[0];

                m_wspSize = 3*m_nquad0*m_nquad1*m_nquad2*m_numElmt;

                const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();
                const Array<OneD, const NekDouble>& z2 = m_stdExp->GetBasis(2)->GetZ();
                
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac2 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac3 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                // calculate 2.0/((1-eta_1)(1-eta_2))
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        for(int k = 0; k < m_nquad2; ++k)
                        {
                            
                            m_fac0[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 4.0/((1-z1[j])*(1-z2[k]));
                            m_fac1[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0*(1+z0[i])/((1-z1[j])*(1-z2[k]));
                            m_fac2[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0/(1-z2[k]);
                            m_fac3[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z1[j])/(1-z2[k]);
                        }
                    }
                }
            
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
                out[0] = output0;  out[1] = output1;    out[2] = output2;

                for(int i = 0; i < m_dim; ++i)
                {
                    Diff[i] = wsp + i*ntot;
                }

                // dEta0
                Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                            m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                            m_nquad0,0.0,&Diff[0][0],m_nquad0);

                // dEta2
                for(int  i = 0; i < m_numElmt; ++i)
                {
                    Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                                1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                                m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                                m_nquad0*m_nquad1);
                }
                
                for(int  i = 0; i < m_numElmt; ++i)
                {
                
                    // dEta1 
                    for (int j = 0; j < m_nquad2; ++j)
                    {
                        Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                    1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0, m_Deriv1, m_nquad1, 0.0, 
                                    &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0);
                    }
                    
                    // dxi2 = (1 + eta_1)/(1 -eta_2)*dEta1 + dEta2
                    Vmath::Vvtvp(nPhys,m_fac3.get(),1,Diff[1].get()+i*nPhys,1,Diff[2].get()
                                 +i*nPhys,1,Diff[2].get()+i*nPhys,1);

                    // dxi1 =  2/(1 - eta_2) dEta1
                    Vmath::Vmul(nPhys,m_fac2.get(),1,Diff[1].get()+i*nPhys,1,Diff[1].get()+i*nPhys,1);
                    
                    // dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi1
                    Vmath::Vvtvp(nPhys,m_fac1.get(),1,Diff[0].get()+i*nPhys,1,Diff[1].get()+i*nPhys,
                                 1,Diff[1].get()+i*nPhys,1);
                    
                    // dxi2 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) dEta0 + dxi2
                    Vmath::Vvtvp(nPhys,m_fac1.get(),1,Diff[0].get()+i*nPhys,1,Diff[2].get()+i*nPhys,
                                 1,Diff[2].get()+i*nPhys,1);

                    // dxi0 = 4.0/((1-eta_1)(1-eta_2)) dEta0
                    Vmath::Vmul(nPhys,m_fac0.get(),1,Diff[0].get()+i*nPhys,1,Diff[0].get()+i*nPhys,1);
                    
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
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Tet)

            Array<TwoD, const NekDouble> m_derivFac;
            int m_dim;
            int m_coordim;
            const int m_nquad0;
            const int m_nquad1;
            const int m_nquad2;
            NekDouble *m_Deriv0;
            NekDouble *m_Deriv1;
            NekDouble *m_Deriv2;
            Array<OneD, NekDouble> m_fac0;
            Array<OneD, NekDouble> m_fac1;
            Array<OneD, NekDouble> m_fac2;
            Array<OneD, NekDouble> m_fac3;
        };

        OperatorKey PhysDeriv_SumFac_Tet::m_typeArr[] =
        {
            GetOperatorFactory().RegisterCreatorFunction(
                OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eSumFac,false),
                PhysDeriv_SumFac_Tet::create, "PhysDeriv_SumFac_Tet")
        };



        /*
         * ----------------------------------------------------------
         * IProductWRTDerivBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTDerivBase_SumFac_Tet : public Operator
        {
        public:
            IProductWRTDerivBase_SumFac_Tet(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                            CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
                  m_base0   (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1   (m_stdExp->GetBasis(1)->GetBdata()),
                  m_base2   (m_stdExp->GetBasis(2)->GetBdata()),
                  m_derbase0(m_stdExp->GetBasis(0)->GetDbdata()),
                  m_derbase1(m_stdExp->GetBasis(1)->GetDbdata()),
                  m_derbase2(m_stdExp->GetBasis(2)->GetDbdata())
                
            {
                m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
                m_wspSize = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
                m_derivFac = GeomData->GetDerivFactors(pCollExp);


                const Array<OneD, const NekDouble>& z0 = m_stdExp->GetBasis(0)->GetZ();
                const Array<OneD, const NekDouble>& z1 = m_stdExp->GetBasis(1)->GetZ();
                const Array<OneD, const NekDouble>& z2 = m_stdExp->GetBasis(2)->GetZ();
                
                m_fac0 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac1 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac2 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                m_fac3 = Array<OneD, NekDouble>(m_nquad0*m_nquad1*m_nquad2);
                // calculate 2.0/((1-eta_1)(1-eta_2))
                for (int i = 0; i < m_nquad0; ++i)
                {
                    for(int j = 0; j < m_nquad1; ++j)
                    {
                        for(int k = 0; k < m_nquad2; ++k)
                        {
                            
                            m_fac0[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 4.0/((1-z1[j])*(1-z2[k]));
                            m_fac1[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z0[i])*0.5;
                            m_fac2[i+j*m_nquad0+k*m_nquad0*m_nquad1] = 2.0/(1-z2[k]);
                            m_fac3[i+j*m_nquad0+k*m_nquad0*m_nquad1] = (1+z1[j])*0.5;
                        }
                    }
                }
            
                if(m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
                {
                    m_sortTopEdge = true;
                }
                else
                {
                    m_sortTopEdge = false;
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
                

                // calculate (dphi/dx,in[0]) = ((dphi/dxi_0 dxi_0/dx + dphi/dxi_1 dxi_1/dx + dphi/dxi_2 dxi_2/dx),in[0])
                //     +     (dphi/dy,in[1]) = ((dphi/dxi_0 dxi_0/dy + dphi/dxi_1 dxi_1/dy + dphi/dxi_2 dxi_2/dy),in[1]) 
                //     +     (dphi/dz,in[2]) = ((dphi/dxi_0 dxi_0/dz + dphi/dxi_1 dxi_1/dz + dphi/dxi_2 dxi_2/dz),in[1]) 
                //
                // Note dphi/dxi_0  = 
                //             dphi/deta_0 4/((1-eta_1)(1-eta2))
                // 
                //      dphi/dxi_1  = 
                //             dphi/deta_0 2(1+eta_0)/((1-eta_1)(1-eta_2)) + dphi/deta_1 2/(1-eta_2)
                // 
                //      dphi/dxi_2  = 
                //             dphi/deta_0 2(1+eta_0)/((1-eta_1)(1-eta_2)) + dphi/deta_1 (1+eta_1)/(1-eta_2)  + dphi/deta_2
                //
                // and so the full inner products are 
                //
                // (dphi/dx,in[0]) + (dphi/dy,in[1]) + (dphi/dz,in[2])  
                //    = (dphi/deta_0, fac0 (tmp0 + fac1(tmp1 + tmp2)))
                //    + (dphi/deta_1, fac2 (tmp1 + fac3 tmp2))
                //    + (dphi/deta_2, tmp2)
                //
                // tmp0 = (dxi_0/dx in[0] + dxi_0/dy in[1] + dxi_0/dz in[2])
                // tmp1 = (dxi_1/dx in[0] + dxi_1/dy in[1] + dxi_1/dz in[2])
                // tmp2 = (dxi_2/dx in[0] + dxi_2/dy in[1] + dxi_2/dz in[2])

                // fac0 = 4/((1-eta_1)(1-eta2))
                // fac1 = (1+eta_0)/2
                // fac2 = 2/(1-eta_2)
                // fac3 = (1+eta_1)/2

                for(int i = 0; i < 3; ++i)
                {
                    Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1,tmp[i],1);
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
                    // add tmp[1] + tmp[2] 
                    Vmath::Vadd(nPhys,tmp[1].get()+i*nPhys,1,tmp[2].get()+i*nPhys,1,wsp1.get(),1);

                    // scale wsp1 by fac1 and add to tmp0
                    Vmath::Vvtvp(nPhys,&m_fac1[0],1,wsp1.get(),1,
                                 tmp[0].get()+i*nPhys,1,tmp[0].get()+i*nPhys,1);

                    // scale tmp[0] by fac0 
                    Vmath::Vmul(nPhys,&m_fac0[0],1,tmp[0].get()+i*nPhys,1,
                                tmp[0].get()+i*nPhys,1);
                    
                    // scale tmp[2] by fac3 and add to tmp1
                    Vmath::Vvtvp(nPhys,&m_fac3[0],1,tmp[2].get()+i*nPhys,1,
                                 tmp[1].get()+i*nPhys,1,tmp[1].get()+i*nPhys,1);

                    // scale tmp[1] by fac2
                    Vmath::Vmul(nPhys,&m_fac2[0],1,tmp[1].get()+i*nPhys,1,
                                tmp[1].get()+i*nPhys,1);
                }
                
                
                // calculate Iproduct WRT Std Deriv
                TetIProduct(m_sortTopEdge, m_numElmt,
                            m_nquad0,   m_nquad1,  m_nquad2,
                            m_nmodes0,  m_nmodes1, m_nmodes2,
                            m_derbase0, m_base1,   m_base2,
                            m_jac,tmp[0],output,wsp1);

                TetIProduct(m_sortTopEdge, m_numElmt,
                            m_nquad0,  m_nquad1,   m_nquad2,
                            m_nmodes0, m_nmodes1,  m_nmodes2,
                            m_base0,   m_derbase1, m_base2,
                            m_jac,tmp[1],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
                
                TetIProduct(m_sortTopEdge, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_derbase2,
                            m_jac,tmp[2],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Tet)
            
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
            Array<OneD, NekDouble> m_fac2;
            Array<OneD, NekDouble> m_fac3;
            bool m_sortTopEdge;
        };
        
        OperatorKey IProductWRTDerivBase_SumFac_Tet::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eTetrahedron,
                                                eIProductWRTDerivBase, eSumFac,false),
                                    IProductWRTDerivBase_SumFac_Tet::create, 
                                    "IProductWRTDerivBase_SumFac_Tet");


    }
}
