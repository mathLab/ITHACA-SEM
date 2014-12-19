///////////////////////////////////////////////////////////////////////////////
//
// File: HexOperators.cpp
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
        //#define ELMTLOOP
        class BwdTrans_SumFac_Hex : public Operator
        {
        public:
            BwdTrans_SumFac_Hex(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (pCollExp[0]->GetNumPoints(0)),
                  m_nquad1  (pCollExp[0]->GetNumPoints(1)),
                  m_nquad2  (pCollExp[0]->GetNumPoints(2)),
                  m_nmodes0 (pCollExp[0]->GetBasisNumModes(0)),
                  m_nmodes1 (pCollExp[0]->GetBasisNumModes(1)),
                  m_nmodes2 (pCollExp[0]->GetBasisNumModes(2)),
                  m_base0   (pCollExp[0]->GetBasis(0)->GetBdata()),
                  m_base1   (pCollExp[0]->GetBasis(1)->GetBdata()),
                  m_base2   (pCollExp[0]->GetBasis(2)->GetBdata()),
                  m_colldir0(pCollExp[0]->GetBasis(0)->Collocation()),
                  m_colldir1(pCollExp[0]->GetBasis(1)->Collocation()),
                  m_colldir2(pCollExp[0]->GetBasis(2)->Collocation())
            {
#ifdef ELMTLOOP
                m_wspSize = 2*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
#else
                m_wspSize =  m_numElmt*m_nmodes0*(m_nmodes1*m_nquad2 + m_nquad1*m_nquad2);
#endif

            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                if(m_colldir0 && m_colldir1 && m_colldir2)
                {
                    Vmath::Vcopy(m_numElmt*m_nmodes0*m_nmodes1*m_nmodes2,input.get(),1,output.get(),1);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                    
                    // Assign second half of workspace for 2nd DGEMM operation.
                    int totmodes  = m_nmodes0*m_nmodes1*m_nmodes2;

#ifndef ELMTLOOP
                    Array<OneD, NekDouble> wsp2 = wsp + m_nmodes0*m_nmodes1*m_nquad2*m_numElmt; 

                    //loop over elements  and do bwd trans wrt c
                    for(int n = 0; n < m_numElmt; ++n)
                    {
                        Blas::Dgemm('N','T', m_nquad2, m_nmodes0*m_nmodes1,  m_nmodes2,
                                    1.0, m_base2.get(), m_nquad2, &input[n*totmodes], 
                                    m_nmodes0*m_nmodes1, 0.0, 
                                    &wsp[n*m_nquad2], m_nquad2*m_numElmt);
                    } 
                    
                    // trans wrt b 
                    Blas::Dgemm('N','T', m_nquad1, m_nquad2*m_numElmt*m_nmodes0,  
                                m_nmodes1, 1.0, m_base1.get(), m_nquad1,
                                wsp.get(),m_nquad2*m_numElmt*m_nmodes0,  
                                0.0, wsp2.get(), m_nquad1);

                    // trans wrt a
                    Blas::Dgemm('N','T', m_nquad0, m_nquad1*m_nquad2*m_numElmt,  
                                m_nmodes0, 1.0, m_base0.get(), m_nquad0,
                                wsp2.get(), m_nquad1*m_nquad2*m_numElmt,  
                                0.0, output.get(), m_nquad0);
#else
                    int totpoints = m_nquad0*m_nquad1*m_nquad2;
                    if(m_numElmt < m_nmodes0 || 1) // note sure what criterion we should use to swap around these strategies
                    {
                        Array<OneD, NekDouble> wsp2 = wsp + m_nmodes1*m_nmodes2*m_nquad0;

                        //loop over elements 
                        for(int n = 0; n < m_numElmt; ++n)
                        {
                            // BwdTrans in each direction using DGEMM
                            Blas::Dgemm('T','T', m_nmodes1*m_nmodes2, m_nquad0, m_nmodes0,
                                        1.0, &input[n*totmodes],   m_nmodes0,  m_base0.get(),   m_nquad0,
                                        0.0, &wsp[0], m_nmodes1*m_nmodes2);
                            
                            Blas::Dgemm('T','T', m_nquad0*m_nmodes2,  m_nquad1, m_nmodes1,
                                        1.0, &wsp[0],  m_nmodes1,  m_base1.get(), m_nquad1,
                                        0.0, &wsp2[0], m_nquad0*m_nmodes2);

                            Blas::Dgemm('T','T', m_nquad0*m_nquad1,   m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0], m_nmodes2, m_base2.get(), m_nquad2,
                                        0.0, &output[n*totpoints],  m_nquad0*<m_nquad1);
                        }
                    }
                    else
                    {
                        Array<OneD, NekDouble> wsp2 = wsp + m_numElmt*(max(totpoints,totmodes));
                        
                        // large degmm but copy at end. 
                        Blas::Dgemm('T','T', m_nmodes1*m_nmodes2*m_numElmt, m_nquad0, m_nmodes0,
                                    1.0, &input[0],   m_nmodes0,  m_base0.get(),   m_nquad0,
                                    0.0, &wsp[0],    m_nmodes1*m_nmodes2*m_numElmt);
                        
                        Blas::Dgemm('T','T', m_nmodes2*m_numElmt*m_nquad0,  m_nquad1, m_nmodes1,
                                    1.0, &wsp[0],   m_nmodes1, m_base1.get(),   m_nquad1,
                                    0.0, &wsp2[0],  m_nmodes2*m_numElmt*m_nquad0);

                        if(m_numElmt > 1)
                        {
                            Blas::Dgemm('T','T', m_numElmt*m_nquad0*m_nquad1, m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0],  m_nmodes2,  m_base2.get(),   m_nquad2,
                                        0.0, &wsp[0],  m_numElmt*m_nquad0*m_nquad1);
                            
                            for(int i = 0; i < totpoints; ++i)
                            {
                                Vmath::Vcopy(m_numElmt,&wsp[i*m_numElmt],1,&output[i],totpoints);
                            }
                        }
                        else
                        {
                            Blas::Dgemm('T','T', m_numElmt*m_nquad0*m_nquad1, m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0],  m_nmodes2,  m_base2.get(),   m_nquad2,
                                        0.0, &output[0],  m_numElmt*m_nquad0*m_nquad1);
                        }
                    }
#endif
                }
            }
            
            OPERATOR_CREATE(BwdTrans_SumFac_Hex)
            
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
            const bool m_colldir0;
            const bool m_colldir1;
            const bool m_colldir2;
        };
        
        OperatorKey BwdTrans_SumFac_Hex::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eSumFac,false),
                                    BwdTrans_SumFac_Hex::create, "BwdTrans_SumFac_Hex");
        

        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Hex : public Operator
        {
        public:
            IProductWRTBase_SumFac_Hex(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
                  m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
                  m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
                  m_colldir2(m_stdExp->GetBasis(2)->Collocation()),
                  m_base0    (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1    (m_stdExp->GetBasis(1)->GetBdata()),
                  m_base2    (m_stdExp->GetBasis(2)->GetBdata())
                
            {
                m_jac = GeomData->GetJacWithStdWeights(pCollExp);
                m_wspSize = 3*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD, NekDouble> &output,
                                    Array<OneD, NekDouble> &output1,
                                    Array<OneD, NekDouble> &output2,
                                    Array<OneD, NekDouble> &wsp)
            {

                ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                
                HexIProduct(m_colldir0,m_colldir1,m_colldir2, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_base2,
                            m_jac,input,output,wsp);
            }
            
            OPERATOR_CREATE(IProductWRTBase_SumFac_Hex)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const int  m_nmodes2;
            const bool m_colldir0;
            const bool m_colldir1;
            const bool m_colldir2;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_base2;
        };
        
        OperatorKey IProductWRTBase_SumFac_Hex::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron, eIProductWRTBase, eSumFac,false),
                                    IProductWRTBase_SumFac_Hex::create, "IProductWRTBase_SumFac_Hex");

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */
        
        class PhysDeriv_SumFac_Hex : public Operator
        {
        public:
            PhysDeriv_SumFac_Hex(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
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

#if 1 
                Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                            m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                            m_nquad0,0.0,&Diff[0][0],m_nquad0);
#endif
                for(int  i = 0; i < m_numElmt; ++i)
                {
#if 0
                    Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2,
                                m_nquad0,1.0, m_Deriv0,m_nquad0,&input[i*nPhys],
                                m_nquad0,0.0,&Diff[0][i*nPhys],m_nquad0);
#endif

                    for (int j = 0; j < m_nquad2; ++j)
                    {
                        Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                    1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0, m_Deriv1, m_nquad1, 0.0, 
                                    &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0);
                    }

                    Blas::Dgemm('N','T',m_nquad0*m_nquad1,m_nquad2,m_nquad2,
                                1.0, &input[i*nPhys],m_nquad0*m_nquad1,
                                m_Deriv2,m_nquad2, 0.0,&Diff[2][i*nPhys],
                                m_nquad0*m_nquad1);
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
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Hex)

            Array<TwoD, const NekDouble> m_derivFac;
            int m_dim;
            int m_coordim;
            const int m_nquad0;
            const int m_nquad1;
            const int m_nquad2;
            NekDouble *m_Deriv0;
            NekDouble *m_Deriv1;
            NekDouble *m_Deriv2;
        };

        OperatorKey PhysDeriv_SumFac_Hex::m_typeArr[] =
        {
            GetOperatorFactory().RegisterCreatorFunction(
                OperatorKey(LibUtilities::eHexahedron, ePhysDeriv, eSumFac,false),
                PhysDeriv_SumFac_Hex::create, "PhysDeriv_SumFac_Hex")
        };

        /*
         * ----------------------------------------------------------
         * IProductWRTDerivBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTDerivBase_SumFac_Hex : public Operator
        {
        public:
            IProductWRTDerivBase_SumFac_Hex(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pCollExp, GeomData),
                  m_nquad0  (m_stdExp->GetNumPoints(0)),
                  m_nquad1  (m_stdExp->GetNumPoints(1)),
                  m_nquad2  (m_stdExp->GetNumPoints(2)),
                  m_nmodes0 (m_stdExp->GetBasisNumModes(0)),
                  m_nmodes1 (m_stdExp->GetBasisNumModes(1)),
                  m_nmodes2 (m_stdExp->GetBasisNumModes(2)),
                  m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
                  m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
                  m_colldir2(m_stdExp->GetBasis(2)->Collocation()),
                  m_base0    (m_stdExp->GetBasis(0)->GetBdata()),
                  m_base1    (m_stdExp->GetBasis(1)->GetBdata()),
                  m_base2    (m_stdExp->GetBasis(2)->GetBdata()),
                  m_derbase0    (m_stdExp->GetBasis(0)->GetDbdata()),
                  m_derbase1    (m_stdExp->GetBasis(1)->GetDbdata()),
                  m_derbase2    (m_stdExp->GetBasis(2)->GetDbdata())
                
            {
                m_jac     = GeomData->GetJacWithStdWeights(pCollExp);
                m_wspSize = 6*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
                m_derivFac = GeomData->GetDerivFactors(pCollExp);
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
                
                // calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
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

                // calculate Iproduct WRT Std Deriv
                HexIProduct(false,m_colldir1,m_colldir2, m_numElmt,
                            m_nquad0,   m_nquad1,  m_nquad2,
                            m_nmodes0,  m_nmodes1, m_nmodes2,
                            m_derbase0, m_base1,   m_base2,
                            m_jac,tmp[0],output,wsp1);

                HexIProduct(m_colldir0,false,m_colldir2, m_numElmt,
                            m_nquad0,  m_nquad1,   m_nquad2,
                            m_nmodes0, m_nmodes1,  m_nmodes2,
                            m_base0,   m_derbase1, m_base2,
                            m_jac,tmp[1],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
                
                HexIProduct(m_colldir0,m_colldir1,false, m_numElmt,
                            m_nquad0,  m_nquad1,  m_nquad2,
                            m_nmodes0, m_nmodes1, m_nmodes2,
                            m_base0,   m_base1,   m_derbase2,
                            m_jac,tmp[2],tmp[0],wsp1);
                Vmath::Vadd(m_numElmt*nmodes,tmp[0],1,output,1,output,1);
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_SumFac_Hex)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nquad2;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const int  m_nmodes2;
            const bool m_colldir0;
            const bool m_colldir1;
            const bool m_colldir2;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
            Array<OneD, const NekDouble> m_base2;
            Array<OneD, const NekDouble> m_derbase0;
            Array<OneD, const NekDouble> m_derbase1;
            Array<OneD, const NekDouble> m_derbase2;
            Array<TwoD, const NekDouble> m_derivFac;
        };
        
        OperatorKey IProductWRTDerivBase_SumFac_Hex::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron,
                                                eIProductWRTDerivBase, eSumFac,false),
                                    IProductWRTDerivBase_SumFac_Hex::create, 
                                    "IProductWRTDerivBase_SumFac_Hex");

    }
}
