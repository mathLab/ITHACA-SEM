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

namespace Nektar 
{
    namespace Collections 
    {

        /*
         * ----------------------------------------------------------
         * BwdTrans operators
         * ----------------------------------------------------------
         */       
        class BwdTrans_SumFac_Hex : public Operator
        {
        public:
            BwdTrans_SumFac_Hex(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_nmodes2 (pExp->GetBasisNumModes(2)),
                  m_colldir0(pExp->GetBasis(0)->Collocation()),
                  m_colldir1(pExp->GetBasis(1)->Collocation()),
                  m_colldir2(pExp->GetBasis(2)->Collocation())
            {
                m_wspSize = 2*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                Array<OneD, const NekDouble> base0  = m_stdExp->GetBasis(0)->GetBdata();
                Array<OneD, const NekDouble> base1  = m_stdExp->GetBasis(1)->GetBdata();
                Array<OneD, const NekDouble> base2  = m_stdExp->GetBasis(2)->GetBdata();


                if(m_colldir0 && m_colldir1 && m_colldir2)
                {
                    Vmath::Vcopy(m_numElmt*m_nmodes0*m_nmodes1*m_nmodes2,input.get(),1,output.get(),1);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                    
                    // Assign second half of workspace for 2nd DGEMM operation.
                    int totmodes  = m_nmodes0*m_nmodes1*m_nmodes2;
                    int totpoints = m_nquad0*m_nquad1*m_nquad2;

                    if(m_numElmt < m_nmodes0) // note sure what criterion we shoudl use to swap around these strategies
                    {
                        Array<OneD, NekDouble> wsp2 = wsp + m_nmodes1*m_nmodes2*m_nquad0;

                        //loop over elements 
                        for(int n = 0; n < m_numElmt; ++n)
                        {
                            // BwdTrans in each direction using DGEMM
                            Blas::Dgemm('T','T', m_nmodes1*m_nmodes2, m_nquad0, m_nmodes0,
                                        1.0, &input[n*totmodes],   m_nmodes0,  base0.get(),   m_nquad0,
                                        0.0, &wsp[0], m_nmodes1*m_nmodes2);
                            
                            Blas::Dgemm('T','T', m_nquad0*m_nmodes2,  m_nquad1, m_nmodes1,
                                        1.0, &wsp[0],  m_nmodes1,  base1.get(), m_nquad1,
                                        0.0, &wsp2[0], m_nquad0*m_nmodes2);

                            Blas::Dgemm('T','T', m_nquad0*m_nquad1,   m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0], m_nmodes2, base2.get(), m_nquad2,
                                        0.0, &output[n*totpoints],  m_nquad0*m_nquad1);
                        }
                    }
                    else
                    {
                        Array<OneD, NekDouble> wsp2 = wsp + m_numElmt*(max(totpoints,totmodes));
                        
                        // large degmm but copy at end. 
                        Blas::Dgemm('T','T', m_nmodes1*m_nmodes2*m_numElmt, m_nquad0, m_nmodes0,
                                    1.0, &input[0],   m_nmodes0,  base0.get(),   m_nquad0,
                                    0.0, &wsp[0],    m_nmodes1*m_nmodes2*m_numElmt);
                        
                        Blas::Dgemm('T','T', m_nmodes2*m_numElmt*m_nquad0,  m_nquad1, m_nmodes1,
                                    1.0, &wsp[0],   m_nmodes1, base1.get(),   m_nquad1,
                                    0.0, &wsp2[0],  m_nmodes2*m_numElmt*m_nquad0);

                        if(m_numElmt > 1)
                        {
                            Blas::Dgemm('T','T', m_numElmt*m_nquad0*m_nquad1, m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0],  m_nmodes2,  base2.get(),   m_nquad2,
                                        0.0, &wsp[0],  m_numElmt*m_nquad0*m_nquad1);
                            
                            for(int i = 0; i < totpoints; ++i)
                            {
                                Vmath::Vcopy(m_numElmt,&wsp[i*m_numElmt],1,&output[i],totpoints);
                            }
                        }
                        else
                        {
                            Blas::Dgemm('T','T', m_numElmt*m_nquad0*m_nquad1, m_nquad2, m_nmodes2,
                                        1.0, &wsp2[0],  m_nmodes2,  base2.get(),   m_nquad2,
                                        0.0, &output[0],  m_numElmt*m_nquad0*m_nquad1);
                        }
                    }
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
            const bool m_colldir0;
            const bool m_colldir1;
            const bool m_colldir2;
        };
        
        OperatorKey BwdTrans_SumFac_Hex::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eSumFac),
                                    BwdTrans_SumFac_Hex::create, "BwdTrans_SumFac_Hex");
        
        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Hex : public Operator
        {
        public:
            IProductWRTBase_SumFac_Hex(StdRegions::StdExpansionSharedPtr pExp,
                                       vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                       CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nquad2  (pExp->GetNumPoints(2)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_nmodes2 (pExp->GetBasisNumModes(2)),
                  m_colldir0(pExp->GetBasis(0)->Collocation()),
                  m_colldir1(pExp->GetBasis(1)->Collocation()),
                  m_colldir2(pExp->GetBasis(2)->Collocation())
            {
                m_jac = GeomData->GetJacWithStdWeights(pExp,pGeom);

                m_base0 = GeomData->GetBase(0,pExp);
                m_base1 = GeomData->GetBase(1,pExp);
                m_base2 = GeomData->GetBase(2,pExp);
                m_wspSize = 3*m_numElmt*(max(m_nquad0*m_nquad1*m_nquad2,m_nmodes0*m_nmodes1*m_nmodes2));
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD, NekDouble> &output,
                                    Array<OneD, NekDouble> &output1,
                                    Array<OneD, NekDouble> &output2,
                                    Array<OneD, NekDouble> &wsp)
            {
                int totmodes  = m_nmodes0*m_nmodes1*m_nmodes2;
                int totpoints = m_nquad0 *m_nquad1 *m_nquad2;

                
                if(m_colldir0 && m_colldir1 && m_colldir2)
                {
                
                    Vmath::Vmul(m_numElmt*totpoints,m_jac,1,input,1,output,1);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                    
                    Vmath::Vmul(m_numElmt*totpoints,m_jac,1,input,1,wsp,1);

                    // Assign second half of workspace for 2nd DGEMM operation.
                    Array<OneD, NekDouble> wsp1 = wsp  + totpoints*m_numElmt; 
                    if(m_numElmt < m_nmodes0) // note sure what criterion we should use to swap around these strategies
                    {
                        Array<OneD, NekDouble> wsp2 = wsp1 + m_nmodes0*m_nquad1*m_nquad2;

                        //loop over elements 
                        for(int n = 0; n < m_numElmt; ++n)
                        {
                            if(m_colldir0)
                            {

                                for(int i = 0; i < m_nmodes0; ++i)
                                {
                                    Vmath::Vcopy(m_nquad1*m_nquad2,&wsp[n*totpoints] + i,m_nquad0,
                                                 wsp1.get()+m_nquad1*m_nquad2*i,1);
                                }
                            }
                            else
                            {
                                Blas::Dgemm('T', 'N', m_nquad1*m_nquad2, m_nmodes0, m_nquad0,
                                            1.0, &wsp[n*totpoints],  m_nquad0,
                                            m_base0.get(), m_nquad0,
                                            0.0, wsp1.get(),  m_nquad1*m_nquad2);
                            }

                            
                            if(m_colldir1)
                            {
                                // reshuffle data for next operation.
                                for(int i = 0; i < m_nmodes1; ++i)
                                {
                                    Vmath::Vcopy(m_nquad2*m_nmodes0,wsp1.get()+i,m_nquad1,
                                                wsp2.get()+m_nquad2*m_nmodes0*i,1);
                                }
                            }
                            else
                            {
                                Blas::Dgemm('T', 'N', m_nquad2*m_nmodes0, m_nmodes1, m_nquad1,
                                            1.0, wsp1.get(),  m_nquad1,
                                            m_base1.get(), m_nquad1,
                                            0.0, wsp2.get(),  m_nquad2*m_nmodes0);
                            }
                            
                            if(m_colldir2)
                            {
                                // reshuffle data for next operation.
                                for(int i = 0; i < m_nmodes2; ++i)
                                {
                                    Vmath::Vcopy(m_nmodes0*m_nmodes1,wsp2.get()+i,m_nquad2,
                                                &output[n*totmodes]+m_nmodes0*m_nmodes1*i,1);
                                }
                            }
                            else
                            {
                                Blas::Dgemm('T', 'N', m_nmodes0*m_nmodes1, m_nmodes2, m_nquad2,
                                        1.0, wsp2.get(),  m_nquad2,
                                        m_base2.get(), m_nquad2,
                                        0.0, &output[n*totmodes], m_nmodes0*m_nmodes1);
                            }
                        }
                    }
                    else
                    {
                        Array<OneD, NekDouble> wsp2 = wsp1 + m_numElmt*(max(totpoints,totmodes));
                        
                        if(m_colldir0)
                        {
                            for(int i = 0; i < m_nquad0; ++i)
                            {
                                Vmath::Vcopy(m_nquad1*m_nquad2*m_numElmt,&wsp[i],m_nquad0,
                                            &wsp1[i*m_nquad1*m_nquad2*m_numElmt],1);
                            }
                        }
                        else
                        {
                            // large degmm but copy at end. 
                            Blas::Dgemm('T','N', m_nquad1*m_nquad2*m_numElmt, m_nmodes0, m_nquad0,
                                        1.0, &wsp[0],  m_nquad0,  m_base0.get(),   m_nquad0,
                                        0.0, &wsp1[0], m_nquad1*m_nquad2*m_numElmt);
                        }
                        
                        if(m_colldir1)
                        {
                            for(int i = 0; i < m_nquad1; ++i)
                            {
                                Vmath::Vcopy(m_nquad2*m_numElmt*m_nmodes0,&wsp1[i],m_nquad1,
                                            &wsp2[i*m_nquad2*m_numElmt*m_nmodes0],1);
                            }
                        }
                        else
                        {
                            Blas::Dgemm('T','N', m_nquad2*m_numElmt*m_nmodes0,  m_nmodes1, m_nquad1,
                                        1.0, &wsp1[0],   m_nquad1, m_base1.get(),   m_nquad1,
                                        0.0, &wsp2[0],  m_nquad2*m_numElmt*m_nmodes0);
                        }
                        

                        if(m_numElmt > 1)
                        {
                            if(m_colldir2)
                            {
                                for(int i = 0; i < m_nquad2; ++i)
                                {
                                    Vmath::Vcopy(m_numElmt*m_nmodes0*m_nmodes1,&wsp2[i],m_nquad2,
                                                &wsp1[i*m_numElmt*m_nmodes0*m_nmodes1],1);
                                }
                            }
                            else
                            {
                                
                                Blas::Dgemm('T','N', m_numElmt*m_nmodes0*m_nmodes1, m_nmodes2, m_nquad2,
                                            1.0, &wsp2[0],  m_nquad2,  m_base2.get(),   m_nquad2,
                                            0.0, &wsp1[0],  m_numElmt*m_nmodes0*m_nmodes1);
                            }       
                            
                            for(int i = 0; i < totmodes; ++i)
                            {
                                Vmath::Vcopy(m_numElmt,&wsp1[i*m_numElmt],1,&output[i],totmodes);
                            }
                            
                        }
                        else
                        {
                            if(m_colldir2)
                            {
                                for(int i = 0; i < m_nquad2; ++i)
                                {
                                    Vmath::Vcopy(m_nmodes0*m_nmodes1,&wsp2[i],m_nquad2,
                                                &output[i*m_nmodes0*m_nmodes1],1);
                                }
                            }
                            else
                            {
                                Blas::Dgemm('T','N', m_numElmt*m_nmodes0*m_nmodes1, m_nmodes2, m_nquad2,
                                            1.0, &wsp2[0],  m_nquad2,  m_base2.get(),   m_nquad2,
                                            0.0, &output[0],  m_numElmt*m_nmodes0*m_nmodes1);
                            }
                        }
                    }
                } 
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
            RegisterCreatorFunction(OperatorKey(LibUtilities::eHexahedron, eIProductWRTBase, eSumFac),
                                    IProductWRTBase_SumFac_Hex::create, "IProductWRTBase_SumFac_Hex");

        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */
        
        class PhysDeriv_SumFac : public Operator
        {
        public:
            PhysDeriv_SumFac(StdRegions::StdExpansionSharedPtr pExp,
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
                out[0] = output0;  out[1] = output1;    out[2] = output2;

                for(int i = 0; i < m_dim; ++i)
                {
                    Diff[i] = wsp + i*ntot;
                }

                Blas::Dgemm('N','N', m_nquad0,m_nquad1*m_nquad2*m_numElmt,
                            m_nquad0,1.0, m_Deriv0,m_nquad0,&input[0],
                            m_nquad0,0.0,&Diff[0][0],m_nquad0);
                
                for(int  i = 0; i < m_numElmt; ++i)
                {
                    for (int j = 0; j < m_nquad2; ++j)
                    {
                        Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1,
                                    1.0, &input[i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0, m_Deriv1, m_nquad1, 0.0, 
                                    &Diff[1][i*nPhys+j*m_nquad0*m_nquad1],
                                    m_nquad0);
                    }
                }

                for(int  i = 0; i < m_numElmt; ++i)
                {
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
            
            OPERATOR_CREATE(PhysDeriv_SumFac)

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

        OperatorKey PhysDeriv_SumFac::m_typeArr[] =
        {
            GetOperatorFactory().RegisterCreatorFunction(
                OperatorKey(LibUtilities::eHexahedron, ePhysDeriv, eSumFac),
                PhysDeriv_SumFac::create, "PhysDeriv_SumFac_Hex")
        };
    }
}
