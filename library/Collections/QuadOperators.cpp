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
        class BwdTrans_SumFac_Quad : public Operator
        {
        public:
            BwdTrans_SumFac_Quad(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_colldir0(pExp->GetBasis(0)->Collocation()),
                  m_colldir1(pExp->GetBasis(1)->Collocation())
            {
                m_wspSize = m_nquad0*m_nmodes1*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                if(m_colldir0 && m_colldir1)
                {
                    Vmath::Vcopy(m_numElmt*m_nmodes0*m_nmodes1,input.get(),1,output.get(),1);
                }
                else if(m_colldir0)
                {
                    Array<OneD, const NekDouble> base1  = m_stdExp->GetBasis(1)->GetBdata();
                    for(int i = 0; i < m_numElmt; ++i)
                    {
                        Blas::Dgemm('N','T', m_nquad0, m_nquad1,m_nmodes1, 1.0, &input[i*m_nquad0*m_nmodes1], m_nquad0, 
                                    base1.get(), m_nquad1, 0.0, &output[i*m_nquad0*m_nquad1], m_nquad0);
                    }
                }
                else if(m_colldir1)
                {
                    Array<OneD, const NekDouble> base0  = m_stdExp->GetBasis(0)->GetBdata();
                    Blas::Dgemm('N','N', m_nquad0,m_nmodes1*m_numElmt,m_nmodes0,1.0, base0.get(),
                                m_nquad0, &input[0], m_nmodes0,0.0,&output[0], m_nquad0);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                    
                    Array<OneD, const NekDouble> base0  = m_stdExp->GetBasis(0)->GetBdata();
                    Array<OneD, const NekDouble> base1  = m_stdExp->GetBasis(1)->GetBdata();

                    // Those two calls correpsond to the operation
                    // out = B0*in*Transpose(B1); 
                    Blas::Dgemm('N','N', m_nquad0,m_nmodes1*m_numElmt,m_nmodes0,1.0, base0.get(),
                                m_nquad0, &input[0], m_nmodes0,0.0,&wsp[0], m_nquad0);
            
                    for(int i = 0; i < m_numElmt; ++i)
                    {
                        Blas::Dgemm('N','T', m_nquad0, m_nquad1,m_nmodes1, 1.0, &wsp[i*m_nquad0*m_nmodes1], m_nquad0, 
                                    base1.get(), m_nquad1, 0.0, &output[i*m_nquad0*m_nquad1], m_nquad0);
                    }
                } 
            }
            
            OPERATOR_CREATE(BwdTrans_SumFac_Quad)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const bool m_colldir0;
            const bool m_colldir1;
        };
        
        OperatorKey BwdTrans_SumFac_Quad::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, 
                                                eSumFac),
                                    BwdTrans_SumFac_Quad::create, "BwdTrans_SumFac_Quad");


        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Quad : public Operator
        {
        public:
            IProductWRTBase_SumFac_Quad(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nquad1  (pExp->GetNumPoints(1)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_nmodes1 (pExp->GetBasisNumModes(1)),
                  m_colldir0(pExp->GetBasis(0)->Collocation()),
                  m_colldir1(pExp->GetBasis(1)->Collocation())
            {
                m_jac     = GeomData->GetJacWithStdWeights(pExp,pGeom);
                m_base0   = GeomData->GetBase(0,pExp);
                m_base1   = GeomData->GetBase(1,pExp);
                m_wspSize = 2*m_numElmt*(max(m_nquad0*m_nquad1,m_nmodes0*m_nmodes1));
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                int totmodes  = m_nmodes0*m_nmodes1;
                int totpoints = m_nquad0 *m_nquad1;
                
                Vmath::Vmul(m_numElmt*totpoints,m_jac,1,input,1,wsp,1);

                if(m_colldir0 && m_colldir1)
                {
                    Vmath::Vcopy(m_numElmt*totmodes,wsp.get(),1,output.get(),1);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements() == m_wspSize, "Incorrect workspace size");
                    
                    Array<OneD, NekDouble> wsp1 = wsp  + max(totpoints,totmodes)*m_numElmt; 

                    if(m_colldir0)
                    {
                        for(int i = 0; i < m_nquad0; ++i)
                        {
                            Vmath::Vcopy(m_nquad1*m_numElmt,&wsp[i],m_nquad0,
                                         &wsp1[i*m_nquad1*m_numElmt],1);
                        }
                    }
                    else
                    {
                        Blas::Dgemm('T','N', m_nquad1*m_numElmt,m_nmodes0,m_nquad0,1.0,
                                &wsp[0],m_nquad0, m_base0.get(), m_nquad0, 
                                    0.0,&wsp1[0], m_nquad1*m_numElmt);
                    }

                    
                    if(m_numElmt > 1)
                    {

                        if(m_colldir1)
                        {
                            for(int i = 0; i < m_nquad1; ++i)
                            {
                                Vmath::Vcopy(m_numElmt*m_nmodes0,&wsp1[i],m_nquad1,
                                             &wsp[i*m_numElmt*m_nmodes0],1);
                            }
                        }
                        else
                        {
                            
                            Blas::Dgemm('T','N', m_numElmt*m_nmodes0,  m_nmodes1, m_nquad1,
                                        1.0, &wsp1[0], m_nquad1, m_base1.get(),   m_nquad1,
                                        0.0, &wsp[0], m_numElmt*m_nmodes0);
                        }
                        
                        for(int i = 0; i < totmodes; ++i)
                        {
                            Vmath::Vcopy(m_numElmt,&wsp[i*m_numElmt],1,&output[i],totmodes);
                        }
                    }
                    else
                    {
                        if(m_colldir1)
                        {
                            for(int i = 0; i < m_nquad1; ++i)
                            {
                                Vmath::Vcopy(m_numElmt*m_nmodes0,&wsp1[i],m_nquad1,
                                             &output[i*m_numElmt*m_nmodes0],1);
                            }
                        }
                        else
                        {
                            Blas::Dgemm('T','N', m_nmodes0,  m_nmodes1, m_nquad1,
                                        1.0, &wsp1[0], m_nquad1, m_base1.get(),   m_nquad1,
                                        0.0, &output[0], m_nmodes0);
                        }
                    }
                } 
            }
            
            OPERATOR_CREATE(IProductWRTBase_SumFac_Quad)
            
            protected:
            const int  m_nquad0;
            const int  m_nquad1;
            const int  m_nmodes0;
            const int  m_nmodes1;
            const bool m_colldir0;
            const bool m_colldir1;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
            Array<OneD, const NekDouble> m_base1;
        };
        
        OperatorKey IProductWRTBase_SumFac_Quad::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eQuadrilateral, eIProductWRTBase, eSumFac),IProductWRTBase_SumFac_Quad::create, "IProductWRTBase_SumFac_Quad");



        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */       
        class PhysDeriv_SumFac_Quad : public Operator
        {
        public:
            PhysDeriv_SumFac_Quad(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator (pExp, pGeom, GeomData),
                  m_nquad0 (pExp->GetNumPoints(0)),
                  m_nquad1 (pExp->GetNumPoints(1))
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_coordim = pExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);

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
                
                Blas::Dgemm('N', 'N', m_nquad0, m_nquad1*m_numElmt, 
                            m_nquad0, 1.0, m_Deriv0, m_nquad0, 
                            input.get(), m_nquad0, 0.0,
                            diff0.get(), m_nquad0);
                
                int cnt = 0;
                for (int i = 0; i < m_numElmt; ++i, cnt += nqtot)
                {
                    Blas::Dgemm('N', 'T', m_nquad0, m_nquad1, m_nquad1, 1.0,
                                input.get() + cnt, m_nquad0, 
                                m_Deriv1, m_nquad1, 0.0,
                                diff1.get() + cnt, m_nquad0);
                }
                
                Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[1], 1, diff1, 1, output0, 1, output0, 1);
                Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output1, 1);
                Vmath::Vvtvp (nqcol, m_derivFac[3], 1, diff1, 1, output1, 1, output1, 1);
                
                if (m_coordim == 3)
                {
                    Vmath::Vmul  (nqcol, m_derivFac[4], 1, diff0, 1, output2, 1);
                    Vmath::Vvtvp (nqcol, m_derivFac[5], 1, diff1, 1, output2, 1, output2, 1);
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Quad)
            
            protected:
            int m_coordim;
            const int m_nquad0;
            const int m_nquad1;
            Array<TwoD, const NekDouble> m_derivFac;
            NekDouble *m_Deriv0;
            NekDouble *m_Deriv1;
        };
        
        OperatorKey PhysDeriv_SumFac_Quad::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eQuadrilateral, ePhysDeriv, 
                                                eSumFac),
                                    PhysDeriv_SumFac_Quad::create, "PhysDeriv_SumFac_Quad");


    }
}
