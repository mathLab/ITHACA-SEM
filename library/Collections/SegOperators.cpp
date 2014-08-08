///////////////////////////////////////////////////////////////////////////////
//
// File: SegOperators.cpp
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
        class BwdTrans_SumFac_Seg : public Operator
        {
        public:
            BwdTrans_SumFac_Seg(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_colldir0(pExp->GetBasis(0)->Collocation())
            {
                m_wspSize = 0;
                m_base0 = GeomData->GetBase(0,pExp);
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                if(m_colldir0)
                {
                    Vmath::Vcopy(m_numElmt*m_nmodes0,input.get(),1,output.get(),1);
                }
                else
                { 
                    // out = B0*in; 
                    Blas::Dgemm('N','N', m_nquad0,m_numElmt,m_nmodes0,1.0, m_base0.get(),
                                m_nquad0, &input[0], m_nmodes0,0.0,&output[0], m_nquad0);
                }
            }
            
            OPERATOR_CREATE(BwdTrans_SumFac_Seg)
            
            protected:
            const int  m_nquad0;
            const int  m_nmodes0;
            const bool m_colldir0;
            Array<OneD, const NekDouble> m_base0;
        };
        
        OperatorKey BwdTrans_SumFac_Seg::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eSegment, eBwdTrans, 
                                                eSumFac),
                                    BwdTrans_SumFac_Seg::create, "BwdTrans_SumFac_Seg");



        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */       
        class IProductWRTBase_SumFac_Seg : public Operator
        {
        public:
            IProductWRTBase_SumFac_Seg(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator  (pExp, pGeom, GeomData),
                  m_nquad0  (pExp->GetNumPoints(0)),
                  m_nmodes0 (pExp->GetBasisNumModes(0)),
                  m_colldir0(pExp->GetBasis(0)->Collocation())
            {
                m_wspSize = m_numElmt*m_nquad0;
                m_jac = GeomData->GetJacWithStdWeights(pExp,pGeom);
                m_base0 = GeomData->GetBase(0,pExp);
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                if(m_colldir0)
                {
                    Vmath::Vmul(m_numElmt*m_nquad0,m_jac,1,input,1,output,1);
                }
                else
                { 
                    Vmath::Vmul(m_numElmt*m_nquad0,m_jac,1,input,1,wsp,1);

                    // out = B0*in; 
                    Blas::Dgemm('T','N', m_nmodes0,m_numElmt,m_nquad0,1.0, m_base0.get(), m_nquad0,
                                &wsp[0], m_nquad0, 0.0,&output[0], m_nmodes0);
                }
            }
            
            OPERATOR_CREATE(IProductWRTBase_SumFac_Seg)
            
            protected:
            const int  m_nquad0;
            const int  m_nmodes0;
            const bool m_colldir0;
            Array<OneD, const NekDouble> m_jac;
            Array<OneD, const NekDouble> m_base0;
        };
        
        OperatorKey IProductWRTBase_SumFac_Seg::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eSegment, eIProductWRTBase, 
                                                eSumFac),
                                    IProductWRTBase_SumFac_Seg::create, "IProductWRTBase_SumFac_Seg");


        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */       
        class PhysDeriv_SumFac_Seg : public Operator
        {
        public:
            PhysDeriv_SumFac_Seg(StdRegions::StdExpansionSharedPtr pExp,
                                  vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                  CoalescedGeomDataSharedPtr GeomData)
                : Operator (pExp, pGeom, GeomData),
                  m_nquad0 (pExp->GetNumPoints(0))
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_coordim = pExp->GetCoordim();

                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);

                m_Deriv0 = &((pExp->GetBasis(0)->GetD())->GetPtr())[0];
                m_wspSize = m_nquad0*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                const int nqcol   = m_nquad0*m_numElmt;
                
                ASSERTL1(wsp.num_elements() == m_wspSize,
                         "Incorrect workspace size");
                ASSERTL1(input.num_elements() >= nqcol,
                         "Incorrect input size");
                
                Array<OneD, NekDouble> diff0(nqcol, wsp);
                
                Blas::Dgemm('N', 'N', m_nquad0, m_numElmt, 
                            m_nquad0, 1.0, m_Deriv0, m_nquad0, 
                            input.get(), m_nquad0, 0.0,
                            diff0.get(), m_nquad0);
                
                Vmath::Vmul  (nqcol, m_derivFac[0], 1, diff0, 1, output0, 1);
                
                if (m_coordim == 2)
                {
                    Vmath::Vmul  (nqcol, m_derivFac[1], 1, diff0, 1, output1, 1);
                }
                else if (m_coordim == 3)
                {
                    Vmath::Vmul  (nqcol, m_derivFac[1], 1, diff0, 1, output1, 1);
                    Vmath::Vmul  (nqcol, m_derivFac[2], 1, diff0, 1, output2, 1);
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_SumFac_Seg)
            
            protected:
            int m_coordim;
            const int m_nquad0;
            Array<TwoD, const NekDouble> m_derivFac;
            NekDouble *m_Deriv0;
        };
        
        OperatorKey PhysDeriv_SumFac_Seg::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eSegment,
                                                ePhysDeriv, eSumFac),
                     PhysDeriv_SumFac_Seg::create, "PhysDeriv_SumFac_Seg");




    }
}
