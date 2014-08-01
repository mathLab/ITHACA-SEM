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

    }
}
