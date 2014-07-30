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
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                
                if(m_colldir0 )
                {
                    Vmath::Vcopy(m_numElmt*m_nmodes0,input.get(),1,output.get(),1);
                }
                else
                { 
                    Array<OneD, const NekDouble> base0  = m_stdExp->GetBasis(0)->GetBdata();

                    // out = B0*in; 
                    Blas::Dgemm('N','N', m_nquad0,m_numElmt,m_nmodes0,1.0, base0.get(),
                                m_nquad0, &input[0], m_nmodes0,0.0,&output[0], m_nquad0);
                }
            }
            
            OPERATOR_CREATE(BwdTrans_SumFac_Seg)
            
            protected:
            const int  m_nquad0;
            const int  m_nmodes0;
            const bool m_colldir0;
        };
        
        OperatorKey BwdTrans_SumFac_Seg::m_type = GetOperatorFactory().
            RegisterCreatorFunction(OperatorKey(LibUtilities::eSegment, eBwdTrans, 
                                                eSumFac),
                                    BwdTrans_SumFac_Seg::create, "BwdTrans_SumFac_Seg");



    }
}
