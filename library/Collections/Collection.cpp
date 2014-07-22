///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.cpp
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
// and/or sell copies of the Software, and to permit persons to whom the
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>

namespace Nektar {

    namespace Collections {


        CoalescedGeomData::CoalescedGeomData(void)
        {
        }

        CoalescedGeomData::~CoalescedGeomData(void)
        {
        }

        const Array<OneD, const NekDouble> &CoalescedGeomData::GetJac(const LibUtilities::PointsKeyVector &ptsKeys,
                                                               vector<SpatialDomains::GeometrySharedPtr> &pGeom)
        {
            if(m_oneDGeomData.count(eJac) == 0)
            {

                int nElmts = pGeom.size();
                
                // set up Cached Jacobians to be continuous 
                int npts = 1;
                for (int i = 0; i < ptsKeys.size(); ++i)
                {
                    npts   *= ptsKeys[i].GetNumPoints();
                }
                

                Array<OneD, NekDouble> newjac(npts*nElmts);
            
                //copy Jacobians into a continuous list and set new chatched value
                int cnt = 0;
                for(int i = 0; i < nElmts; ++i)
                {
                    const Array<OneD, const NekDouble> jac= pGeom[i]->GetGeomFactors()->GetJac(ptsKeys);
                    
                    if (pGeom[i]->GetGeomFactors()->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vcopy(npts, &jac[0], 1, &newjac[cnt], 1);
                    }
                    else
                    {
                        Vmath::Fill(npts, jac[0], &newjac[cnt], 1);
                    }

                    cnt += npts;
                }
                
                m_oneDGeomData[eJac] = newjac; 
            }
            
            return m_oneDGeomData[eJac];
        }
        
        
        Collection::Collection(StdRegions::StdExpansionSharedPtr pExp,
                   vector<SpatialDomains::GeometrySharedPtr> pGeom)
            : m_stdExp(pExp), m_geom(pGeom)
        {
            
            OperatorKey bwdStdMat    (pExp->DetShapeType(), eBwdTrans, eStdMat);
            OperatorKey bwdIterPerExp(pExp->DetShapeType(), eBwdTrans, eIterPerExp);
            OperatorKey bwdSumFac    (pExp->DetShapeType(), eBwdTrans, eSumFac);
            
            OperatorKey iproductWRTBaseStdMat(
                pExp->DetShapeType(), eIProductWRTBase, eStdMat);
            OperatorKey iproductWRTBaseIterPerExp(
                pExp->DetShapeType(), eIProductWRTBase, eIterPerExp);

            OperatorKey derivSumFac(
                pExp->DetShapeType(), ePhysDeriv, eSumFac);
            OperatorKey derivIterPerExp(
                pExp->DetShapeType(), ePhysDeriv, eIterPerExp);

            m_geomData = MemoryManager<CoalescedGeomData>::AllocateSharedPtr();
            
            m_ops[eBwdTrans]        = GetOperatorFactory().CreateInstance(
                                                          bwdSumFac, pExp, pGeom, m_geomData);
            //m_ops[eBwdTrans]        = GetOperatorFactory().CreateInstance(
            //                                            bwdStdMat, pExp, pGeom, m_geomData);
            m_ops[eIProductWRTBase] = GetOperatorFactory().CreateInstance(
                                                          iproductWRTBaseIterPerExp, pExp, pGeom,m_geomData);

            //m_ops[eBwdTrans] = GetOperatorFactory().CreateInstance(
            //bwdStdMat, pExp, pGeom);
            //m_ops[ePhysDeriv] = GetOperatorFactory().CreateInstance(
            //    derivIterPerExp, pExp, pGeom);
        }
    }
}
    
