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

#include <SpatialDomains/GeomFactors.h>
#include <Collections/CoalescedGeomData.h>

#include <LocalRegions/Expansion.h>

using namespace std;

namespace Nektar {
namespace Collections {

CoalescedGeomData::CoalescedGeomData(void)
{
}

CoalescedGeomData::~CoalescedGeomData(void)
{
}

const Array<OneD, const NekDouble> &CoalescedGeomData::GetJac(
        vector<StdRegions::StdExpansionSharedPtr> &pCollExp)
{

    if(m_oneDGeomData.count(eJac) == 0)
    {

        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();
        int nElmts = pCollExp.size();

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
            const StdRegions::StdExpansion * sep = &(*pCollExp[i]);
            const LocalRegions::Expansion  * lep = dynamic_cast<const LocalRegions::Expansion*>( sep );
          
            const Array<OneD, const NekDouble> jac = lep->GetMetricInfo()->GetJac( ptsKeys );

            if( lep->GetMetricInfo()->GetGtype() == SpatialDomains::eDeformed )
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



const Array<OneD, const NekDouble> &CoalescedGeomData::GetJacWithStdWeights(
        vector<StdRegions::StdExpansionSharedPtr> &pCollExp)
{
    if(m_oneDGeomData.count(eJacWithStdWeights) == 0)
    {
        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();
        int nElmts = pCollExp.size();

        // set up Cached Jacobians to be continuous
        int npts = 1;
        for (int i = 0; i < ptsKeys.size(); ++i)
        {
            npts   *= ptsKeys[i].GetNumPoints();
        }


        Array<OneD, NekDouble> newjac(npts*nElmts), tmp;

        //copy Jacobians into a continuous list and set new chatched value
        int cnt = 0;
        for(int i = 0; i < nElmts; ++i)
        {
            const StdRegions::StdExpansion * sep = &(*pCollExp[i]);
            const LocalRegions::Expansion  * lep = dynamic_cast<const LocalRegions::Expansion*>( sep );

            const Array<OneD, const NekDouble> jac = lep->GetMetricInfo()->GetJac(ptsKeys);

            if( lep->GetMetricInfo()->GetGtype() == SpatialDomains::eDeformed )
            {
                Vmath::Vcopy(npts, &jac[0], 1, &newjac[cnt], 1);
            }
            else
            {
                Vmath::Fill(npts, jac[0], &newjac[cnt], 1);
            }

            pCollExp[0]->MultiplyByStdQuadratureMetric(newjac + cnt,
                                                       tmp = newjac + cnt);
            cnt += npts;
        }

        m_oneDGeomData[eJacWithStdWeights] = newjac;
    }

    return m_oneDGeomData[eJacWithStdWeights];
}


const Array<TwoD, const NekDouble> &CoalescedGeomData::GetDerivFactors(
        vector<StdRegions::StdExpansionSharedPtr> &pCollExp)
{
    if(m_twoDGeomData.count(eDerivFactors) == 0)
    {
        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();

        int nElmts = pCollExp.size();
        const int coordim = pCollExp[0]->GetCoordim();
        int dim = ptsKeys.size();

        // set up Cached Jacobians to be continuous
        int npts = 1;
        for (int i = 0; i < dim; ++i)
        {
            npts   *= ptsKeys[i].GetNumPoints();
        }


        Array<TwoD, NekDouble> newDFac(dim*coordim,npts*nElmts);

        //copy Jacobians into a continuous list and set new chatched value
        int cnt = 0;
        for(int i = 0; i < nElmts; ++i)
        {
            const StdRegions::StdExpansion * sep = &(*pCollExp[i]);
            const LocalRegions::Expansion  * lep = dynamic_cast<const LocalRegions::Expansion*>( sep );

            const Array<TwoD, const NekDouble> Dfac = lep->GetMetricInfo()->GetDerivFactors( ptsKeys );

            if( lep->GetMetricInfo()->GetGtype() == SpatialDomains::eDeformed)
            {
                for (int j = 0; j < dim*coordim; ++j)
                {
                    Vmath::Vcopy(npts, &Dfac[j][0], 1, &newDFac[j][cnt], 1);
                }
            }
            else
            {
                for (int j = 0; j < dim*coordim; ++j)
                {
                    Vmath::Fill(npts, Dfac[j][0], &newDFac[j][cnt], 1);
                }
            }
            cnt += npts;
        }

        m_twoDGeomData[eDerivFactors] = newDFac;
    }

    return m_twoDGeomData[eDerivFactors];
}

}
}
