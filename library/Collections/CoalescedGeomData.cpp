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

        if(IsDeformed(pCollExp))
        {
            Array<OneD, NekDouble> newjac(npts*nElmts);
            
            //copy Jacobians into a continuous list and set new chatched value
            int cnt = 0;
            for(int i = 0; i < nElmts; ++i)
            {
                const StdRegions::StdExpansion * sep = &(*pCollExp[i]);
                const LocalRegions::Expansion
                    *lep = dynamic_cast<const LocalRegions::Expansion*>( sep );
                
                const Array<OneD, const NekDouble>
                    jac = lep->GetMetricInfo()->GetJac( ptsKeys );
                
                Vmath::Vcopy(npts, &jac[0], 1, &newjac[cnt], 1);
                
                cnt += npts;
            }
            
            m_oneDGeomData[eJac] = newjac;
        }
        else
        {
            Array<OneD, NekDouble> newjac(nElmts);
            //copy Jacobians into a continuous list 
            for(int i = 0; i < nElmts; ++i)
            {
                const StdRegions::StdExpansion * sep = &(*pCollExp[i]);
                const LocalRegions::Expansion  * lep =
                    dynamic_cast<const LocalRegions::Expansion*>( sep );
                
                const Array<OneD, const NekDouble> jac =
                    lep->GetMetricInfo()->GetJac( ptsKeys );
                
                newjac[i] = jac[0]; 
                
            }
            m_oneDGeomData[eJac] = newjac;
        }
    }

    return m_oneDGeomData[eJac];
}

const std::shared_ptr<VecVec_t> CoalescedGeomData::GetJacInterLeave(
                         vector<StdRegions::StdExpansionSharedPtr> &pCollExp,
                         int nElmt)
{

    if(m_oneDGeomDataInterLeave.count(eJac) == 0)
    {
        const Array<OneD, const NekDouble> jac = GetJac(pCollExp);
        int jacsize = jac.size();

        ASSERTL1(nElmt % vec_t::width == 0,
                 "Number of elements not divisible by vector "
                 "width, padding not yet implemented.");
        int  nBlocks = nElmt / vec_t::width;
        
        VecVec_t newjac; 

        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();
            
        // set up Cached Jacobians to be continuous
        int nq = 1;
        for (int i = 0; i < ptsKeys.size(); ++i)
        {
            nq  *= ptsKeys[i].GetNumPoints();
        }
        
        if(IsDeformed(pCollExp))
        {

            newjac.resize(nBlocks*nq);

            alignas(vec_t::alignment) NekDouble tmp[vec_t::width];

            for (size_t block = 0; block < nBlocks; ++block)
            {
                size_t nblock_width = block*nq*vec_t::width;
                for(size_t q = 0; q < nq; q++)
                {
                    for (int j = 0; j < vec_t::width; ++j)
                    {
                        if(nblock_width+ nq*j + q < jacsize)
                        {
                            tmp[j] = jac[nblock_width + nq*j + q];
                        }
                        else
                        {
                            tmp[j] = 0.0; 
                        }
                    }

                    //Order is [block][quadpt]
                    newjac[block*nq + q].load(&tmp[0]);
                }
            }
        }
        else
        {
            newjac.resize(nBlocks);

            alignas(vec_t::alignment) NekDouble tmp[vec_t::width];
            for (size_t i = 0; i < nBlocks; ++i)
            {
                for (int j = 0; j < vec_t::width; ++j)
                {
                    if(vec_t::width*i+j < jacsize)
                    {
                        tmp[j] = jac[vec_t::width*i+j];
                    }
                    else
                    {
                        tmp[j] = 0.0;
                    }
                }

                newjac[i].load(&tmp[0]);
            }
        }

        m_oneDGeomDataInterLeave[eJac] =
            MemoryManager<VecVec_t>::AllocateSharedPtr(newjac);
    }

    return m_oneDGeomDataInterLeave[eJac];
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
            const LocalRegions::Expansion
                *lep = dynamic_cast<const LocalRegions::Expansion*>( sep );

            const Array<TwoD, const NekDouble>
                Dfac = lep->GetMetricInfo()->GetDerivFactors( ptsKeys );

            if(IsDeformed(pCollExp))
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
                    newDFac[j][i] = Dfac[j][0];
                }
            }
            cnt += npts;
        }

        m_twoDGeomData[eDerivFactors] = newDFac;
    }

    return m_twoDGeomData[eDerivFactors];
}

const std::shared_ptr<VecVec_t> CoalescedGeomData::GetDerivFactorsInterLeave(
                         vector<StdRegions::StdExpansionSharedPtr> &pCollExp,
                         int nElmt)
{
    if(m_twoDGeomDataInterLeave.count(eDerivFactors) == 0)
    {
        ASSERTL1(nElmt % vec_t::width == 0,
                 "Number of elements not divisible by vector "
                 "width, padding not yet implemented.");

        int  nBlocks = nElmt / vec_t::width;
        
        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();
        const int coordim = pCollExp[0]->GetCoordim();
        int dim = ptsKeys.size();


        unsigned int n_df = coordim*dim;
        alignas(vec_t::alignment) NekDouble vec[vec_t::width];

        const Array<TwoD, const NekDouble> df = GetDerivFactors(pCollExp);
        int dfsize = df.GetColumns();

        VecVec_t newdf; 

        int nq = 1;
        for (int i = 0; i < dim; ++i)
        {
            nq  *= ptsKeys[i].GetNumPoints();
        }
            
        if(IsDeformed(pCollExp))
        {
            newdf.resize(nBlocks * n_df *nq);
            auto *df_ptr = &newdf[0];
            for (int e = 0; e < nBlocks; ++e)
            {
                for (int q = 0; q < nq; q++)
                {
                    for (int dir = 0; dir < n_df; ++dir, ++df_ptr)
                    {
                        for (int j = 0; j < vec_t::width; ++j)
                        {
                            // manage padding 
                            if((vec_t::width*e + j)*nq + q < dfsize)
                            {
                                vec[j] = df[dir][(vec_t::width*e + j)*nq + q];
                            }
                            else
                            {
                                vec[j] = 0.0;
                            }
                        }
                        (*df_ptr).load(&vec[0]);
                    }
                }
            }
        }
        else
        {
            newdf.resize(nBlocks * n_df);
            for (int e = 0; e < nBlocks; ++e)
            {
                for (int dir = 0; dir < n_df; ++dir)
                {
                    for (int j = 0; j < vec_t::width; ++j)
                    {
                        // padding
                        if(vec_t::width*e + j < dfsize)
                        {
                            vec[j] = df[dir][vec_t::width*e + j];
                        }
                        else
                        {
                            vec[j] = 0.0;
                        }
                    }
                    // Must have all vec_t::width elemnts aligned to do a load.
                    newdf[e*n_df + dir].load(&vec[0]);
                }
            }
        }
        
        m_twoDGeomDataInterLeave[eDerivFactors] =
            MemoryManager<VecVec_t>::AllocateSharedPtr(newdf);
    }
    
    return m_twoDGeomDataInterLeave[eDerivFactors];
}


bool CoalescedGeomData::IsDeformed(
    vector<StdRegions::StdExpansionSharedPtr> &pCollExp)
{
    if (!m_isDeformedSet)
    {
        LibUtilities::PointsKeyVector ptsKeys = pCollExp[0]->GetPointsKeys();
        const StdRegions::StdExpansion * sep = &(*pCollExp[0]);
        const LocalRegions::Expansion  * lep =
            dynamic_cast<const LocalRegions::Expansion*>( sep );

        const Array<OneD, const NekDouble> jac =
            lep->GetMetricInfo()->GetJac(ptsKeys);

        m_deformed = lep->GetMetricInfo()->GetGtype() ==
            SpatialDomains::eDeformed;
    }

    return m_deformed;
}

}
}
