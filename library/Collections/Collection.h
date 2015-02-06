///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.h
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

#ifndef NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H
#define NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H

#include <Collections/Operator.h>

#include <boost/unordered_map.hpp>

namespace Nektar
{
namespace Collections
{

enum GeomData
{
    eJac,
    eJacWithStdWeights,
    eDerivFactors
};

class CoalescedGeomData
{
public:
    CoalescedGeomData(void);

    virtual ~CoalescedGeomData(void);

    const Array<OneD, const NekDouble> &GetJac(
            vector<StdRegions::StdExpansionSharedPtr> &pColLExp);

    const Array<OneD, const NekDouble> &GetJacWithStdWeights(
            vector<StdRegions::StdExpansionSharedPtr> &pColLExp);

    const Array<TwoD, const NekDouble> &GetDerivFactors(
            vector<StdRegions::StdExpansionSharedPtr> &pColLExp);

private:
    map<GeomData,Array<OneD, NekDouble> > m_oneDGeomData;
    map<GeomData,Array<TwoD, NekDouble> > m_twoDGeomData;
};

typedef boost::shared_ptr<CoalescedGeomData>   CoalescedGeomDataSharedPtr;

static CoalescedGeomDataSharedPtr GeomDataNull;


/**
 * @brief Collection
 */
class Collection
{
public:

    Collection(vector<StdRegions::StdExpansionSharedPtr> pColLExp,
               OperatorImpMap &impTypes);

    void ApplyOperator(const OperatorType                 &op,
                       const Array<OneD, const NekDouble> &inarray,
                             Array<OneD,       NekDouble> &output)
    {
        Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
        (*m_ops[op])(inarray, output, NullNekDouble1DArray,
                     NullNekDouble1DArray, wsp);
    }


    void ApplyOperator(const OperatorType                 &op,
                       const Array<OneD, const NekDouble> &inarray,
                             Array<OneD,       NekDouble> &output0,
                             Array<OneD,       NekDouble> &output1)
    {
        Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
        (*m_ops[op])(inarray, output0, output1, NullNekDouble1DArray, wsp);
    }

    void ApplyOperator(const OperatorType                 &op,
                       const Array<OneD, const NekDouble> &inarray,
                             Array<OneD,       NekDouble> &output0,
                             Array<OneD,       NekDouble> &output1,
                             Array<OneD,       NekDouble> &output2)
    {
        Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
        (*m_ops[op])(inarray, output0, output1, output2, wsp);
    }

protected:
    StdRegions::StdExpansionSharedPtr                     m_stdExp;
    vector<SpatialDomains::GeometrySharedPtr>             m_geom;
    boost::unordered_map<OperatorType, OperatorSharedPtr> m_ops;
    CoalescedGeomDataSharedPtr                            m_geomData;

};

typedef std::vector<Collection> CollectionVector;
typedef boost::shared_ptr<CollectionVector> CollectionVectorSharedPtr;

}
}

#endif
