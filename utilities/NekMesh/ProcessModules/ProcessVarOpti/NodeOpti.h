////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI
#define UTILITIES_NEKMESH_NODEOPTI

#include <ostream>
#include <mutex>

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/BasicUtils/Thread.h>

#include "ProcessVarOpti.h"

namespace Nektar
{
namespace Utilities
{

class NodeOptiJob;

class NodeOpti
{
    // Typedef for derivative storage, we use boost::multi_array so we can pass
    // this to functions easily
    typedef boost::multi_array<NekDouble, 4> DerivArray;

public:
    NodeOpti(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
             ResidualSharedPtr r,
             std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
             optiType o, int dim)
        : m_node(n), m_res(r), m_derivUtils(d), m_opti(o)
    {
        // filter element types within d vector
        for (int i = 0; i < e.size(); i++)
        {
            m_data[e[i]->GetEl()->GetShapeType()].push_back(e[i]);
        }

        // Set up storage for GetFunctional to avoid reallocation on each call.
        size_t storageCount = 0;

        // Count total storage needed.
        for (auto &typeIt : m_data)
        {
            const int pts    = m_derivUtils[typeIt.first]->pts;
            const int nElmt  = typeIt.second.size();

            storageCount = std::max(storageCount,
                                    dim * m_derivUtils[typeIt.first]->ptsStd *
                                    typeIt.second.size());

            m_derivs.insert(std::make_pair(
                                typeIt.first,
                                DerivArray(boost::extents[dim][nElmt][dim][pts])));
        }

        m_tmpStore.resize(storageCount);
    }

    virtual ~NodeOpti(){};

    void CalcMinJac();

    virtual void Optimise() = 0;
    NodeOptiJob *GetJob();

    template <int DIM>
    NekDouble GetFunctional(NekDouble &minJacNew, bool gradient = true);

    template <int DIM> void MinEigen(NekDouble &val);

protected:
    NodeSharedPtr m_node;
    std::mutex mtx;
    std::map<LibUtilities::ShapeType, std::vector<ElUtilSharedPtr> > m_data;
    std::vector<NekDouble> m_grad;
    std::vector<NekDouble> m_tmpStore;
    std::unordered_map<LibUtilities::ShapeType, DerivArray, EnumHash> m_derivs;


    template <int DIM> int IsIndefinite();

    NekDouble m_minJac;
    ResidualSharedPtr m_res;
    std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> m_derivUtils;
    optiType m_opti;

    static NekDouble c1()
    {
        return 1e-3;
    }
    static NekDouble gradTol()
    {
        return 1e-8;
    }
    static NekDouble alphaTol()
    {
        return 1e-8;
    }
};

typedef std::shared_ptr<NodeOpti> NodeOptiSharedPtr;
typedef LibUtilities::NekFactory<
    int, NodeOpti, NodeSharedPtr, std::vector<ElUtilSharedPtr>,
    ResidualSharedPtr, std::map<LibUtilities::ShapeType, DerivUtilSharedPtr>,
    optiType>
    NodeOptiFactory;

NodeOptiFactory &GetNodeOptiFactory();

class NodeOpti3D3D : public NodeOpti // 1D optimsation in 3D space
{
public:
    NodeOpti3D3D(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
                 ResidualSharedPtr r,
                 std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
                 optiType o)
        : NodeOpti(n, e, r, d, o, 3)
    {
    }

    ~NodeOpti3D3D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n, std::vector<ElUtilSharedPtr> e, ResidualSharedPtr r,
        std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d, optiType o)
    {
        return NodeOptiSharedPtr(new NodeOpti3D3D(n, e, r, d, o));
    }

private:
};

class NodeOpti2D2D : public NodeOpti // 1D optimsation in 3D space
{
public:
    NodeOpti2D2D(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
                 ResidualSharedPtr r,
                 std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
                 optiType o)
        : NodeOpti(n, e, r, d, o, 2)
    {
    }

    ~NodeOpti2D2D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n, std::vector<ElUtilSharedPtr> e, ResidualSharedPtr r,
        std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d, optiType o)
    {
        return NodeOptiSharedPtr(new NodeOpti2D2D(n, e, r, d, o));
    }

private:
};

class NodeOptiJob : public Thread::ThreadJob
{
public:
    NodeOptiJob(NodeOpti *no) : node(no)
    {
    }

    void Run()
    {
        node->Optimise();
    }

private:
    NodeOpti *node;
};
}
}

#endif
