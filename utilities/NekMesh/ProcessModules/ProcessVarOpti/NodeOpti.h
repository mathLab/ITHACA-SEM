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
//  License for the specific language governing rights and limitations under
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

#include <LibUtilities/BasicUtils/Thread.h>

#include "ProcessVarOpti.h"

namespace Nektar
{
namespace Utilities
{

class NodeOptiJob;

class NodeOpti
{
public:
    NodeOpti(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
             ResidualSharedPtr r,
             std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
             optiType o)
        : m_node(n), m_res(r), m_derivUtils(d), m_opti(o)
    {
        // filter element types within d vector
        for (int i = 0; i < e.size(); i++)
        {
            m_data[e[i]->GetEl()->GetShapeType()].push_back(e[i]);
        }
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
    boost::mutex mtx;
    std::map<LibUtilities::ShapeType, std::vector<ElUtilSharedPtr> > m_data;
    Array<OneD, NekDouble> m_grad;

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

typedef boost::shared_ptr<NodeOpti> NodeOptiSharedPtr;
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
        : NodeOpti(n, e, r, d, o)
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
        : NodeOpti(n, e, r, d, o)
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
