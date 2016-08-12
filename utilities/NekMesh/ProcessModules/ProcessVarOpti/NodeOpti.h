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

#include "../../Module.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/BasicUtils/Thread.h>

namespace Nektar
{
namespace Utilities
{

const NekDouble dir[13][3] = {
    {  0.0,  0.0,  0.0 },  // 0  (x   , y   , z   )
    {  1.0,  0.0,  0.0 },  // 1  (x+dx, y   , z   )
    {  1.0,  1.0,  0.0 },  // 2  (x+dx, y+dy, z   )
    {  0.0,  1.0,  0.0 },  // 3  (x   , y+dy, z   )
    { -1.0,  0.0,  0.0 },  // 4  (x-dx, y   , z   )
    { -1.0, -1.0,  0.0 },  // 5  (x-dx, y-dy, z   )
    {  0.0, -1.0,  0.0 },  // 6  (x   , y-dy, z   )
    { -1.0,  0.0, -1.0 },  // 7  (x-dx, y   , z-dz)
    {  0.0,  0.0, -1.0 },  // 8  (x   , y   , z-dz)
    {  0.0,  0.0,  1.0 },  // 9  (x   , y   , z+dz)
    {  1.0,  0.0,  1.0 },  // 10 (x+dx, y   , z+dz)
    {  0.0,  1.0,  1.0 },  // 11 (x   , y+dy, z+dz)
    {  0.0, -1.0, -1.0 }   // 12 (x   , y-dy, z-dz)
};

class NodeOptiJob;

class NodeOpti
{
public:
    NodeOpti(NodeSharedPtr n,
             std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
             ResidualSharedPtr r, DerivUtilSharedPtr d,
             optimiser o)
        : node(n), nodeIds(e.first), data(e.second), res(r), derivUtil(d), opti(o)
    {
    }

    virtual ~NodeOpti(){};

    virtual void Optimise() = 0;
    NodeOptiJob *GetJob();

protected:
    virtual Array<OneD, NekDouble> GetGrad(bool analytic = false)
    {
        return Array<OneD,NekDouble>();
    }
    template<int DIM> NekDouble GetFunctional(bool analytic = false);
    NodeSharedPtr node;
    boost::mutex mtx;
    std::vector<int> nodeIds;
    std::vector<ElUtilSharedPtr> data;
    Array<OneD, NekDouble> grad;

    void CalcDX();
    void CalcMinJac();

    NekDouble dx;
    NekDouble minJac;
    ResidualSharedPtr res;
    DerivUtilSharedPtr derivUtil;
    optimiser opti;
};

typedef boost::shared_ptr<NodeOpti> NodeOptiSharedPtr;
typedef LibUtilities::NekFactory<int,
                                 NodeOpti,
                                 NodeSharedPtr,
                                 std::pair<std::vector<int>,
                                           std::vector<ElUtilSharedPtr> >,
                                 ResidualSharedPtr,
                                 DerivUtilSharedPtr,
                                 optimiser> NodeOptiFactory;

NodeOptiFactory &GetNodeOptiFactory();


class NodeOpti3D3D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti3D3D(NodeSharedPtr n,
                 std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 optimiser o)
                 : NodeOpti(n,e,r,d,o)
    {
    }

    ~NodeOpti3D3D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n,
        std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
        ResidualSharedPtr r, DerivUtilSharedPtr d,
        optimiser o)
    {
        return NodeOptiSharedPtr(new NodeOpti3D3D(n, e, r, d, o));
    }

private:
    Array<OneD, NekDouble> GetGrad(bool analytic = false);
};

class NodeOpti2D2D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti2D2D(NodeSharedPtr n,
                 std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 optimiser o)
                 : NodeOpti(n,e,r,d,o)
    {
    }

    ~NodeOpti2D2D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n, std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
        ResidualSharedPtr r, DerivUtilSharedPtr d,
        optimiser o)
    {
        return NodeOptiSharedPtr(new NodeOpti2D2D(n, e, r, d, o));
    }

private:
    Array<OneD, NekDouble> GetGrad(bool analytic = false);
};

class NodeOptiJob : public Thread::ThreadJob
{
public:
    NodeOptiJob(NodeOpti* no) : node(no) {}

    void Run()
    {
        node->Optimise();
    }
private:
    NodeOpti* node;
};

}
}

#endif
