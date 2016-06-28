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

#include "../../Module.h"
#include "ProcessVarOpti.h"

#include <NekMeshUtils/CADSystem/CADCurve.h>

namespace Nektar
{
namespace Utilities
{

class NodeOptiJob;

class NodeOpti
{
public:
    NodeOpti(NodeSharedPtr n, vector<ElDataSharedPtr> e,
             ResidualSharedPtr r, DerivUtilSharedPtr d,
             PtsHelperSharedPtr p, optimiser o)
             : node(n), data(e), res(r), derivUtil(d), ptsHelp(p), opti(o)
    {
    }

    virtual ~NodeOpti(){};

    virtual void Optimise(){};
    NodeOptiJob *GetJob();
protected:
    virtual Array<OneD, NekDouble> GetGrad()
    {
        return Array<OneD,NekDouble>();
    }
    template<int DIM> NekDouble GetFunctional();
    NodeSharedPtr node;
    vector<ElDataSharedPtr> data;

    ResidualSharedPtr res;
    DerivUtilSharedPtr derivUtil;
    PtsHelperSharedPtr ptsHelp;
    optimiser opti;
};

class NodeOpti1D3D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti1D3D(NodeSharedPtr n, vector<ElDataSharedPtr> e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 PtsHelperSharedPtr p, optimiser o, CADCurveSharedPtr c)
                 : NodeOpti(n,e,r,d,p,o), curve(c)
    {
    }

    ~NodeOpti1D3D(){};

    void Optimise();

private:
    Array<OneD, NekDouble> GetGrad();
    CADCurveSharedPtr curve;
};

class NodeOpti2D3D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti2D3D(NodeSharedPtr n, vector<ElDataSharedPtr> e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 PtsHelperSharedPtr p, optimiser o, CADSurfSharedPtr s)
                 : NodeOpti(n,e,r,d,p,o), surf(s)
    {
    }

    ~NodeOpti2D3D(){};

    void Optimise();

private:
    Array<OneD, NekDouble> GetGrad();
    CADSurfSharedPtr surf;
};

class NodeOpti3D3D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti3D3D(NodeSharedPtr n, vector<ElDataSharedPtr> e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 PtsHelperSharedPtr p, optimiser o)
                 : NodeOpti(n,e,r,d,p,o)
    {
    }

    ~NodeOpti3D3D(){};

    void Optimise();

private:
    Array<OneD, NekDouble> GetGrad();
};

class NodeOpti2D2D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti2D2D(NodeSharedPtr n, vector<ElDataSharedPtr> e,
                 ResidualSharedPtr r, DerivUtilSharedPtr d,
                 PtsHelperSharedPtr p, optimiser o)
                 : NodeOpti(n,e,r,d,p,o)
    {
    }

    ~NodeOpti2D2D(){};

    void Optimise();

private:
    Array<OneD, NekDouble> GetGrad();
};

}
}

#endif
