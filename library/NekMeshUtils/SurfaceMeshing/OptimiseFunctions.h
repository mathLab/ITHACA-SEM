////////////////////////////////////////////////////////////////////////////////
//
//  File: Curvemesh.h
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
//  Description: object for individual curve meshes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_SURFACEMESH_OPTIMISEFUNCTIONS_H
#define NEKTAR_MESHUTILS_SURFACEMESH_OPTIMISEFUNCTIONS_H

#include <LocalRegions/MatrixKey.h>
#include <NekMeshUtils/CADSystem/CADObject.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/Optimisation/OptimiseObj.h>

namespace Nektar
{
namespace NekMeshUtils
{

class OptiEdge : public OptiObj
{
public:
    friend class MemoryManager<OptiEdge>;

    OptiEdge(Array<OneD, NekDouble> a,
             Array<OneD, NekDouble> dis,
             CADObjectSharedPtr ob)
    {
        all = a;
        z   = dis;
        o   = ob;
    };

    ~OptiEdge(){};

    NekDouble F(Array<OneD, NekDouble> xitst);
    DNekMat dF(Array<OneD, NekDouble> xitst);
    Array<OneD, NekDouble> Getxi();
    Array<OneD, NekDouble> Getli();
    Array<OneD, NekDouble> Getui();
    void Update(Array<OneD, NekDouble> xinew);

    Array<OneD, NekDouble> GetSolution()
    {
        return all;
    };

private:
    CADObjectSharedPtr o;
    Array<OneD, NekDouble> z;
    Array<OneD, NekDouble> all;
};
typedef std::shared_ptr<OptiEdge> OptiEdgeSharedPtr;

}
}
#endif
