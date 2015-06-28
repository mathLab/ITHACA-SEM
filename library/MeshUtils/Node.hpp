////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_NODE_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_NODE_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace MeshUtils {
    
    class Node
    {
    public:
        
        Node(NekDouble x, NekDouble y, NekDouble z) :
                    m_x(x), m_y(y), m_z(z)
        {
        };
        
        void SetCurve(int i, NekDouble t)
        {
            CADCurve.push_back(i);
            CurveT.push_back(t);
        }

        Array<OneD, NekDouble> GetLoc()
        {
            Array<OneD, NekDouble> out(3);
            out[0]=m_x;
            out[1]=m_y;
            out[2]=m_z;
            return out;
        }
        
        
    private:
        
        NekDouble m_x, m_y, m_z;
        std::vector<int> CADCurve;
        std::vector<NekDouble> CurveT;
        std::vector<int> CADSurf;
        std::vector<std::vector<NekDouble> > SurfUV;
        
    };
        
}
}

#endif
