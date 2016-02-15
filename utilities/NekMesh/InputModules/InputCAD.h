////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.h
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
//  Description: Create mesh from CAD.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_INPUTCAD
#define UTILITIES_NEKMESH_INPUTCAD

#include "../Module.h"

namespace Nektar
{
namespace Utilities
{

class InputCAD : public InputModule
{
public:
    InputCAD(MeshSharedPtr m);
    virtual ~InputCAD();
    virtual void Process();

    /// Creates an instance of this class
    static ModuleSharedPtr create(MeshSharedPtr m)
    {
        return MemoryManager<InputCAD>::AllocateSharedPtr(m);
    }
    /// %ModuleKey for class.
    static ModuleKey className;

private:
    NekDouble m_minDelta, m_maxDelta, m_eps, m_blthick;
    int m_order;
    string m_CADName, m_udsName;
    bool m_makeBL, m_writeoctree;
};
}
}

#endif
