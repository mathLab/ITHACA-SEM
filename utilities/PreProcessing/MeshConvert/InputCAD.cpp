////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <LibUtilities/CADSystem/CADSystem.h>

#include "MeshElements.h"
#include "InputCAD.h"

namespace Nektar
{
namespace Utilities
{
    
    
    ModuleKey InputCAD::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "CAD"), InputCAD::create,
        "Reads CAD geometry and will generate the mesh file.");
    
    /**
     * @brief Set up InputCAD object.
     */
    InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
    {
        m_config["min"] = ConfigOption(false,"-1",
                "minimum delta to be in mesh");
        m_config["max"] = ConfigOption(false,"-1",
                "maximum delta to be in mesh");
        m_config["eps"] = ConfigOption(false,"-1",
                "sensitivity to curvature");
        m_config["order"] = ConfigOption(false,"-1",
                "order of the mesh to be produced");
    }
    
    InputCAD::~InputCAD()
    {
        
    }
    
    
    void InputCAD::Process()
    {
        string CADName = m_config["infile"].as<string>();
        NekDouble m_minDelta = m_config["min"].as<NekDouble>();
        NekDouble m_maxDelta = m_config["max"].as<NekDouble>();
        NekDouble m_eps = m_config["eps"].as<NekDouble>();
        int m_order = m_config["order"].as<int>();
        
        ASSERTL0(!(m_minDelta == -1 || m_maxDelta == -1 ||
                   m_eps == -1 || m_order == -1),
                 "User parameters required");
        
        
        LibUtilities::CADSystemSharedPtr CAD =
            MemoryManager<LibUtilities::CADSystem>::AllocateSharedPtr(CADName);
        
        cout << CAD->GetName() << endl;

        ASSERTL0(CAD->LoadCAD(),
                 "Failed to load CAD");
        
        if(m_mesh->m_verbose)
        {
            cout << "min delta: " << m_minDelta << " max delta: " << m_maxDelta
                 << " esp: " << m_eps << " order: " << m_order << endl;
            CAD->Report();
        }
        
    }
    
    
}
}