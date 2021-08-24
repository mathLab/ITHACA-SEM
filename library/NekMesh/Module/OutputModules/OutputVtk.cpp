////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.cpp
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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMesh/MeshElements/Element.h>
#include <LibUtilities/BasicUtils/VtkUtil.hpp>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellType.h>

#include "OutputVtk.h"

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar
{
namespace NekMesh
{

ModuleKey OutputVtk::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "vtk"), OutputVtk::create, "Writes a VTK file.");

OutputVtk::OutputVtk(MeshSharedPtr m) : OutputModule(m)
{
    m_config["uncompress"] = ConfigOption(true, "0", "Uncompress xml sections");
    m_config["legacy"]        = ConfigOption(true, "0", "Output in legacy format");
}

OutputVtk::~OutputVtk()
{
}

int OutputVtk::GetVtkCellType(std::string pType)
{
    if (pType == "S")
    {
        return VTK_LINE;
    }
    else if (pType == "T")
    {
        return VTK_TRIANGLE;
    }
    else if (pType == "Q")
    {
        return VTK_QUAD;
    }
    else if (pType == "A")
    {
        return VTK_TETRA;
    }
    else if (pType == "P")
    {
        return VTK_PYRAMID;
    }
    else if (pType == "R")
    {
        return VTK_WEDGE;
    }
    else if (pType == "H")
    {
        return VTK_HEXAHEDRON;
    }
    else
    {
        ASSERTL0(false, "Element type not supported.");
        return 0;
    }
}

void OutputVtk::Process()
{
    m_log(VERBOSE) << "Writing VTK file '"
                   << m_config["outfile"].as<string>() << "'." << endl;

    vtkUnstructuredGrid *vtkMesh   = vtkUnstructuredGrid::New();
    vtkPoints           *vtkPoints = vtkPoints::New();

    std::set<NodeSharedPtr> tmp(m_mesh->m_vertexSet.begin(),
                                m_mesh->m_vertexSet.end());

    for (auto &n : tmp)
    {
        vtkPoints->InsertPoint(n->m_id, n->m_x, n->m_y, n->m_z);
    }

    vtkIdType p[8];
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];
    for (int i = 0; i < elmt.size(); ++i)
    {
        int vertexCount = elmt[i]->GetVertexCount();
        for (int j = 0; j < vertexCount; ++j)
        {
            p[j] = elmt[i]->GetVertex(j)->m_id;
        }
        // Adjust vertex order to the vtk convention
        if (elmt[i]->GetTag() == "R")
        {
            std::swap(p[2], p[4]);
        }
        vtkMesh->InsertNextCell(GetVtkCellType(elmt[i]->GetTag()),
                                vertexCount, &p[0]);
    }

    vtkMesh->SetPoints(vtkPoints);

    // Write out the new mesh in XML or legacy format
    if (m_config["legacy"].beenSet)
    {
        vtkUnstructuredGridWriter *vtkMeshWriter = vtkUnstructuredGridWriter::New();
        vtkMeshWriter->SetFileName(m_config["outfile"].as<string>().c_str());

#if VTK_MAJOR_VERSION <= 5
        vtkMeshWriter->SetInput(vtkMesh);
#else
        vtkMeshWriter->SetInputData(vtkMesh);
#endif
        vtkMeshWriter->Update();
    }
    else // XML format
    {
        vtkXMLUnstructuredGridWriter *vtkMeshWriter = vtkXMLUnstructuredGridWriter::New();
        vtkMeshWriter->SetFileName(m_config["outfile"].as<string>().c_str());

#if VTK_MAJOR_VERSION <= 5
        vtkMeshWriter->SetInput(vtkMesh);
#else
        vtkMeshWriter->SetInputData(vtkMesh);
#endif
        if (m_config["uncompress"].beenSet)
        {
            vtkMeshWriter->SetDataModeToAscii();
        }
        vtkMeshWriter->Update();
    }
}
}
}
