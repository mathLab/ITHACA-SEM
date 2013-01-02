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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "MeshElements.h"
#include "OutputVtk.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputVtk::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "vtk"), OutputVtk::create,
                "Writes a VTK file.");

        OutputVtk::OutputVtk(MeshSharedPtr m) : OutputModule(m)
        {
            
        }

        OutputVtk::~OutputVtk()
        {

        }
        
        void OutputVtk::Process()
        {
            if (m->verbose)
            {
                cout << "OutputVtk: Writing file..." << endl;
            }

            vtkPolyData *vtkMesh = vtkPolyData::New();
            vtkPoints *vtkPoints = vtkPoints::New();
            vtkCellArray *vtkPolys = vtkCellArray::New();

            std::set<NodeSharedPtr>::iterator it;

            std::set<NodeSharedPtr> tmp(
                    m->vertexSet.begin(),
                    m->vertexSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                NodeSharedPtr n = *it;
                vtkPoints->InsertPoint(n->id, n->x, n->y, n->z);
            }

            vtkIdType p[8];
            vector<ElementSharedPtr> &elmt = m->element[m->expDim];
            for(int i = 0; i < elmt.size(); ++i)
            {
                int vertexCount = elmt[i]->GetVertexCount();
                for (int j = 0; j < vertexCount; ++j)
                {
                    p[j] = elmt[i]->GetVertex(j)->id;
                }
                vtkPolys->InsertNextCell(vertexCount, &p[0]);
            }

            vtkMesh->SetPoints(vtkPoints);
            vtkMesh->SetPolys(vtkPolys);

            // Write out the new mesh
            vtkPolyDataWriter *vtkMeshWriter = vtkPolyDataWriter::New();
            vtkMeshWriter->SetFileName(config["outfile"].as<string>().c_str());
            vtkMeshWriter->SetInput(vtkMesh);
            vtkMeshWriter->Update();
        }        
    }
}
