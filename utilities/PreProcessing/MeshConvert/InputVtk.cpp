////////////////////////////////////////////////////////////////////////////////
//
//  File: InputVtk.cpp
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
//  Description: VTK converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "MeshElements.h"
#include "InputVtk.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputVtk::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "vtk"), InputVtk::create,
                "Reads VTK format.");

        InputVtk::InputVtk(MeshSharedPtr m) : InputModule(m)
        {

        }

        InputVtk::~InputVtk()
        {

        }


        /**
         * Gmsh file contains a list of nodes and their coordinates, along with
         * a list of elements and those nodes which define them. We read in and
         * store the list of nodes in #m_node and store the list of elements in
         * #m_element. Each new element is supplied with a list of entries from
         * #m_node which defines the element. Finally some mesh statistics are
         * printed.
         *
         * @param   pFilename           Filename of Gmsh file to read.
         */
        void InputVtk::Process()
        {
            if (m->verbose)
            {
                cout << "InputVtk: Start reading file..." << endl;
            }

            vtkPolyDataReader *vtkMeshReader = vtkPolyDataReader::New();
            vtkMeshReader->SetFileName(config["infile"].as<string>().c_str());
            vtkMeshReader->Update();
            vtkPolyData *vtkMesh = vtkMeshReader->GetOutput();

            vtkPoints *vtkPoints = vtkMesh->GetPoints();

            const int numCellTypes = 2;
            vtkCellArray* vtkCells[numCellTypes];
            vtkCells[0] = vtkMesh->GetPolys();
            vtkCells[1] = vtkMesh->GetStrips();

            vtkIdType npts;
            vtkIdType *pts = 0;
            double p[3];

            for (int i = 0; i < vtkPoints->GetNumberOfPoints(); ++i)
            {
                vtkPoints->GetPoint(i, p);

                if ((p[0] * p[0]) > 0.000001 && m->spaceDim < 1)
                {
                    m->spaceDim = 1;
                }
                if ((p[1] * p[1]) > 0.000001 && m->spaceDim < 2)
                {
                    m->spaceDim = 2;
                }
                if ((p[2] * p[2]) > 0.000001 && m->spaceDim < 3)
                {
                    m->spaceDim = 3;
                }

                m->node.push_back(boost::shared_ptr<Node>(new Node(i, p[0], p[1], p[2])));
            }

            for (int c = 0; c < numCellTypes; ++c)
            {
                vtkCells[c]->InitTraversal();
                for (int i = 0; vtkCells[c]->GetNextCell(npts, pts); ++i)
                {
                    for (int j = 0; j < npts - 2; ++j)
                    {
                        // Create element tags
                        vector<int> tags;
                        tags.push_back(0); // composite
                        tags.push_back(eTriangle); // element type

                        // Read element node list
                        vector<NodeSharedPtr> nodeList;
                        for (int k = j; k < j + 3; ++k)
                        {
                            nodeList.push_back(m->node[pts[k]]);
                        }

                        // Create element
                        ElmtConfig conf(eTriangle,1,false,false);
                        ElementSharedPtr E = GetElementFactory().
                            CreateInstance(eTriangle,conf,nodeList,tags);

                        // Determine mesh expansion dimension
                        if (E->GetDim() > m->expDim) {
                            m->expDim = E->GetDim();
                        }
                        m->element[E->GetDim()].push_back(E);
                    }
                }
            }

            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
