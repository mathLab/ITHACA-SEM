////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "../MeshElements.h"
#include "ProcessJac.h"

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

#include <vector>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessJac::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "jac"), ProcessJac::create,
                "Process elements based on values of Jacobian.");

        ProcessJac::ProcessJac(MeshSharedPtr m) : ProcessModule(m)
        {
            m_config["extract"] = ConfigOption(
                true, "0", "Extract non-valid elements from mesh.");
            m_config["removecurveifsingular"] = ConfigOption(
                true, "0", "remove curve nodes if element is singular.");
            m_config["list"]    = ConfigOption(
                true, "0", "Print list of elements having negative Jacobian.");
        }

        ProcessJac::~ProcessJac()
        {

        }

        void ProcessJac::Process()
        {
            if (m_mesh->m_verbose)
            {
                cout << "ProcessJac: Calculating Jacobians... " << endl;
            }

            bool extract = m_config["extract"].as<bool>();
            bool printList = m_config["list"].as<bool>();
            bool RemoveCurveIfSingular = m_config["removecurveifsingular"].as<bool>();
            vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

            if (extract)
            {
                m_mesh->m_element[m_mesh->m_expDim].clear();
            }

            if (printList)
            {
                cout << "Elements with negative Jacobian:" << endl;
            }

            int nNeg = 0;

            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                // Create elemental geometry.
                SpatialDomains::GeometrySharedPtr geom =
                    el[i]->GetGeom(m_mesh->m_spaceDim);

                // Generate geometric factors.
                SpatialDomains::GeomFactorsSharedPtr gfac =
                    geom->GetGeomFactors();

                // Get the Jacobian and, if it is negative, print a warning
                // message.
                if (!gfac->IsValid())
                {
                    nNeg++;

                    if (printList)
                    {
                        cout << "  - " << el[i]->GetId() << " ("
                             << StdRegions::ElementTypeMap[el[i]->GetConf().m_e]
                             << ")" << endl;
                    }

                    if(RemoveCurveIfSingular)
                    {
                        int nSurf = el[i]->GetFaceCount(); 
                        if(nSurf == 5) // prism mesh 
                        {
                            // find edges ndoes and blend to far side.
                            for(int e = 0; e < el[i]->GetEdgeCount(); ++e)
                            {
                                EdgeSharedPtr ed = el[i]->GetEdge(e);
                                EdgeSet::iterator it; 
                                // find edge in m_edgeSet; 
                                if((it = m_mesh->m_edgeSet.find(ed)) != m_mesh->m_edgeSet.end())
                                {
                                    if((*it)->m_edgeNodes.size())
                                    {
                                        vector<NodeSharedPtr> zeroNodes;
                                        (*it)->m_edgeNodes = zeroNodes;
                                    }
                                }
                            }
                        }
                    }

                    if (extract)
                    {
                        m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
                    }
                }
            }

            if (extract)
            {
                m_mesh->m_element[m_mesh->m_expDim-1].clear();
                ProcessVertices();
                ProcessEdges();
                ProcessFaces();
                ProcessElements();
                ProcessComposites();
            }

            if (printList || m_mesh->m_verbose)
            {
                cout << "Total negative Jacobians: " << nNeg << endl;
            }
            else if (nNeg > 0)
            {
                cout << "WARNING: Detected " << nNeg << " element"
                     << (nNeg == 1 ? "" : "s") << " with negative Jacobian."
                     << endl;
            }
        }
    }
}
