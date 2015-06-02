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

#include "MeshElements.h"
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
            m_config["deformprismifsingular"] = ConfigOption(
                true, "0", "deform internal prism face if singular.");
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
            bool DeformPrismIfSingular = m_config["deformprismifsingular"].as<bool>();

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

                    if(DeformPrismIfSingular)
                    {
                        int nSurf = el[i]->GetFaceCount(); 
                        
                        if(nSurf == 5) // prism mesh 
                        {
                            int triedges[3][2] = {{0,2},{5,6},{4,7}};

                            // find edges ndoes and blend to far side.
                            for(int e = 0; e < 3; ++e)
                            {
                                EdgeSharedPtr e1 = el[i]->GetEdge(triedges[e][0]);
                                EdgeSharedPtr e2 = el[i]->GetEdge(triedges[e][1]);
                                EdgeSharedPtr efrom,eto;
                                
                                if(e1->m_edgeNodes.size() == 0)
                                {
                                    efrom = e2;
                                    eto   = e1;
                                }
                                else if(e2->m_edgeNodes.size() == 0)
                                {
                                    efrom = e1;
                                    eto   = e2;
                                }
                                else
                                {
                                    continue;
                                }
                                
                                // determine vector fron node 1 to 2
                                NodeSharedPtr nfrom1 = efrom->m_n1;
                                NodeSharedPtr nfrom2 = efrom->m_n2;

                                // determine vector fron node 1 to 2 of 
                                NodeSharedPtr nto1 = eto->m_n1;
                                NodeSharedPtr nto2 = eto->m_n2;

                                Node n1to2   = *nto2   - *nto1;
                                Node n1from2 = *nfrom2 - *nfrom1;
                                
                                int ne = efrom->m_edgeNodes.size();
                                vector<int> mapnode(ne);

                                //edges are reversed. 
                                if(n1to2.dot(n1from2) < 0)
                                {
                                    for(int j = 0; j < ne; ++j)
                                    {
                                        mapnode[j] = ne-j-1;
                                    }
                                    nto1 = eto->m_n2;
                                    nto2 = eto->m_n1;
                                }
                                else
                                {
                                    for(int j = 0; j < ne; ++j)
                                    {
                                        mapnode[j] = j;
                                    }
                                }                                


                                SpatialDomains::GeometrySharedPtr edgeom = efrom->GetGeom(m_mesh->m_spaceDim);

                                StdRegions::StdExpansionSharedPtr xmap = edgeom->GetXmap();
                                Array<OneD, Array<OneD,  NekDouble> > xc(m_mesh->m_spaceDim);
                                int nq = xmap->GetTotPoints();
                                for(int j = 0; j < m_mesh->m_spaceDim; ++j)
                                {
                                    xc[j] = Array<OneD, NekDouble>(nq);
                                }
                                
                                Array<OneD, NekDouble> coeffs(xmap->GetNcoeffs());

                                edgeom->FillGeom();
                                // extract coefficients from modal
                                // space replace the node values and
                                // backward transform.

                                // X coord
                                Vmath::Vcopy(xmap->GetNcoeffs(),
                                             edgeom->GetCoeffs(0),1,
                                             coeffs,1);
                                
                                // replace vertices with to values; 
                                coeffs[0] = nto1->m_x;
                                coeffs[1] = nto2->m_x;
                                    
                                // bwdtrans
                                xmap->BwdTrans(coeffs,xc[0]);

                                // Y coord
                                Vmath::Vcopy(xmap->GetNcoeffs(),
                                             edgeom->GetCoeffs(1),1,
                                             coeffs,1);
                                
                                // replace vertices with to values; 
                                coeffs[0] = nto1->m_y;
                                coeffs[1] = nto2->m_y;
                                    
                                // bwdtrans
                                xmap->BwdTrans(coeffs,xc[1]);

                                // Z coord
                                Vmath::Vcopy(xmap->GetNcoeffs(),
                                             edgeom->GetCoeffs(2),1,
                                             coeffs,1);
                                
                                // replace vertices with to values; 
                                coeffs[0] = nto1->m_z;
                                coeffs[1] = nto2->m_z;
                                    
                                // bwdtrans
                                xmap->BwdTrans(coeffs,xc[2]);
                                
                                
                                for(int j = 0; j < ne; ++j)
                                {
                                    NodeSharedPtr n = boost::shared_ptr<Node>(new Node());
                                    
                                    (*n).m_x = xc[0][mapnode[j]+1];
                                    (*n).m_y = xc[1][mapnode[j]+1];
                                    (*n).m_z = xc[2][mapnode[j]+1];
                                    
                                    eto->m_edgeNodes.push_back(n);
                                }
                                
                            }
                            
                        }
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

            if(DeformPrismIfSingular)
            {
                ProcessEdges();
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
