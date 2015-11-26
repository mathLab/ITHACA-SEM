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

#include <NekMeshUtils/MeshElements/MeshElements.h>

#include <SpatialDomains/MeshGraph.h>

#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

#include "ProcessJac.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

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
    m_config["file"]    = ConfigOption(
        false, "-1", "Write a file with a histogram of scaled Jacobians");
    m_config["roca"]    = ConfigOption(
        false, "-1", "Write a file with a histogram of scaled Jacobians");
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
    string fn = m_config["file"].as<string>();
    string rc = m_config["roca"].as<string>();
    bool file = !boost::iequals(fn,"-1");
    bool roca = !boost::iequals(rc,"-1");

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

    Array<OneD, int> bin(20,0);

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
                     << LibUtilities::ShapeTypeMap[el[i]->GetConf().m_e]
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

        if(roca)
        {
            StdRegions::StdExpansionSharedPtr chi = geom->GetXmap();
            LibUtilities::PointsKeyVector p = chi->GetPointsKeys();

            SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

            // get 4th order triangle points
            int order = 5;
            LibUtilities::PointsKey nodalKey(order, LibUtilities::eNodalTriEvenlySpaced);
            Array<OneD, NekDouble> x, y;
            LibUtilities::PointsManager()[nodalKey]->GetPoints(x,y);

            int numNodalPoints = x.num_elements();

            Array<ThreeD, NekDouble> jacmat(3,3,numNodalPoints,0.0);

            for (int i = 0; i < numNodalPoints; ++i)
            {
                for(int j = 0; j < 3; j++)
                {
                    for(int k = 0; k < 3; k++)
                    {
                        Array<OneD, NekDouble> coords(2);
                        coords[0] = x[i];
                        coords[1] = y[i];
                        jacmat[k][j][i] = chi->PhysEvaluate(coords, deriv[k][j]);
                    }

                }
            }

        }

        if(file)
        {
            LibUtilities::PointsKeyVector pkey =
                        geom->GetXmap()->GetPointsKeys();
            Array<OneD, NekDouble> jac = gfac->GetJac(pkey);
            int nq = jac.num_elements();

            NekDouble sjac = Vmath::Vmin(nq,jac,1) / Vmath::Vmax(nq,jac,1);
            if(!gfac->IsValid()) sjac = sjac*-1.0;

            if(sjac < -0.9) bin[0]++;
            else if(sjac < -0.8 && sjac >= -0.9) bin[1]++;
            else if(sjac < -0.7 && sjac >= -0.8) bin[2]++;
            else if(sjac < -0.6 && sjac >= -0.7) bin[3]++;
            else if(sjac < -0.5 && sjac >= -0.6) bin[4]++;
            else if(sjac < -0.4 && sjac >= -0.5) bin[5]++;
            else if(sjac < -0.3 && sjac >= -0.4) bin[6]++;
            else if(sjac < -0.2 && sjac >= -0.3) bin[7]++;
            else if(sjac < -0.1 && sjac >= -0.2) bin[8]++;
            else if(sjac < 0.0 && sjac >= -0.1) bin[9]++;
            else if(sjac < 0.1 && sjac >= 0.0) bin[10]++;
            else if(sjac < 0.2 && sjac >= 0.1) bin[11]++;
            else if(sjac < 0.3 && sjac >= 0.2) bin[12]++;
            else if(sjac < 0.4 && sjac >= 0.3) bin[13]++;
            else if(sjac < 0.5 && sjac >= 0.4) bin[14]++;
            else if(sjac < 0.6 && sjac >= 0.5) bin[15]++;
            else if(sjac < 0.7 && sjac >= 0.6) bin[16]++;
            else if(sjac < 0.8 && sjac >= 0.7) bin[17]++;
            else if(sjac < 0.9 && sjac >= 0.8) bin[18]++;
            else if(sjac <= 1.0 && sjac >= 0.9) bin[19]++;
        }
    }

    if(file)
    {
        ofstream out;
        out.open(fn.c_str());
        for(int i = 0; i < 20; i++)
        {
            out << -1.0+i*0.1+0.05 << " " << bin[i] << endl;
        }
        out.close();
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
