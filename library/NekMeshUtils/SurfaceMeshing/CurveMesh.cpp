////////////////////////////////////////////////////////////////////////////////
//
//  File: CurveMesh.cpp
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
//  Description: curvemesh object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void CurveMesh::Mesh()
{
    // this algorithm is mostly based on the work in chapter 19

    m_bounds      = m_cadcurve->GetBounds();
    m_curvelength = m_cadcurve->GetTotLength();
    m_numSamplePoints =
        int(m_curvelength / m_mesh->m_octree->GetMinDelta()) + 10;
    ds = m_curvelength / (m_numSamplePoints - 1);

    GetSampleFunction();

    Ae = 0.0;

    for (int i = 0; i < m_numSamplePoints - 1; i++)
    {
        Ae += ds * (1.0 / m_dst[i][0] + 1.0 / m_dst[i + 1][0]) / 2.0;
    }

    Ne = round(Ae);

    if (Ne + 1 < 2)
    {
        meshsvalue.resize(2);
        meshsvalue[0] = 0.0;
        meshsvalue[1] = m_curvelength;
        Ne            = 1;
    }
    else
    {

        GetPhiFunction();

        meshsvalue.resize(Ne + 1);
        meshsvalue[0]  = 0.0;
        meshsvalue[Ne] = m_curvelength;

        for (int i = 1; i < Ne; i++)
        {
            int iterationcounter = 0;
            bool iterate         = true;
            int k                = i;
            NekDouble ski        = meshsvalue[i - 1];
            NekDouble lastSki;
            while (iterate)
            {
                iterationcounter++;
                NekDouble rhs = EvaluateDS(ski) / Ae * (EvaluatePS(ski) - k);
                lastSki       = ski;
                ski           = ski - rhs;
                if (abs(lastSki - ski) < 1E-8)
                {
                    iterate = false;
                }

                ASSERTL0(iterationcounter < 1000000, "iteration failed");
            }

            meshsvalue[i] = ski;
        }
    }

    NekDouble t;
    Array<OneD, NekDouble> loc;

    vector<CADVertSharedPtr> verts = m_cadcurve->GetVertex();
    vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> > s = m_cadcurve->GetAdjSurf();

    NodeSharedPtr n = verts[0]->GetNode();
    t               = m_bounds[0];
    n->SetCADCurve(m_id, m_cadcurve, t);
    loc = n->GetLoc();
    for (int j = 0; j < s.size(); j++)
    {
        if (verts[0]->IsDegen() == s[j].first->GetId())
        {
            // if the degen has been set for this node the node
            // already knows its corrected location
            continue;
        }

        Array<OneD, NekDouble> uv = s[j].first->locuv(loc);
        n->SetCADSurf(s[j].first->GetId(), s[j].first, uv);
    }
    m_meshpoints.push_back(n);

    for (int i = 1; i < meshsvalue.size() - 1; i++)
    {
        t                = m_cadcurve->tAtArcLength(meshsvalue[i]);
        loc              = m_cadcurve->P(t);
        NodeSharedPtr n2 = boost::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, loc[0], loc[1], loc[2]));
        n2->SetCADCurve(m_id, m_cadcurve, t);
        for (int j = 0; j < s.size(); j++)
        {
            Array<OneD, NekDouble> uv = s[j].first->locuv(loc);
            n2->SetCADSurf(s[j].first->GetId(), s[j].first, uv);
        }
        m_meshpoints.push_back(n2);
    }

    n = verts[1]->GetNode();
    t = m_bounds[1];
    n->SetCADCurve(m_id, m_cadcurve, t);
    loc = n->GetLoc();
    for (int j = 0; j < s.size(); j++)
    {
        if (verts[1]->IsDegen() == s[j].first->GetId())
        {
            // if the degen has been set for this node the node
            // already knows its corrected location
            continue;
        }

        Array<OneD, NekDouble> uv = s[j].first->locuv(loc);
        n->SetCADSurf(s[j].first->GetId(), s[j].first, uv);
    }
    m_meshpoints.push_back(n);

    ASSERTL0(Ne + 1 == m_meshpoints.size(),
             "incorrect number of points in curve mesh");

    // make edges and add them to the edgeset for the face mesher to use
    for (int i = 0; i < m_meshpoints.size() - 1; i++)
    {
        EdgeSharedPtr e = boost::shared_ptr<Edge>(
            new Edge(m_meshpoints[i], m_meshpoints[i + 1]));
        e->m_parentCAD = m_cadcurve;
        m_mesh->m_edgeSet.insert(e);
        m_meshedges.push_back(e);
    }

    if (m_mesh->m_verbose)
    {
        cout << "\r                                                            "
                "    "
                "                             ";
        cout << scientific << "\r\t\tCurve " << m_id << endl
             << "\t\t\tLength: " << m_curvelength << endl
             << "\t\t\tNodes: " << m_meshpoints.size() << endl
             << "\t\t\tSample points: " << m_numSamplePoints << endl
             << endl;
    }
}

void CurveMesh::GetPhiFunction()
{
    m_ps.resize(m_numSamplePoints);
    vector<NekDouble> newPhi;
    newPhi.resize(2);

    newPhi[0] = 0.0;
    newPhi[1] = 0.0;

    m_ps[0] = newPhi;

    NekDouble runningInt = 0.0;

    for (int i = 1; i < m_numSamplePoints; i++)
    {
        runningInt += (1.0 / m_dst[i - 1][0] + 1.0 / m_dst[i][0]) / 2.0 * ds;
        newPhi[0] = Ne / Ae * runningInt;
        newPhi[1] = m_dst[i][1];
        m_ps[i]   = newPhi;
    }
}

NekDouble CurveMesh::EvaluateDS(NekDouble s)
{
    int a = 0;
    int b = 0;

    ASSERTL1(!(s < 0) && !(s > m_curvelength),"s out of bounds");

    if (s == 0)
    {
        return m_dst[0][0];
    }
    else if (s == m_curvelength)
    {
        return m_dst[m_numSamplePoints - 1][0];
    }

    for (int i = 0; i < m_numSamplePoints - 1; i++)
    {
        if (m_dst[i][1] < s && m_dst[i + 1][1] >= s)
        {
            a = i;
            b = i + 1;
            break;
        }
    }

    NekDouble s1 = m_dst[a][1];
    NekDouble s2 = m_dst[b][1];
    NekDouble d1 = m_dst[a][0];
    NekDouble d2 = m_dst[b][0];

    NekDouble m = (d2 - d1) / (s2 - s1);
    NekDouble c = d2 - m * s2;

    ASSERTL0(m * s + c == m * s + c, "DS"); // was getting nans here

    return m * s + c;
}

NekDouble CurveMesh::EvaluatePS(NekDouble s)
{
    int a = 0;
    int b = 0;

    ASSERTL1(!(s < 0) && !(s > m_curvelength),"s out of bounds");

    if (s == 0)
    {
        return m_ps[0][0];
    }
    else if (s == m_curvelength)
    {
        return m_ps[m_numSamplePoints - 1][0];
    }

    for (int i = 0; i < m_numSamplePoints - 1; i++)
    {
        if (m_ps[i][1] < s && m_ps[i + 1][1] >= s)
        {
            a = i;
            b = i + 1;
            break;
        }
    }

    if (a == b)
    {
        cout << endl;
        cout << a << " " << b << endl;
        cout << s << endl;
        exit(-1);
    }

    NekDouble s1 = m_ps[a][1];
    NekDouble s2 = m_ps[b][1];
    NekDouble d1 = m_ps[a][0];
    NekDouble d2 = m_ps[b][0];

    NekDouble m = (d2 - d1) / (s2 - s1);
    NekDouble c = d2 - m * s2;

    ASSERTL0(m * s + c == m * s + c, "PS");

    return m * s + c;
}

void CurveMesh::GetSampleFunction()
{
    m_dst.resize(m_numSamplePoints);

    vector<NekDouble> dsti;
    dsti.resize(3);

    for (int i = 0; i < m_numSamplePoints; i++)
    {
        dsti[1]     = i * ds;
        NekDouble t = m_cadcurve->tAtArcLength(dsti[1]);

        Array<OneD, NekDouble> loc = m_cadcurve->P(t);

        /*NekDouble ts =
            m_bl.Evaluate(m_blID, loc[0], loc[1], loc[2], 0.0);

        if (ts > 0.0)
        {
            Array<OneD, NekDouble> N = m_cadcurve->N(t);
            Array<OneD, NekDouble> Nwrt = m_cadcurve->NormalWRT(t, 0);

            if(N[0]*N[0] + N[1]*N[1] + N[2]*N[2] < 1e-6)
            {
                dsti[0] = m_mesh->m_octree->Query(loc);
            }
            else if ( N[0]*Nwrt[0] + N[1]*Nwrt[1] + N[2]*Nwrt[2] > 0)
            {
                //concave
                dsti[0] = m_mesh->m_octree->Query(loc);
            }
            else
            {
                NekDouble R = 1.0 / m_cadcurve->Curvature(t);
                if(R > 2.0*t)
                {
                    R = 2.0*t;
                }
                Array<OneD, NekDouble> tloc(3);
                tloc[0] = loc[0] + ts * Nwrt[0];
                tloc[1] = loc[1] + ts * Nwrt[1];
                tloc[2] = loc[2] + ts * Nwrt[2];

                NekDouble d = m_mesh->m_octree->Query(tloc);

                dsti[0] = d * R / (R + ts);
            }
        }
        else
        {*/
            dsti[0] = m_mesh->m_octree->Query(loc);
        //}

        dsti[2] = t;

        m_dst[i] = dsti;
    }
}
}
}
