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

void CurveMesh::ReMesh()
{
    m_meshpoints.clear();
    m_dst.clear();
    m_ps.clear();
    meshsvalue.clear();
    for(int i = 0; i < m_meshedges.size(); i++)
    {
        m_mesh->m_edgeSet.erase(m_meshedges[i]);
    }
    m_meshedges.clear();

    Mesh(true);
}

void CurveMesh::Mesh(bool forceThree)
{
    // this algorithm is mostly based on the work in chapter 19

    m_bounds      = m_cadcurve->GetBounds();
    m_curvelength = m_cadcurve->GetTotLength();
    m_numSamplePoints =
        int(m_curvelength / m_mesh->m_octree->GetMinDelta()) + 10;
    ds = m_curvelength / (m_numSamplePoints - 1);

    // compute the offset due to adjacent BLs
    NekDouble totalOffset = 0.0;
    for (map<unsigned, NekDouble>::iterator ie = m_endoffset.begin();
         ie != m_endoffset.end(); ++ie)
    {
        totalOffset += ie->second;
    }
    ASSERTL0(m_curvelength > totalOffset,
             "Boundary layers too thick for adjacent curve");

    GetSampleFunction();

    Ae = 0.0;

    for (int i = 0; i < m_numSamplePoints - 1; i++)
    {
        Ae += ds * (1.0 / m_dst[i][0] + 1.0 / m_dst[i + 1][0]) / 2.0;
    }

    Ne = round(Ae);

    if (Ne + 1 < 2 + m_endoffset.size())
    {
        Ne = 1 + m_endoffset.size();

        meshsvalue.resize(Ne + 1);
        meshsvalue[0] = 0.0;
        meshsvalue[1] = m_curvelength;

        if (m_endoffset.count(0))
        {
            meshsvalue[1] = m_endoffset[0];
        }
        if (m_endoffset.count(1))
        {
            meshsvalue[Ne - 1] = m_curvelength - m_endoffset[1];
        }
    }
    else if(Ne + 1 == 2 && forceThree)
    {
        Ne++;
        meshsvalue.resize(Ne + 1);
        meshsvalue[0] = 0.0;
        meshsvalue[1] = m_curvelength/ 2.0;
        meshsvalue[2] = m_curvelength;
    }
    else
    {

        GetPhiFunction();

        meshsvalue.resize(Ne + 1);
        meshsvalue[0]  = 0.0;
        meshsvalue[Ne] = m_curvelength;

        // force the second and/or the second to last point(s) if an offset is
        // defined
        if (m_endoffset.count(0))
        {
            meshsvalue[1] = m_endoffset[0];
        }
        if (m_endoffset.count(1))
        {
            meshsvalue[Ne - 1] = m_curvelength - m_endoffset[1];
        }

        for (int i = 1 + m_endoffset.count(0); i < Ne - m_endoffset.count(1);
             i++)
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
    vector<pair<weak_ptr<CADSurf>, CADOrientation::Orientation> > s =
        m_cadcurve->GetAdjSurf();

    NodeSharedPtr n = verts[0]->GetNode();
    t               = m_bounds[0];
    n->SetCADCurve(m_cadcurve, t);
    loc = n->GetLoc();
    for (int j = 0; j < s.size(); j++)
    {
        if (verts[0]->IsDegen() == s[j].first.lock()->GetId())
        {
            // if the degen has been set for this node the node
            // already knows its corrected location
            continue;
        }

        Array<OneD, NekDouble> uv = s[j].first.lock()->locuv(loc);
        n->SetCADSurf(s[j].first.lock(), uv);
    }
    m_meshpoints.push_back(n);

    for (int i = 1; i < meshsvalue.size() - 1; i++)
    {
        t                = m_cadcurve->tAtArcLength(meshsvalue[i]);
        loc              = m_cadcurve->P(t);
        NodeSharedPtr n2 = std::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, loc[0], loc[1], loc[2]));
        n2->SetCADCurve(m_cadcurve, t);
        for (int j = 0; j < s.size(); j++)
        {
	    Array<OneD, NekDouble> uv = s[j].first.lock()->locuv(loc);
            n2->SetCADSurf(s[j].first.lock(), uv);
        }
        m_meshpoints.push_back(n2);
    }

    n = verts[1]->GetNode();
    t = m_bounds[1];
    n->SetCADCurve(m_cadcurve, t);
    loc = n->GetLoc();
    for (int j = 0; j < s.size(); j++)
    {
        if (verts[1]->IsDegen() == s[j].first.lock()->GetId())
        {
            // if the degen has been set for this node the node
            // already knows its corrected location
            continue;
        }

        Array<OneD, NekDouble> uv = s[j].first.lock()->locuv(loc);
        n->SetCADSurf(s[j].first.lock(), uv);
    }
    m_meshpoints.push_back(n);

    ASSERTL0(Ne + 1 == m_meshpoints.size(),
             "incorrect number of points in curve mesh");

    // make edges and add them to the edgeset for the face mesher to use
    for (int i = 0; i < m_meshpoints.size() - 1; i++)
    {
        EdgeSharedPtr e = std::shared_ptr<Edge>(
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

    ASSERTL1(!(s < 0)&& !(s > m_curvelength), "s out of bounds");

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

    ASSERTL1(!(s < 0) && !(s > m_curvelength), "s out of bounds");

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

        bool found = false;

        // if inside the BL, dsti[0] set to the BL thickness, i.e. the offset
        if (m_endoffset.count(0))
        {
            if (dsti[1] < m_endoffset[0])
            {
                dsti[0] = m_endoffset[0];
                found   = true;
            }
        }
        if (m_endoffset.count(1) && !found)
        {
            if (dsti[1] > m_curvelength - m_endoffset[1])
            {
                dsti[0] = m_endoffset[1];
                found   = true;
            }
        }
        // else, dsti[0] is found from the octree
        if (!found)
        {
            dsti[0] = m_mesh->m_octree->Query(loc);
        }

        dsti[2] = t;

        m_dst[i] = dsti;
    }
}

void CurveMesh::PeriodicOverwrite(CurveMeshSharedPtr from)
{
    // clear current mesh points and remove edges from edgeset
    m_meshpoints.clear();
    for (int i = 0; i < m_meshedges.size(); i++)
    {
        m_mesh->m_edgeSet.erase(m_meshedges[i]);
    }
    m_meshedges.clear();

    ///////

    int tid = from->GetId();
    Array<OneD, NekDouble> T =
        m_mesh->m_cad->GetPeriodicTranslationVector(tid, m_id);

    CADCurveSharedPtr c1 = m_mesh->m_cad->GetCurve(tid);

    bool reversed = c1->GetOrienationWRT(1) == m_cadcurve->GetOrienationWRT(1);

    vector<NodeSharedPtr> nodes = from->GetMeshPoints();

    vector<pair<weak_ptr<CADSurf>, CADOrientation::Orientation> > surfs =
        m_cadcurve->GetAdjSurf();

    for (int i = 1; i < nodes.size() - 1; i++)
    {
        Array<OneD, NekDouble> loc = nodes[i]->GetLoc();
        NodeSharedPtr nn = NodeSharedPtr(
            new Node(m_mesh->m_numNodes++, loc[0] + T[0], loc[1] + T[1], 0.0));

        for (int j = 0; j < surfs.size(); j++)
        {
            Array<OneD, NekDouble> uv = surfs[j].first.lock()->locuv(nn->GetLoc());
            nn->SetCADSurf(surfs[j].first.lock(), uv);
        }

        NekDouble t;
        m_cadcurve->loct(nn->GetLoc(), t);
        nn->SetCADCurve(m_cadcurve, t);

        m_meshpoints.push_back(nn);
    }

    // Reverse internal nodes of the vector if necessary
    if (reversed)
    {
        reverse(m_meshpoints.begin(), m_meshpoints.end());
    }

    vector<CADVertSharedPtr> verts = m_cadcurve->GetVertex();

    m_meshpoints.insert(m_meshpoints.begin(), verts[0]->GetNode());
    m_meshpoints.push_back(verts[1]->GetNode());
    // dont need to realign cad for vertices

    // make edges and add them to the edgeset for the face mesher to use
    for (int i = 0; i < m_meshpoints.size() - 1; i++)
    {
        EdgeSharedPtr e = std::shared_ptr<Edge>(
            new Edge(m_meshpoints[i], m_meshpoints[i + 1]));
        e->m_parentCAD = m_cadcurve;
        m_mesh->m_edgeSet.insert(e);
        m_meshedges.push_back(e);
    }
}
}
}
