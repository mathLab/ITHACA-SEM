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

#include <MeshUtils/SurfaceMeshing/CurveMesh.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

void CurveMesh::Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                     std::map<int, MeshEdgeSharedPtr> &Edges)
{
    m_bounds = m_cadcurve->Bounds();
    m_curvelength = m_cadcurve->GetTotLength();

    m_numSamplePoints = int(m_curvelength/m_octree->GetMinDelta())+5;

    ds = m_curvelength/(m_numSamplePoints-1);

    GetSampleFunction();

    Ae = 0.0;

    for(int i = 0; i < m_numSamplePoints-1; i++)
    {
        Ae+=ds*(1.0/m_dst[i][0]+1.0/m_dst[i+1][0])/2.0;
    }

    Ne=round(Ae);

    if(Ne+1<2)
    {
        meshsvalue.resize(2);
        meshsvalue[0]=0.0;
        meshsvalue[1]=m_curvelength;
        Ne=1;
    }
    else
    {

        GetPhiFunction();

        meshsvalue.resize(Ne+1);
        meshsvalue[0]=0.0;
        meshsvalue[Ne]=m_curvelength;

        for(int i = 1; i < Ne; i++)
        {
            int iterationcounter=0;
            bool iterate = true;
            int k = i;
            NekDouble ski = meshsvalue[i-1];
            NekDouble lastSki;
            while(iterate)
            {
                iterationcounter++;
                NekDouble rhs =EvaluateDS(ski)/Ae*(EvaluatePS(ski)-k);
                lastSki=ski;
                ski=ski-rhs;
                if(abs(lastSki-ski)<1E-10)
                {
                    iterate = false;
                }

                ASSERTL0(iterationcounter<1000000, "iteration failed");
            }

            meshsvalue[i]=ski;
        }
    }

    NekDouble t;
    Array<OneD, NekDouble> loc;

    vector<int> endNodeIds = m_cadcurve->GetVertex();

    //need a better fix to this
    if(m_id == 19)
    {
        reverse(endNodeIds.begin(),endNodeIds.end());
    }

    t = m_bounds[0];
    Nodes[endNodeIds[0]]->SetCurve(m_id,t);
    m_meshpoints.push_back(endNodeIds[0]);

    for(int i = 1; i < meshsvalue.size()-1; i++)
    {
        t = m_cadcurve->tAtArcLength(meshsvalue[i]);
        loc = m_cadcurve->P(t);
        MeshNodeSharedPtr n2 = boost::shared_ptr<MeshNode>(
                          new MeshNode(Nodes.size(),loc[0],loc[1],loc[2]));
        n2->SetCurve(m_id,t);
        Nodes[Nodes.size()] = n2;
        m_meshpoints.push_back(Nodes.size()-1);
    }

    t = m_bounds[1];
    Nodes[endNodeIds[1]]->SetCurve(m_id,t);
    m_meshpoints.push_back(endNodeIds[1]);

    for(int i = 0; i < m_meshpoints.size()-1; i++)
    {
        Edges[Edges.size()] = MemoryManager<MeshEdge>::
                AllocateSharedPtr(Edges.size(),
                                  m_meshpoints[i],
                                  m_meshpoints[i+1]);
        Nodes[m_meshpoints[i]]->SetEdge(Edges.size()-1);
        Nodes[m_meshpoints[i+1]]->SetEdge(Edges.size()-1);
        Edges[Edges.size()-1]->SetCurve(m_id);
    }




}

void CurveMesh::Report()
{
    cout << scientific << "\tLength: " << m_curvelength << "\tPoints: " << Ne+1 << endl;
}

void CurveMesh::GetPhiFunction()
{
    m_ps.resize(m_numSamplePoints);
    vector<NekDouble> newPhi;
    newPhi.resize(2);

    newPhi[0]=0.0;
    newPhi[1]=0.0;

    m_ps[0]=newPhi;

    NekDouble runningInt=0.0;

    for(int i = 1; i < m_numSamplePoints; i++)
    {
        runningInt+=(1.0/m_dst[i-1][0]+1.0/m_dst[i][0])/2.0*ds;
        newPhi[0]=Ne/Ae*runningInt;
        newPhi[1]=m_dst[i][1];
        m_ps[i]=newPhi;
    }
}

NekDouble CurveMesh::EvaluateDS(NekDouble s)
{
    int a=0;
    int b=0;

    if(s==0)
    {
        return m_dst[0][0];
    }
    else if(s==m_curvelength)
    {
        return m_dst[m_numSamplePoints-1][0];
    }

    for(int i = 0; i < m_numSamplePoints-1; i++)
    {
        if(m_dst[i][1]<s && m_dst[i+1][1]>=s)
        {
            a=i;
            b=i+1;
            break;
        }
    }

    NekDouble s1 = m_dst[a][1];
    NekDouble s2 = m_dst[b][1];
    NekDouble d1 = m_dst[a][0];
    NekDouble d2 = m_dst[b][0];

    NekDouble m  = (d2-d1)/(s2-s1);
    NekDouble c = d2- m*s2;

    ASSERTL0(m*s+c==m*s+c,"DS"); //was getting nans here

    return m*s+c;
}

NekDouble CurveMesh::EvaluatePS(NekDouble s)
{
    int a=0;
    int b=0;

    if(s==0)
    {
        return m_ps[0][0];
    }
    else if(s==m_curvelength)
    {
        return m_ps[m_numSamplePoints-1][0];
    }

    for(int i = 0; i < m_numSamplePoints-1; i++)
    {
        if(m_ps[i][1]<s && m_ps[i+1][1]>=s)
        {
            a=i;
            b=i+1;
            break;
        }
    }

    if(a==b)
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

    NekDouble m  = (d2-d1)/(s2-s1);
    NekDouble c = d2- m*s2;

    ASSERTL0(m*s+c==m*s+c,"PS");

    return m*s+c;
}

void CurveMesh::GetSampleFunction()
{
    m_dst.resize(m_numSamplePoints);
    Array<OneD, NekDouble> loc = m_cadcurve->P(m_bounds[0]);

    vector<NekDouble> dsti;
    dsti.resize(3);

    dsti[0]=m_octree->Query(loc);
    dsti[1]=0.0;
    dsti[2]=m_bounds[0];

    m_dst[0] = dsti;

    for(int i = 1; i < m_numSamplePoints; i++)
    {
        dsti[1]=i*ds;
        NekDouble t = m_cadcurve->tAtArcLength(dsti[1]);

        loc = m_cadcurve->P(t);

        dsti[0]=m_octree->Query(loc);
        dsti[2]=t;

        m_dst[i]=dsti;
    }
}

}
}
