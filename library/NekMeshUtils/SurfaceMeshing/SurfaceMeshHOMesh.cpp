////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <list>
#include <algorithm>
#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <NekMeshUtils/Optimisation/Guass-Seidel.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @todo write a optimistaion algorithm for edges on curves, not stricly
 * needed because OCC curves are not distorted, but should be included for
 * completness
 */
void SurfaceMesh::HOSurf()
{
    if(m_mesh->m_verbose)
        cout << endl << "\tHigh-Order Surface meshing" << endl;

    //this bit of code sets up information for the standard edge and face.
    //and a mapping for node ordering for spring optimistaion

    LibUtilities::PointsKey ekey(m_mesh->m_nummode,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;

    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    /*LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                 LibUtilities::eNodalTriFekete);
    Array<OneD, NekDouble> u,v;

    int TotNumPoints = LibUtilities::PointsManager()[pkey]->
                                                    GetTotNumPoints();
    int numInteriorPoints = (m_mesh->m_nummode-3)*(m_mesh->m_nummode-2)/2;

    LibUtilities::PointsManager()[pkey]->GetPoints(u,v);

    DNekMat c (3,3,1.0);
    c(0,0) = u[0];
    c(1,0) = v[0];
    c(2,0) = 1.0;
    c(0,1) = u[1];
    c(1,1) = v[1];
    c(2,1) = 1.0;
    c(0,2) = u[2];
    c(1,2) = v[2];
    c(2,2) = 1.0;
    c.Invert();

    DNekMat p (3,numInteriorPoints,1.0);
    for(int j = 0; j < numInteriorPoints; j++)
    {
        p(0,j) = u[TotNumPoints-numInteriorPoints+j];
        p(1,j) = v[TotNumPoints-numInteriorPoints+j];
        p(2,j) = 1.0;
    }

    map<pair<int,int>, int> nodeorder;
    pair<int, int> id;
    id.first = 1;
    id.second  = 1;
    nodeorder[id] = 0;
    id.first = m_mesh->m_nummode;
    id.second  = 1;
    nodeorder[id] = 1;
    id.first = 1;
    id.second  = m_mesh->m_nummode;
    nodeorder[id] = 2;

    for(int i = 1; i < m_mesh->m_nummode -1; i++)
    {
        id.second = 1;
        id.first  = i+1;
        nodeorder[id] = i+2;
    }
    for(int i = 0; i < m_mesh->m_nummode-2; i++)
    {
        id.first = m_mesh->m_nummode-1 - i;
        id.second = 2 + i;
        nodeorder[id] = m_mesh->m_nummode + 1 + i;
    }
    for(int i = 0; i < m_mesh->m_nummode-2; i++)
    {
        id.first = 1;
        id.second = m_mesh->m_nummode-1 - i;
        nodeorder[id] = m_mesh->m_nummode + m_mesh->m_nummode-1 + i;
    }
    int i = 2;
    int j = 2;
    int limit = m_mesh->m_nummode-1 - 1;
    for(int k = 0; k < numInteriorPoints; k++)
    {
        id.first = i;
        id.second = j;
        nodeorder[id] = 3*(m_mesh->m_nummode-1) + k;
        i++;
        if(i > limit)
        {
            limit--;
            j++;
            i=2;
        }
    }

    map<pair<int,int>, vector<NekDouble> > nodeweight;

    //calculates spring weights from the standard triangle
    for(j = 2; j <= m_mesh->m_nummode-1-1; j++)
    {
        for(i = 2; i <= m_mesh->m_nummode -j; i++)
        {
            NekDouble Zu,Zv,Z1u,Z1v,Z2u,Z2v,Z3u,Z3v;

            Z1u = 0.0; Z1v = 0.0;
            id.first = i-1; id.second = j;
            Z1u += 0.5*u[nodeorder[id]]; Z1v += 0.5*v[nodeorder[id]];
            id.first = i; id.second = j-1;
            Z1u += 0.5*u[nodeorder[id]]; Z1v += 0.5*v[nodeorder[id]];

            Z2u = 0.0; Z2v = 0.0;
            id.first = i+1; id.second = j-1;
            Z2u += 0.5*u[nodeorder[id]]; Z2v += 0.5*v[nodeorder[id]];
            id.first = i+1; id.second = j;
            Z2u += 0.5*u[nodeorder[id]]; Z2v += 0.5*v[nodeorder[id]];

            Z3u = 0.0; Z3v = 0.0;
            id.first = i; id.second = j+1;
            Z3u += 0.5*u[nodeorder[id]]; Z3v += 0.5*v[nodeorder[id]];
            id.first = i-1; id.second = j+1;
            Z3u += 0.5*u[nodeorder[id]]; Z3v += 0.5*v[nodeorder[id]];

            id.first = i; id.second = j;
            Zu = u[nodeorder[id]]; Zv = v[nodeorder[id]];

            vector<NekDouble> K(3);
            K[0] = (Zu - Z2u)*(Zv - Z3v) - (Zv - Z2v)*(Zu - Z3u);
            K[1] = -1.0*((Zu - Z1u)*(Zv - Z3v) - (Zv - Z1v)*(Zu - Z3u));
            K[2] = (Zu - Z1u)*(Zv - Z2v) - (Zv - Z1v)*(Zu - Z2u);

            nodeweight[id] = K;
        }
    }*/

    //because edges are listed twice need a method to not repeat over them
    EdgeSet completedEdges;

    //loop over all the faces in the surface mesh, check all three edges for high
    //order info, if nothing high-order the edge.
    //skip edges which are entirely on planar surfaces
    //if all three edges have no high-order information skip the face

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        if(m_mesh->m_verbose)
        {
            //LibUtilities::PrintProgressbar(i,m_mesh->m_element[2].size(),
            //                               "\t\tSurface elements");
        }

        if(m_mesh->m_element[2][i]->GetConf().m_e == LibUtilities::eQuadrilateral)
        {
            //not setup for high-order quads yet
            continue;
        }

        FaceSharedPtr f = m_mesh->m_element[2][i]->GetFaceLink();
        vector<EdgeSharedPtr> surfedges = m_mesh->m_element[2][i]->GetEdgeList();
        int surf = m_mesh->m_element[2][i]->CADSurfId;

        vector<EdgeSharedPtr> edges = f->m_edgeList;
        for(int j = 0; j < edges.size(); j++)
        {
            //test insert the edge into completedEdges
            //if the edge already exists move on
            //if not figure out its high-order information

            EdgeSet::iterator test = completedEdges.find(edges[j]);

            if(test != completedEdges.end())
            {
                continue;
            }

            EdgeSharedPtr e = edges[j];

            //the edges in the element are different to those in the face
            //the cad information is stored in the element edges which are not
            //in the m_mesh->m_edgeSet group.
            //need to link them together and copy the cad information to be
            //able to identify how to make it high-order
            bool foundsurfaceedge = false;
            for(int k = 0; k < surfedges.size(); k++)
            {
                if(surfedges[k] == e)
                {
                    e->onCurve = surfedges[k]->onCurve;
                    e->CADCurveId = surfedges[k]->CADCurveId;
                    e->CADCurve = surfedges[k]->CADCurve;
                    foundsurfaceedge = true;
                }
            }
            ASSERTL0(foundsurfaceedge,"cannot find corresponding surface edge");


            if(e->onCurve)
            {
                int cid = e->CADCurveId;
                CADCurveSharedPtr c = e->CADCurve;
                NekDouble tb = e->m_n1->GetCADCurveInfo(cid);
                NekDouble te = e->m_n2->GetCADCurveInfo(cid);

                //distrobute points along curve as inital guess
                Array<OneD, NekDouble> ti(m_mesh->m_nummode);
                for(int k = 0; k < m_mesh->m_nummode; k++)
                {
                    ti[k] = tb*(1.0 -  gll[k])/2.0 +
                            te*(1.0 +  gll[k])/2.0;
                }
                Array<OneD, NekDouble> x(m_mesh->m_nummode - 2);
                for(int k = 1; k < m_mesh->m_nummode -1; k++)
                {
                    x[k-1] = ti[k];
                }

                DNekMat J, H;
                CurveEdgeJac(ti, gll, c, J, H);

                bool repeat = true;
                while(repeat)
                {
                    NekDouble Norm = 0;
                    for(int k = 0; k < m_mesh->m_nummode - 2; k++)
                    {
                        Norm += J(k,0)*J(k,0);
                    }
                    Norm = sqrt(Norm);

                    if(Norm < 1E-6)
                    {
                        repeat = false;
                        break;
                    }

                    gsOptimise(x, H, J);

                    for(int k = 1; k < m_mesh->m_nummode -1; k++)
                    {
                        ti[k] = x[k-1];
                    }

                    CurveEdgeJac(ti, gll, c, J, H);
                }

                vector<CADSurfSharedPtr> s = c->GetAdjSurf();

                ASSERTL0(s.size() == 2, "Number of common surfs should be 2");

                vector<NodeSharedPtr> honodes(m_mesh->m_nummode-2);
                for(int k = 1; k < m_mesh->m_nummode -1; k++)
                {
                    Array<OneD, NekDouble> loc = c->P(ti[k]);
                    NodeSharedPtr nn = boost::shared_ptr<Node>(new Node(0,loc[0],loc[1],loc[2]));

                    nn->SetCADCurve(cid, c, ti[k]);
                    Array<OneD, NekDouble> uv = s[0]->locuv(loc);
                    nn->SetCADSurf(s[0]->GetId(), s[0], uv);
                    uv = s[1]->locuv(loc);
                    nn->SetCADSurf(s[1]->GetId(), s[1], uv);
                    honodes[k-1] = nn;
                }

                e->m_edgeNodes = honodes;
                e->m_curveType = LibUtilities::eGaussLobattoLegendre;
                completedEdges.insert(e);
            }
            else
            {
                //edge is on surface and needs 2d optimisation
                CADSurfSharedPtr s = m_cad->GetSurf(surf);
                Array<OneD, NekDouble> uvb,uve;
                uvb = e->m_n1->GetCADSurfInfo(surf);
                uve = e->m_n2->GetCADSurfInfo(surf);
                Array<OneD, Array<OneD, NekDouble> > uvi(m_mesh->m_nummode);
                for(int k = 0; k < m_mesh->m_nummode; k++)
                {
                        Array<OneD, NekDouble> uv(2);
                        uv[0] = uvb[0]*(1.0 - gll[k])/2.0 + uve[0]*(1.0 + gll[k])/2.0;
                        uv[1] = uvb[1]*(1.0 - gll[k])/2.0 + uve[1]*(1.0 + gll[k])/2.0;
                        uvi[k] = uv;
                }

                DNekMat J, H;
                FaceEdgeJac(uvi, gll, s, J, H);
                Array<OneD, NekDouble> x(2*(m_mesh->m_nummode - 2));
                for(int k = 1; k < m_mesh->m_nummode -1; k++)
                {
                    x[(k-1)*2+0] = uvi[k][0];
                    x[(k-1)*2+1] = uvi[k][1];
                }

                bool repeat = true;
                bool hit = false;
                int ct = 0;
                /*while(repeat)
                {
                    NekDouble Norm = 0;
                    for(int k = 0; k < 2.0*(m_mesh->m_nummode - 2); k++)
                    {
                        Norm += J(k,0)*J(k,0);
                    }
                    Norm = sqrt(Norm);

                    if(Norm < 1E-6)
                    {
                        repeat = false;
                        break;
                    }
                    else if(!(Norm < 1E-6) && !hit)
                    {
                        hit = true;
                        cout << s->GetId() << " Norm: ";
                    }

                    cout << Norm << " ";

                    gsOptimise(x, H, J);

                    for(int k = 1; k < m_mesh->m_nummode - 1; k++)
                    {
                        uvi[k][0] = x[(k-1)*2+0];
                        uvi[k][1] = x[(k-1)*2+1];
                    }

                    FaceEdgeJac(uvi, gll, s, J, H);

                    ct++;
                    if(ct > 15)
                    {
                        cout << "exceeded iteration" << endl;
                        break;
                    }
                }*/
                if(hit)
                    cout << endl;

                vector<NodeSharedPtr> honodes(m_mesh->m_nummode-2);
                for(int k = 1; k < m_mesh->m_nummode -1; k++)
                {
                    Array<OneD, NekDouble> loc;
                    loc = s->P(uvi[k]);
                    NodeSharedPtr nn = boost::shared_ptr<Node>(new Node(0,loc[0],loc[1],loc[2]));
                    nn->SetCADSurf(s->GetId(), s, uvi[k]);
                    honodes[k-1] = nn;
                }

                e->m_edgeNodes = honodes;
                e->m_curveType = LibUtilities::eGaussLobattoLegendre;
                completedEdges.insert(e);
            }


        }

        /*
        vector<NodeSharedPtr> vertices = f->m_vertexList;
        Array<OneD, NekDouble> uv1,uv2,uv3;
        uv1 = vertices[0]->GetCADSurfInfo(surf);
        uv2 = vertices[1]->GetCADSurfInfo(surf);
        uv3 = vertices[2]->GetCADSurfInfo(surf);

        DNekMat a (3,3,1.0);
        a(0,0) = uv1[0];
        a(1,0) = uv1[1];
        a(2,0) = 1.0;
        a(0,1) = uv2[0];
        a(1,1) = uv2[1];
        a(2,1) = 1.0;
        a(0,2) = uv3[0];
        a(1,2) = uv3[1];
        a(2,2) = 1.0;

        DNekMat M = a*c;
        DNekMat result = M*p;

        CADSurfSharedPtr s = m_cad->GetSurf(surf);
        vector<NodeSharedPtr> honodes(numInteriorPoints);
        for(int i = 0; i < numInteriorPoints; i++)
        {
            Array<OneD, NekDouble> loc;
            Array<OneD, NekDouble> uv(2);
            uv[0] = result(0,i);
            uv[1] = result(1,i);
            loc = s->P(uv);
            NodeSharedPtr nn = boost::shared_ptr<Node>(new Node(0,loc[0],loc[1],loc[2]));
            nn->SetCADSurf(surf, s, uv);
            honodes[i] = nn;
        }

        //construct a vector of all the uv coords of the triangle in nektar order to form boundary conditions
        vector<Array<OneD, NekDouble> > uvList;
        uvList.push_back(uv1); uvList.push_back(uv2); uvList.push_back(uv3);
        for(int j = 0; j < 3; j++)
        {
            vector<NodeSharedPtr> hon = egs[j]->m_edgeNodes;
            if(egs[j]->m_n1 != vertices[j])
                reverse(hon.begin(),hon.end());

            for(int k = 0; k < hon.size(); k++)
            {
                uvList.push_back(hon[k]->GetCADSurfInfo(surf));
            }
        }
        for(int j = 0; j < honodes.size(); j++)
        {
            uvList.push_back(honodes[j]->GetCADSurfInfo(surf));
        }

        bool repeatoverallnodes = true;

        NekDouble tol = 1E-12;

        while(repeatoverallnodes)
        {
            int converged = 0;

            Array<OneD, NekDouble> uvi;

            Array<OneD, NekDouble> bounds = s->GetBounds();
            int node;
            int hocnt = 0;

            for(int j = 2; j <= m_mesh->m_nummode-1 -1; j++)
            {
                for(int i = 2; i <= m_mesh->m_nummode -j; i++)
                {
                    vector<Array<OneD, NekDouble> > bcs;
                    id.first = i-1;
                    id.second = j;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j-1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i+1;
                    id.second = j-1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i+1;
                    id.second = j;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j+1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i-1;
                    id.second = j+1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j;
                    node = nodeorder[id];
                    uvi = uvList[node];
                    vector<NekDouble> weights = nodeweight[id];

                    bool valid;
                    Array<OneD, NekDouble> df = FaceGrad(uvi[0],uvi[1],bcs,
                                                         weights,surf,
                                                         valid);
                    if(!valid)
                    {
                        converged++;
                        continue;
                    }

                    NekDouble a,b;
                    Find1DBounds(a,b,uvi,df,bounds);

                    NekDouble fxi = FaceF(uvi[0],uvi[1],bcs,weights,surf,valid);
                    NekDouble fx= fxi;

                    NekDouble xmin = BrentOpti(a,0,b,fx,tol,surf,
                                               uvi,df,bounds,bcs,weights,
                                               &SurfaceMesh::FaceF);

                    if(fabs(fx - fxi) < tol)
                    {
                        converged++;
                    }
                    else
                    {
                        uvi[0]+=xmin*df[0]; uvi[1]+=xmin*df[1];
                        Array<OneD, NekDouble> loc = s->P(uvi);
                        uvList[node] = uvi;
                        honodes[hocnt]->Move(loc,surf,uvi);
                    }

                    hocnt++;
                }
            }
            if(converged == honodes.size())
            {
                repeatoverallnodes = false;
            }
        }

        f->m_faceNodes = honodes;
        f->m_curveType = LibUtilities::eNodalTriFekete;*/
    }

    //if(m_mesh->m_verbose)
        //cout << endl;
}

}
}
