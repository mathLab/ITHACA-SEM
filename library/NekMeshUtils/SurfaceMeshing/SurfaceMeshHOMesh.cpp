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
#include <NekMeshUtils/Optimisation/BGFS-B.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

map<int, vector<int> > nodeToSixAround(int nq)
{
    map<pair<int,int>, int> nodeorder;
    map<int, pair<int, int> > nodeorderRev;
    pair<int, int> id;

    id.first = 0;
    id.second  = 0;
    nodeorder[id] = 0;
    nodeorderRev[0] = id;

    id.first = nq-1;
    id.second  = 0;
    nodeorder[id] = 1;
    nodeorderRev[1] = id;


    id.first = 0;
    id.second  = nq-1;
    nodeorder[id] = 2;
    nodeorderRev[2] = id;

    for(int i = 0; i < nq -2; i++)
    {
        id.second = 0;
        id.first  = i+1;
        nodeorder[id] = i+3;
        nodeorderRev[i+3] = id;
    }
    for(int i = 0; i < nq -2; i++)
    {
        id.first = nq-2 - i;
        id.second = 1 + i;
        nodeorder[id] = nq + 1 + i;
        nodeorderRev[nq + 1 + i] = id;
    }
    for(int i = 0; i < nq-2; i++)
    {
        id.first = 0;
        id.second = nq-2 - i;
        nodeorder[id] = nq + nq-1 + i;
        nodeorderRev[nq + nq-1 + i] = id;
    }

    int i = 1;
    int j = 1;
    int limit = nq-3;
    for(int k = 0; k < (nq-3)*(nq-2)/2; k++)
    {
        id.first = i;
        id.second = j;
        nodeorder[id] = 3*(nq-1) + k;
        nodeorderRev[3*(nq-1) + k] = id;
        i++;
        if(i > limit)
        {
            limit--;
            j++;
            i=1;
        }
    }

    map<int, vector<int> > ret;

    for(int i = (nq+1)*nq/2 - (nq-3)*(nq-2)/2; i < (nq+1)*nq/2; i++)
    {
        vector<int> ids;
        int pr;

        pair<int,int> p = nodeorderRev[i];
        p.first -= 1;
        pr = nodeorder[p];

        ids.push_back(pr);

        p = nodeorderRev[i];
        p.second -=1;
        pr = nodeorder[p];

        ids.push_back(pr);

        p = nodeorderRev[i];
        p.first +=1;
        pr = nodeorder[p];

        ids.push_back(pr);

        p = nodeorderRev[i];
        p.second +=1;
        pr = nodeorder[p];

        ids.push_back(pr);

        p = nodeorderRev[i];
        p.first +=1;
        p.second -=1;
        pr = nodeorder[p];

        ids.push_back(pr);

        p = nodeorderRev[i];
        p.first -=1;
        p.second +=1;
        pr = nodeorder[p];

        ids.push_back(pr);

        ret[i] = ids;
    }

    return ret;
}

map<int, Array<OneD, NekDouble> > weights(map<int, vector<int> > near, Array<OneD, NekDouble> u, Array<OneD, NekDouble> v)
{
    map<int, Array<OneD, NekDouble> > ret;

    //setup map from right angled reference triangle to equilateral reference triangle
    DNekMat A(3,3,1.0);
    A(0,0) = -1.0; A(1,0) = -1.0;
    A(0,1) =  1.0; A(1,1) = -1.0;
    A(0,2) =  0.0; A(1,2) = -1.0+sqrt(3.0);

    DNekMat B(3,3,1.0);
    B(0,0) = -1.0; B(1,0) = -1.0;
    B(0,1) =  1.0; B(1,1) = -1.0;
    B(0,2) = -1.0; B(1,2) =  1.0;

    B.Invert();

    DNekMat M = A*B;

    DNekMat C(3, u.num_elements(), 1.0);
    for(int i = 0; i < u.num_elements(); i++)
    {
        C(0,i) = u[i];
        C(1,i) = v[i];
    }

    DNekMat pts = M*C;

    map<int, vector<int> >::iterator it;
    for(it = near.begin(); it != near.end(); it++)
    {
        Array<OneD, NekDouble> s(6);

        NekDouble zu = pts(0,it->first);
        NekDouble zv = pts(1,it->first);

        s[0] = sqrt((zu-pts(0,it->second[0]))*(zu-pts(0,it->second[0]))+(zv-pts(1,it->second[0]))*(zv-pts(1,it->second[0])));
        s[1] = sqrt((zu-pts(0,it->second[1]))*(zu-pts(0,it->second[1]))+(zv-pts(1,it->second[1]))*(zv-pts(1,it->second[1])));
        s[2] = sqrt((zu-pts(0,it->second[2]))*(zu-pts(0,it->second[2]))+(zv-pts(1,it->second[2]))*(zv-pts(1,it->second[2])));
        s[3] = sqrt((zu-pts(0,it->second[3]))*(zu-pts(0,it->second[3]))+(zv-pts(1,it->second[3]))*(zv-pts(1,it->second[3])));
        s[4] = sqrt((zu-pts(0,it->second[4]))*(zu-pts(0,it->second[4]))+(zv-pts(1,it->second[4]))*(zv-pts(1,it->second[4])));
        s[5] = sqrt((zu-pts(0,it->second[5]))*(zu-pts(0,it->second[5]))+(zv-pts(1,it->second[5]))*(zv-pts(1,it->second[5])));

        ret[it->first] = s;
    }

    return ret;
}

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

    LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                 LibUtilities::eNodalTriFekete);
    Array<OneD, NekDouble> u,v;

    int np = LibUtilities::PointsManager()[pkey]->GetTotNumPoints();
    int nq = m_mesh->m_nummode;
    int ni = (nq-3)*(nq-2)/2;

    LibUtilities::PointsManager()[pkey]->GetPoints(u,v);

    map<int, vector<int> > near = nodeToSixAround(m_mesh->m_nummode);
    map<int, Array<OneD, NekDouble> > z = weights(near, u, v);

    //because edges are listed twice need a method to not repeat over them
    EdgeSet completedEdges;

    //loop over all the faces in the surface mesh, check all three edges for high
    //order info, if nothing high-order the edge.

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_mesh->m_element[2].size(),
                                           "\t\tSurface elements");
        }

        if(m_mesh->m_element[2][i]->GetConf().m_e == LibUtilities::eQuadrilateral)
        {
            //not setup for high-order quads yet
            continue;
        }

        int surf = m_mesh->m_element[2][i]->CADSurfId;
        CADSurfSharedPtr s = m_cad->GetSurf(surf);

        FaceSharedPtr f = m_mesh->m_element[2][i]->GetFaceLink();
        vector<EdgeSharedPtr> surfedges = m_mesh->m_element[2][i]->GetEdgeList();

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

                DNekMat B(nq-2,nq-2,0.0); //approximate hessian (I to start)
                for(int k = 0; k < nq -2; k++)
                {
                    B(k,k) = 1.0;
                }
                DNekMat H(nq-2,nq-2,0.0); //approximate inverse hessian (I to start)
                for(int k = 0; k < nq -2; k++)
                {
                    H(k,k) = 1.0;
                }

                DNekMat J = EdgeGrad(ti, gll, c);

                bool repeat = true;
                int itct = 0;
                bool hit = false;
                while(repeat)
                {
                    NekDouble Norm = 0;
                    for(int k = 0; k < nq - 2; k++)
                    {
                        Norm += J(k,0)*J(k,0);
                    }
                    Norm = sqrt(Norm);

                    if(Norm < 1E-6)
                    {
                        repeat = false;
                        break;
                    }
                    if(itct > 100)
                    {
                        cout << "failed to optimise on curve" << endl;
                        for(int k = 0; k < nq; k++)
                        {
                            ti[k] = tb*(1.0 -  gll[k])/2.0 +
                                    te*(1.0 +  gll[k])/2.0;
                        }
                        exit(-1);
                        break;
                    }
                    itct++;
                    hit= true;
                    cout << "new iter curve" << endl;
                    cout << endl << "Norm " << Norm << endl;

                    BGFSUpdate(EdgeF, EdgeGrad, ti, gll, c, J, B, H); //ti will be updated

                }
                //if(hit) exit(-1);

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
                Array<OneD, Array<OneD, NekDouble> > uvi(nq);
                for(int k = 0; k < nq; k++)
                {
                        Array<OneD, NekDouble> uv(2);
                        uv[0] = uvb[0]*(1.0 - gll[k])/2.0 + uve[0]*(1.0 + gll[k])/2.0;
                        uv[1] = uvb[1]*(1.0 - gll[k])/2.0 + uve[1]*(1.0 + gll[k])/2.0;
                        uvi[k] = uv;
                }

                Array<OneD, NekDouble> bnds = s->GetBounds();
                Array<OneD, NekDouble> all(2*nq);
                for(int k = 0; k < nq; k++)
                {
                    all[k*2+0] = uvi[k][0];
                    all[k*2+1] = uvi[k][1];
                }

                DNekMat B(2*(nq-2),2*(nq-2),0.0); //approximate hessian (I to start)
                for(int k = 0; k < 2*(nq -2); k++)
                {
                    B(k,k) = 1.0;
                }
                DNekMat H(2*(nq-2),2*(nq-2),0.0); //approximate inverse hessian (I to start)
                for(int k = 0; k < 2*(nq -2); k++)
                {
                    H(k,k) = 1.0;
                }

                DNekMat J = EdgeGrad(all, gll, s);

                bool repeat = true;
                int itct = 0;
                bool hit = false;
                while(repeat)
                {
                    NekDouble Norm = 0;
                    for(int k = 0; k < nq - 2; k++)
                    {
                        Norm += J(k,0)*J(k,0);
                    }
                    Norm = sqrt(Norm);

                    if(Norm < 1E-6)
                    {
                        repeat = false;
                        break;
                    }
                    if(itct > 100)
                    {
                        cout << "failed to optimise on edge" << endl;
                        for(int k = 0; k < nq; k++)
                        {
                            Array<OneD, NekDouble> uv(2);
                            uv[0] = uvb[0]*(1.0 - gll[k])/2.0 + uve[0]*(1.0 + gll[k])/2.0;
                            uv[1] = uvb[1]*(1.0 - gll[k])/2.0 + uve[1]*(1.0 + gll[k])/2.0;
                            uvi[k] = uv;
                        }
                        exit(-1);
                        break;
                    }
                    itct++;
                    hit= true;
                    cout << "new iter edge" << endl;
                    cout << endl << "Norm " << Norm << endl;

                    BGFSUpdate(EdgeF, EdgeGrad, all, gll, s, J, B, H); //all will be updated
                }

                //need to put all backinto uv
                for(int k = 0; k < nq; k++)
                {
                    uvi[k][0] = all[k*2+0];
                    uvi[k][1] = all[k*2+1];
                }

                vector<NodeSharedPtr> honodes(nq-2);
                for(int k = 1; k < nq -1; k++)
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

        ASSERTL0(nq <= 5,"not setup for high-orders yet");

        vector<NodeSharedPtr> vertices = f->m_vertexList;

        SpatialDomains::Geometry2DSharedPtr geom = f->GetGeom(3);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        Array<OneD, NekDouble> xc(nq*nq);
        Array<OneD, NekDouble> yc(nq*nq);
        Array<OneD, NekDouble> zc(nq*nq);

        xmap->BwdTrans(coeffs0,xc);
        xmap->BwdTrans(coeffs1,yc);
        xmap->BwdTrans(coeffs2,zc);


        //build an array of all uvs
        Array<OneD, Array<OneD, NekDouble> > uvi(np);
        int ctr = 0;
        for(int j = 0; j < vertices.size(); j++)
        {
            uvi[ctr++] = vertices[j]->GetCADSurfInfo(surf);
        }
        for(int j = 0; j < edges.size(); j++)
        {
            vector<NodeSharedPtr> ns = edges[j]->m_edgeNodes;
            if(!(edges[j]->m_n1 == vertices[j]))
            {
                reverse(ns.begin(),ns.end());
            }
            for(int k = 0; k < ns.size(); k++)
            {
                uvi[ctr++] = ns[k]->GetCADSurfInfo(surf);
            }
        }
        for(int j = np-ni; j < np; j++)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = u[j];
            xp[1] = v[j];

            Array<OneD, NekDouble> xyz(3);
            xyz[0] = xmap->PhysEvaluate(xp, xc);
            xyz[1] = xmap->PhysEvaluate(xp, yc);
            xyz[2] = xmap->PhysEvaluate(xp, zc);

            Array<OneD, NekDouble> uv(2);
            s->ProjectTo(xyz,uv);
            uvi[ctr++] = uv;
        }

        DNekMat H, J;

        Array<OneD, NekDouble> bnds = s->GetBounds();

        bool repeat = true;
        NekDouble alpha = 1.0;
        int itct = 0;
        /*while(repeat)
        {
            if(itct > 10000)
            {
                cout << "failed to converge on face" << endl;
                break;
            }
            itct++;

            int converged = 0;

            for(int j = np - ni; j < np; j++)
            {
                FaceFaceJac(j, uvi, z, near, s, J, H);

                NekDouble Norm = 0;
                for(int k = 0; k < 2.0; k++)
                {
                    Norm += J(k,0)*J(k,0);
                }
                Norm = sqrt(Norm);

                if(Norm < 1E-8)
                {
                    converged++;
                }

                bool repeatNode = true;
                bool fail = false;
                int noitct = 0;
                while(repeatNode)
                {
                    if(noitct > 100)
                    {
                        break;
                    }
                    noitct++;

                    Array<OneD, NekDouble> dx = gsOptimise(alpha, uvi[j], H, J);

                    if(uvi[j][0] + dx[0] < bnds[0] ||
                       uvi[j][0] + dx[0] > bnds[1] ||
                       uvi[j][1] + dx[1] < bnds[2] ||
                       uvi[j][1] + dx[1] > bnds[3])
                    {
                        fail = true;
                        Array<OneD, NekDouble> ua(2,0.0);
                        for(int l = 0; l < 6; l++)
                        {
                            ua[0] += uvi[near[j][l]][0] / 6.0;
                            ua[1] += uvi[near[j][l]][1] / 6.0;
                        }
                        uvi[j][0] = ua[0];
                        uvi[j][1] = ua[1];
                    }
                    else
                    {
                        uvi[j][0] += dx[0];
                        uvi[j][1] += dx[1];
                    }

                    FaceFaceJac(j, uvi, z, near, s, J, H);

                    Norm = 0;
                    for(int k = 0; k < 2.0; k++)
                    {
                        Norm += J(k,0)*J(k,0);
                    }
                    Norm = sqrt(Norm);

                    if(Norm < 1E-8)
                    {
                        repeatNode = false;
                    }
                }
            }

            if(converged == ni)
            {
                repeat = false;
            }
        }

        vector<NodeSharedPtr> honodes;
        for(int j = np - ni; j < np; j++)
        {
            Array<OneD, NekDouble> loc;
            loc = s->P(uvi[j]);
            NodeSharedPtr nn = boost::shared_ptr<Node>(new Node(0,loc[0],loc[1],loc[2]));
            nn->SetCADSurf(surf, s, uvi[j]);
            honodes.push_back(nn);
        }

        f->m_faceNodes = honodes;
        f->m_curveType = LibUtilities::eNodalTriFekete;*/
    }

    if(m_mesh->m_verbose)
        cout << endl;
}

}
}
