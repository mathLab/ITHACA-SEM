////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/SurfaceMeshing.h>
#include <asa047.hpp>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

map<int, MeshNodeSharedPtr> GlobalNodes;
LibUtilities::CADSurfSharedPtr GlobalCad;
int sn, en, m; //start node end node

    void SurfaceMeshing::Mesh()
    {
        if(m_verbose)
            cout << endl << "Surface meshing" << endl;
        if(m_verbose)
            cout << "\tCurve meshing..." << endl;
        for(int i = 1; i <= m_cad->GetNumCurve(); i++)
        {
            if(m_verbose)
            {
                int pos = 70*i/m_cad->GetNumCurve();
                cout << "\t[";
                for (int j = 0; j < 70; ++j) {
                    if (j < pos) cout << "=";
                    else if (j == pos) cout << ">";
                    else cout << " ";
                }
                cout << "] " << int(float(pos)/(70-1)*100)<< " %\r";
                cout.flush();
            }

            m_curvemeshes[i] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(
                    m_verbose, i, m_cad->GetCurve(i), m_octree);

            m_curvemeshes[i]->Mesh(Nodes, Edges);

        }

        if(m_verbose)
        {
            cout << endl << "\tCurve mesh stats:" << endl;
            for(int i = 1; i <= m_cad->GetNumCurve(); i++)
            {
                cout << "\t\tCurve: " << i;
                m_curvemeshes[i]->Report();
            }
        }

        if(m_verbose)
            cout <<endl << "\tSurface meshing..." << endl;

        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            if(m_verbose)
            {
                int pos = 70*i/m_cad->GetNumSurf();
                cout << "\t[";
                for (int j = 0; j < 70; ++j) {
                    if (j < pos) cout << "=";
                    else if (j == pos) cout << ">";
                    else cout << " ";
                }
                cout << "] " << int(float(pos)/(70-1)*100)<< " %\r";
                cout.flush();
            }
            m_surfacemeshes[i] =
                MemoryManager<SurfaceMesh>::AllocateSharedPtr(i,m_verbose,
                    m_cad->GetSurf(i), m_octree,
                    m_curvemeshes,m_order);

            m_surfacemeshes[i]->Mesh(Nodes,Edges,Tris);

        }

        if(m_verbose)
        {
            cout << endl << "\tSurface mesh stats:" << endl;
            for(int i = 1; i <= m_cad->GetNumSurf(); i++)
            {
                cout << "\t\tSurface: " << i;
                m_surfacemeshes[i]->Report();
            }
        }

        if(m_verbose)
        {
            cout << endl << "\tSurface mesh statistics" << endl;
            cout << "\t\tNodes: " << Nodes.size() << endl;
            cout << "\t\tEdges: " << Edges.size() << endl;
            cout << "\t\tTriangles " << Tris.size() << endl;
            cout << "\t\tEuler-Poincaré characteristic: " << Nodes.size()-
                                                             Edges.size()+
                                                              Tris.size() << endl;
        }

        Validate();

        Optimise();

        nodeinlinearmesh = Nodes.size();
    }

    double EnergyEval(double x[])
    {
        NekDouble dz = 2.0/(m/2);

        Array<OneD, NekDouble> loca,locb;
        Array<OneD, NekDouble> uv(2);

        for(int i = 0; i < m; i++)
        {
            //cout << x[i] << " ";
        }


        double p = 0.0;

        loca = GlobalNodes[sn]->GetLoc();
        uv[0] = x[0]; uv[1] = x[1];
        locb = GlobalCad->P(uv);
        p += 1.0/dz*((loca[0]-locb[0])*(loca[0]-locb[0]) +
                            (loca[1]-locb[1])*(loca[1]-locb[1]) +
                            (loca[2]-locb[2])*(loca[2]-locb[2]) );
        int i;
        for(i = 0; i < m/2 - 1; i++)
        {
            uv[0] = x[i*2+0]; uv[1] = x[i*2+1];
            loca = GlobalCad->P(uv);
            uv[0] = x[(i+1)*2+0]; uv[1] = x[(i+1)*2+1];
            locb = GlobalCad->P(uv);
            p += 1.0/dz*((loca[0]-locb[0])*(loca[0]-locb[0]) +
                                (loca[1]-locb[1])*(loca[1]-locb[1]) +
                                (loca[2]-locb[2])*(loca[2]-locb[2]) );
        }
        uv[0] = x[i*2+0]; uv[1] = x[i*2+1];
        loca = GlobalCad->P(uv);
        locb = GlobalNodes[en]->GetLoc();
        p += 1.0/dz*((loca[0]-locb[0])*(loca[0]-locb[0]) +
                            (loca[1]-locb[1])*(loca[1]-locb[1]) +
                            (loca[2]-locb[2])*(loca[2]-locb[2]) );

        //cout << "\t\t" << p << endl;
        return p;
    }

    void SurfaceMeshing::HOSurf()
    {
        if(m_verbose)
            cout << endl << "\tHigh-Order Surface meshing" << endl;

        if(m_verbose)
            cout << "\t\tEdges..." << endl;

        map<int, MeshEdgeSharedPtr>::iterator eit;\
        int counter = 0;

        for(eit = Edges.begin(); eit != Edges.end(); eit++)
        {
            if(m_verbose)
            {
                int pos = 70*counter/Edges.size();
                cout << "\t\t[";
                for (int j = 0; j < 70; ++j) {
                    if (j < pos) cout << "=";
                    else if (j == pos) cout << ">";
                    else cout << " ";
                }
                cout << "] " << int(float(pos)/(70-1)*100)<< " %\r";
                cout.flush();
            }
            counter++;

            MeshEdgeSharedPtr e = eit->second;
            Array<OneD, int> n = e->GetN();

            if(e->GetCurve() != -1)
            {
                LibUtilities::CADCurveSharedPtr c = m_cad->GetCurve(
                            e->GetCurve());
                //edge is on curve and needs hoing that way
                NekDouble tb = Nodes[n[0]]->GetC(e->GetCurve());
                NekDouble te = Nodes[n[1]]->GetC(e->GetCurve());

                NekDouble dz = 2.0/m_order;

                Array<OneD, NekDouble> ti(m_order+1);

                for(int i = 0; i < m_order+1; i++)
                {
                    ti[i] = tb*((1.0 - (-1.0 + dz*i))/2.0) +
                                   te*((1.0 + (-1.0 + dz*i))/2.0);
                }

                vector<int> Surfs = c->GetAdjSurf();

                ASSERTL0(Surfs.size() == 2, "Number of common surfs should be 2");

                vector<int> honodes(m_order-1);

                LibUtilities::CADSurfSharedPtr s1,s2;
                s1 = m_cad->GetSurf(Surfs[0]);
                s2 = m_cad->GetSurf(Surfs[1]);

                for(int i = 1; i < m_order +1 -1; i++)
                {
                    Array<OneD, NekDouble> loc = c->P(ti[i]);
                    MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                            AllocateSharedPtr(Nodes.size(),loc[0],
                                              loc[1],loc[2]);
                    nn->SetCurve(c->GetID(),ti[i]);
                    Array<OneD, NekDouble> uv = s1->locuv(loc);
                    nn->SetSurf(Surfs[0],uv);
                    uv = s2->locuv(loc);
                    nn->SetSurf(Surfs[1],uv);
                    honodes[i-1] = Nodes.size();
                    Nodes[Nodes.size()] = nn;
                }

                e->SetHONodes(honodes);

                //this needs a 1d optimisation step

            }
            else
            {
                //edge is on surface and needs 2d optimisation
                LibUtilities::CADSurfSharedPtr s = m_cad->GetSurf(e->GetSurf());
                Array<OneD, NekDouble> uvb,uve;
                uvb = Nodes[n[0]]->GetS(e->GetSurf());
                uve = Nodes[n[1]]->GetS(e->GetSurf());

                vector<int> honodes(m_order-1);
                for(int i = 1; i < m_order+1 -1; i++)
                {
                    Array<OneD, NekDouble> loc;
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = uvb[0]+i*(uve[0]-uvb[0])/m_order;
                    uv[1] = uvb[1]+i*(uve[1]-uvb[1])/m_order;
                    loc = s->P(uv);
                    MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                            AllocateSharedPtr(Nodes.size(),loc[0],
                                                loc[1],loc[2]);
                    nn->SetSurf(e->GetSurf(),uv);
                    honodes[i-1] = Nodes.size();
                    Nodes[Nodes.size()] = nn;

                }

                if(s->IsPlane() == true)
                {
                    e->SetHONodes(honodes);
                    continue; //optimimum points on plane are linear
                }
                e->SetHONodes(honodes);
                continue;

                GlobalNodes = Nodes; sn = n[0]; en = n[1]; GlobalCad = m_cad->GetSurf(e->GetSurf());
                m = honodes.size()*2;


                double *start;
                start = new double[honodes.size()*2];
                double *xmin;
                xmin = new double[honodes.size()*2];
                double *step;
                step = new double[honodes.size()*2];

                for(int i = 0; i < honodes.size(); i++)
                {
                    Array<OneD, NekDouble> uv = Nodes[honodes[i]]->GetS(e->GetSurf());
                    start[i*2+0] = uv[0];
                    start[i*2+1] = uv[1];
                    step[i*2+0] = (max(uve[0],uvb[0]) - min(uve[0],uvb[0]));
                    step[i*2+1] = (max(uve[1],uvb[1]) - min(uve[1],uvb[1]));
                }
                double ynew = EnergyEval(start);
                //cout << ynew << endl;
                int icount, ifault, numres;
                nelmin(EnergyEval, honodes.size()*2, start, xmin, &ynew, 1E-5, step,
                       10, 500, &icount, &numres, &ifault);

                int repeats = 0;
                while(ifault == 2)
                {
                    start = xmin;
                    nelmin(EnergyEval, honodes.size()*2, start, xmin, &ynew, 1E-5, step,
                           10, 100, &icount, &numres, &ifault);
                    repeats+=1;
                }

                for(int i = 0; i < honodes.size(); i++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = xmin[i*2+0]; uv[1] = xmin[i*2+1];
                    Array<OneD, NekDouble> l = m_cad->GetSurf(e->GetSurf())->P(uv);
                    Nodes[honodes[i]]->Move(l,uv);
                }

                e->SetHONodes(honodes);
            }
        }

        if(m_verbose)
            cout << endl << "\t\tFaces..." << endl;

        map<int, MeshTriSharedPtr>::iterator trit;

        LibUtilities::PointsKey pkey(m_order+1,
                                     LibUtilities::eNodalTriEvenlySpaced);
        Array<OneD, NekDouble> u,v;

        int TotNumPoints = LibUtilities::PointsManager()[pkey]->
                                                        GetTotNumPoints();
        int numInteriorPoints = (m_order-2)*(m_order-1)/2;

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

        for(trit = Tris.begin(); trit != Tris.end(); trit++)
        {

            Array<OneD, int> n = trit->second->GetN();

            Array<OneD, NekDouble> uv1,uv2,uv3;
            uv1 = Nodes[n[0]]->GetS(trit->second->Getcid());
            uv2 = Nodes[n[1]]->GetS(trit->second->Getcid());
            uv3 = Nodes[n[2]]->GetS(trit->second->Getcid());

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

            vector<int> honodes(numInteriorPoints);
            for(int i = 0; i < numInteriorPoints; i++)
            {
                Array<OneD, NekDouble> loc;
                Array<OneD, NekDouble> uv(2);
                uv[0] = result(0,i);
                uv[1] = result(1,i);
                loc = m_cad->GetSurf(trit->second->Getcid())->P(uv);
                MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                        AllocateSharedPtr(Nodes.size(),loc[0],
                                            loc[1],loc[2]);
                nn->SetSurf(trit->second->Getcid(),uv);
                honodes[i] = Nodes.size();
                Nodes[Nodes.size()] = nn;

            }

            trit->second->SetHONodes(honodes);

        }

    }

    void SurfaceMeshing::Optimise()
    {
        if(m_verbose)
            cout << endl << "\tOptimising linear mesh" << endl;

        for(int q = 0; q <2; q++)
        {
            if(m_verbose)
                cout << "\t\t Edge swap run: " << q+1 << endl;

            map<int, MeshEdgeSharedPtr>::iterator it;
            for(it = Edges.begin(); it != Edges.end(); it++)
            {
                if(it->second->GetCurve() != -1)
                    continue;

                vector<int> t = it->second->GetTri();
                ASSERTL0(t.size() == 2, "each edge should have two tri");
                Array<OneD, int> n = it->second->GetN();
                Array<OneD, int> nt = Tris[t[0]]->GetN();
                int nodecheck = 0;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 2; j++)
                    {
                        if(n[j] == nt[i])
                        {
                            nodecheck++;
                        }
                    }
                }
                ASSERTL0(nodecheck == 2, "edge and tri should have 2 n in commom");

                int A,B,C,D;
                if(nt[0] != n[0] && nt[0] != n[1])
                {
                    C = nt[0];
                    B = nt[1];
                    A = nt[2];
                }
                else if(nt[1] != n[0] && nt[1] != n[1])
                {
                    C = nt[1];
                    B = nt[2];
                    A = nt[0];
                }
                else if(nt[2] != n[0] && nt[2] != n[1])
                {
                    C = nt[2];
                    B = nt[0];
                    A = nt[1];
                }

                nt = Tris[t[1]]->GetN();
                nodecheck = 0;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 2; j++)
                    {
                        if(n[j] == nt[i])
                        {
                            nodecheck++;
                        }
                    }
                }
                ASSERTL0(nodecheck == 2, "edge and tri should have 2 n in commom");

                if(nt[0] != n[0] && nt[0] != n[1])
                {
                    D = nt[0];
                }
                else if(nt[1] != n[0] && nt[1] != n[1])
                {
                    D = nt[1];
                }
                else if(nt[2] != n[0] && nt[2] != n[1])
                {
                    D = nt[2];
                }

                //determine signed area of alternate config
                Array<OneD, NekDouble> ai,bi,ci,di;
                ai = Nodes[A]->GetS(it->second->GetSurf());
                bi = Nodes[B]->GetS(it->second->GetSurf());
                ci = Nodes[C]->GetS(it->second->GetSurf());
                di = Nodes[D]->GetS(it->second->GetSurf());

                NekDouble CDA, CBD;

                CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

                CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                      (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

                if(CDA < 0.0 || CBD < 0.0)
                    continue;

                int nodedefectbefore = 0;
                nodedefectbefore += Nodes[A]->GetEdges().size() > 6 ?
                                    Nodes[A]->GetEdges().size() - 6 :
                                    6 - Nodes[A]->GetEdges().size();
                nodedefectbefore += Nodes[B]->GetEdges().size() > 6 ?
                                    Nodes[B]->GetEdges().size() - 6 :
                                    6 - Nodes[B]->GetEdges().size();
                nodedefectbefore += Nodes[C]->GetEdges().size() > 6 ?
                                    Nodes[C]->GetEdges().size() - 6 :
                                    6 - Nodes[C]->GetEdges().size();
                nodedefectbefore += Nodes[D]->GetEdges().size() > 6 ?
                                    Nodes[D]->GetEdges().size() - 6 :
                                    6 - Nodes[D]->GetEdges().size();



                int nodedefectafter = 0;
                nodedefectafter  += Nodes[A]->GetEdges().size() - 1 > 6 ?
                                    Nodes[A]->GetEdges().size() - 1 - 6 :
                                    6 - (Nodes[A]->GetEdges().size() - 1);
                nodedefectafter  += Nodes[B]->GetEdges().size() - 1 > 6 ?
                                    Nodes[B]->GetEdges().size() - 1 - 6 :
                                    6 - (Nodes[B]->GetEdges().size() - 1);
                nodedefectafter  += Nodes[C]->GetEdges().size() + 1 > 6 ?
                                    Nodes[C]->GetEdges().size() + 1 - 6 :
                                    6 - (Nodes[C]->GetEdges().size() + 1);
                nodedefectafter  += Nodes[D]->GetEdges().size() + 1 > 6 ?
                                    Nodes[D]->GetEdges().size() + 1 - 6 :
                                    6 - (Nodes[D]->GetEdges().size() + 1);

                if(nodedefectafter < nodedefectbefore)
                {
                    Edges[it->first]->Swap(C,D);
                    Nodes[C]->SetEdge(it->first);
                    Nodes[D]->SetEdge(it->first);
                    Nodes[A]->RemoveEdge(it->first);
                    Nodes[B]->RemoveEdge(it->first);
                    Nodes[C]->SetTri(t[1]);
                    Nodes[B]->RemoveTri(t[1]);
                    Nodes[D]->SetTri(t[0]);
                    Nodes[A]->RemoveTri(t[0]);


                    int e1, e2, e3;
                    e1 = Nodes[C]->EdgeInCommon(Nodes[B]);
                    e2 = Nodes[B]->EdgeInCommon(Nodes[D]);
                    Edges[e2]->RemoveTri(t[1]);
                    Edges[e2]->SetTri(t[0]);
                    e3 = it->first;
                    ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                    Tris[t[0]]->Swap(C,B,D);
                    Tris[t[0]]->ResetEdges(e1, e2, e3);

                    Tris[t[1]]->Swap(C,D,A);
                    e1 = it->first;
                    e2 = Nodes[D]->EdgeInCommon(Nodes[A]);
                    e3 = Nodes[A]->EdgeInCommon(Nodes[C]);
                    Edges[e3]->RemoveTri(t[0]);
                    Edges[e3]->SetTri(t[1]);
                    ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                    Tris[t[1]]->ResetEdges(e1, e2 ,e3);


                }
            }
        }

        for(int q = 0; q <2; q++)
        {
            if(m_verbose)
                cout << "\t\t Elastic relaxation run: " << q+1 << endl;

            map<int, MeshEdgeSharedPtr>::iterator it;
            for(it = Edges.begin(); it != Edges.end(); it++)
            {
                if(it->second->GetCurve() != -1)
                    continue;

                vector<int> t = it->second->GetTri();
                ASSERTL0(t.size() == 2, "each edge should have two tri");
                Array<OneD, int> n = it->second->GetN();
                Array<OneD, int> nt = Tris[t[0]]->GetN();

                int A,B,C,D;
                if(nt[0] != n[0] && nt[0] != n[1])
                {
                    C = nt[0];
                    B = nt[1];
                    A = nt[2];
                }
                else if(nt[1] != n[0] && nt[1] != n[1])
                {
                    C = nt[1];
                    B = nt[2];
                    A = nt[0];
                }
                else if(nt[2] != n[0] && nt[2] != n[1])
                {
                    C = nt[2];
                    B = nt[0];
                    A = nt[1];
                }

                nt = Tris[t[1]]->GetN();

                if(nt[0] != n[0] && nt[0] != n[1])
                {
                    D = nt[0];
                }
                else if(nt[1] != n[0] && nt[1] != n[1])
                {
                    D = nt[1];
                }
                else if(nt[2] != n[0] && nt[2] != n[1])
                {
                    D = nt[2];
                }

                //determine signed area of alternate config
                Array<OneD, NekDouble> ai,bi,ci,di;
                ai = Nodes[A]->GetS(it->second->GetSurf());
                bi = Nodes[B]->GetS(it->second->GetSurf());
                ci = Nodes[C]->GetS(it->second->GetSurf());
                di = Nodes[D]->GetS(it->second->GetSurf());

                NekDouble CDA, CBD;

                CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

                CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                      (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

                if(CDA < 0.0 || CBD < 0.0)
                    continue;

                NekDouble minangleb = Nodes[C]->Angle(Nodes[A],Nodes[B]);
                minangleb = min(minangleb,Nodes[B]->Angle(Nodes[C],Nodes[A]));
                minangleb = min(minangleb,Nodes[A]->Angle(Nodes[B],Nodes[C]));

                minangleb = min(minangleb,Nodes[B]->Angle(Nodes[A],Nodes[D]));
                minangleb = min(minangleb,Nodes[D]->Angle(Nodes[B],Nodes[A]));
                minangleb = min(minangleb,Nodes[A]->Angle(Nodes[D],Nodes[B]));

                NekDouble minanglea = Nodes[C]->Angle(Nodes[D],Nodes[B]);
                minanglea = min(minanglea,Nodes[B]->Angle(Nodes[C],Nodes[D]));
                minanglea = min(minanglea,Nodes[D]->Angle(Nodes[B],Nodes[C]));

                minanglea = min(minanglea,Nodes[C]->Angle(Nodes[A],Nodes[D]));
                minanglea = min(minanglea,Nodes[D]->Angle(Nodes[C],Nodes[A]));
                minanglea = min(minanglea,Nodes[A]->Angle(Nodes[D],Nodes[C]));

                if(minanglea > minangleb)
                {
                    Edges[it->first]->Swap(C,D);
                    Nodes[C]->SetEdge(it->first);
                    Nodes[D]->SetEdge(it->first);
                    Nodes[A]->RemoveEdge(it->first);
                    Nodes[B]->RemoveEdge(it->first);
                    Nodes[C]->SetTri(t[1]);
                    Nodes[B]->RemoveTri(t[1]);
                    Nodes[D]->SetTri(t[0]);
                    Nodes[A]->RemoveTri(t[0]);


                    int e1, e2, e3;
                    e1 = Nodes[C]->EdgeInCommon(Nodes[B]);
                    e2 = Nodes[B]->EdgeInCommon(Nodes[D]);
                    Edges[e2]->RemoveTri(t[1]);
                    Edges[e2]->SetTri(t[0]);
                    e3 = it->first;
                    ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                    Tris[t[0]]->Swap(C,B,D);
                    Tris[t[0]]->ResetEdges(e1, e2, e3);

                    Tris[t[1]]->Swap(C,D,A);
                    e1 = it->first;
                    e2 = Nodes[D]->EdgeInCommon(Nodes[A]);
                    e3 = Nodes[A]->EdgeInCommon(Nodes[C]);
                    Edges[e3]->RemoveTri(t[0]);
                    Edges[e3]->SetTri(t[1]);
                    ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                    Tris[t[1]]->ResetEdges(e1, e2 ,e3);

                }
            }
        }

        for(int q = 0; q < 4; q++)
        {
            map<int, MeshNodeSharedPtr>::iterator it;
            for(it = Nodes.begin(); it!=Nodes.end(); it++)
            {
                if(it->second->IsOnACurve())
                    continue;

                NekDouble d = m_octree->Query(it->second->GetLoc());

                std::map<int, Array<OneD, NekDouble> > surf = it->second->GetSurfMap();
                ASSERTL0(surf.size()==1,"node should be interior and only be on one surface");
                std::map<int, Array<OneD, NekDouble> >::iterator sit = surf.begin();
                int surface = sit->first;
                Array<OneD, NekDouble> uvi = sit->second;

                vector<int> es = it->second->GetEdges();
                vector<int> connodes;
                for(int i = 0; i < es.size(); i++)
                {
                    connodes.push_back(Edges[es[i]]->OtherNode(it->first));
                }

                vector<NekDouble> om;
                for(int i = 0; i < connodes.size(); i++)
                {
                    om.push_back(it->second->Distance(Nodes[connodes[i]]) - d);
                }

                NekDouble u0=0.0,v0=0.0,fu=0.0,dfu=0.0,fv=0.0,dfv=0.0;
                for(int i = 0; i < connodes.size(); i++)
                {
                    Array<OneD, NekDouble> uvj = Nodes[connodes[i]]->GetS(surface);
                    u0+=uvj[0]/connodes.size();
                    v0+=uvj[1]/connodes.size();
                }
                for(int i = 0; i < connodes.size();i++)
                {
                    Array<OneD, NekDouble> uvj = Nodes[connodes[i]]->GetS(surface);
                    NekDouble sqr = sqrt((uvj[0]-u0)*(uvj[0]-u0) +
                                         (uvj[1]-v0)*(uvj[1]-v0));
                    fu+=om[i]*(uvj[0]-u0)/sqr;
                    fv+=om[i]*(uvj[1]-v0)/sqr;
                    dfu+=om[i]*sqr*(2*(uvj[0]-u0)*(u0-uvj[0])+
                                        (uvj[1]-v0)*(uvj[1]-v0));
                    dfv+=om[i]*sqr*(2*(uvj[1]-v0)*(uvi[1]-v0)+
                                        (uvj[0]-u0)*(uvj[0]-u0));
                }
                Array<OneD, NekDouble> uv(2);
                uv[0] = u0-fu/dfu; uv[1] = v0-fv/dfv;
                Array<OneD, NekDouble> l2 = m_cad->GetSurf(surface)->P(uv);
                Nodes[it->first]->Move(l2,uv);

            }
        }
    }

    void SurfaceMeshing::Validate()
    {
        if(m_verbose)
            cout << endl << "\tVerifying surface mesh" << endl;

        if(m_cad->GetEPC() != Nodes.size()-Edges.size()+Tris.size())
        {
            if(m_verbose)
                cout << "\t\tFailed" << endl;
            ASSERTL0(false,"Euler-Poincaré characteristics do not match");
        }

        map<int,int> edgecheck;
        map<int, MeshEdgeSharedPtr>::iterator ite;
        for(ite=Edges.begin(); ite!=Edges.end(); ite++)
        {
            edgecheck[ite->first] = 0;
        }

        map<int, MeshTriSharedPtr>::iterator it;
        for(it=Tris.begin(); it!=Tris.end(); it++)
        {
            Array<OneD,int> ed = it->second->GetE();
            for(int i = 0; i < 3; i++)
            {
                edgecheck[ed[i]]++;
            }
        }

        map<int,int>::iterator ch;
        for(ch=edgecheck.begin(); ch!=edgecheck.end(); ch++)
        {
            if(ch->second != 2)
            {
                if(m_verbose)
                    cout << "\t\tFailed" << endl;
                ASSERTL0(false,"edge not listed twice");
            }
        }
        if(m_verbose)
            cout << "\t\tPassed" << endl;
    }


}
}
