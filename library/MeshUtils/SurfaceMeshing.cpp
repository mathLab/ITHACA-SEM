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

#include <string>
#include <fstream>
#include <sstream>

#include <MeshUtils/SurfaceMeshing.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

    void SurfaceMeshing::Mesh()
    {
        for(int i = 1; i <= m_cad->GetNumCurve(); i++)
        {
            if(m_verbose)
                cout << endl << "Meshing Curve: " << i << endl;

            m_curvemeshes[i] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(
                    m_verbose, i, m_cad->GetCurve(i), m_octree);

            m_curvemeshes[i]->Mesh(Nodes, Edges);

        }

        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            if(m_verbose)
                cout << endl << "Surface: " <<  i <<  endl;
            m_surfacemeshes[i] =
                MemoryManager<SurfaceMesh>::AllocateSharedPtr(i,m_verbose,
                    m_cad->GetSurf(i), m_octree,
                    m_curvemeshes,m_order);

            m_surfacemeshes[i]->Mesh(Nodes,Edges,Tris);

        }

        if(m_verbose)
        {
            cout << endl << "Surface mesh statistics" << endl;
            cout << "\tNodes: " << Nodes.size() << endl;
            cout << "\tEdges: " << Edges.size() << endl;
            cout << "\tTriangles " << Tris.size() << endl;
            cout << "\tEuler-Poincaré characteristic: " << Nodes.size()-
                                                           Edges.size()+
                                                           Tris.size() << endl;
        }

        Validate();

        //Optimise();
    }

    void SurfaceMeshing::HOSurf()
    {
        if(m_verbose)
            cout << endl << "High-Order Surface meshing" << endl;

        for(int i = 0; i < Edges.size(); i++)
        {
            MeshEdgeSharedPtr e = Edges[i];
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
                    Array<OneD, NekDouble> loc;
                    c->P(ti[i],loc);
                    MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                            AllocateSharedPtr(Nodes.size(),loc[0],
                                              loc[1],loc[2]);
                    nn->SetCurve(c->GetID(),ti[i]);
                    NekDouble u,v;
                    s1->locuv(u,v,loc);
                    nn->SetSurf(Surfs[0],u,v);
                    s2->locuv(u,v,loc);
                    nn->SetSurf(Surfs[1],u,v);
                    honodes[i-1] = Nodes.size();
                    Nodes[Nodes.size()] = nn;
                }

                e->SetHONodes(honodes);

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
                    loc = s->P(uvb[0]+i*(uve[0]-uvb[0])/m_order,
                               uvb[1]+i*(uve[1]-uvb[1])/m_order);
                    MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                            AllocateSharedPtr(Nodes.size(),loc[0],
                                                loc[1],loc[2]);
                    nn->SetSurf(e->GetSurf(),
                                uvb[0]+i*(uve[0]-uvb[0])/m_order,
                                uvb[1]+i*(uve[1]-uvb[1])/m_order);
                    honodes[i-1] = Nodes.size();
                    Nodes[Nodes.size()] = nn;

                }

                e->SetHONodes(honodes);
            }
        }
    }

    void SurfaceMeshing::Optimise()
    {
        map<int, MeshNodeSharedPtr>::iterator it;
        for(it = Nodes.begin(); it!=Nodes.end(); it++)
        {
            if(it->second->IsOnACurve())
                continue;

            NekDouble d = m_octree->Query(it->second->GetLoc());

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

            //NekDouble u0=0.0,v0=0.0,fu=0.0,dfu=0.0,fv=0.0,dfv=0.0;

        }
    }

    void SurfaceMeshing::Validate()
    {
        if(m_verbose)
            cout << endl << "Verifying surface mesh" << endl;

        if(m_cad->GetEPC() != Nodes.size()-Edges.size()+Tris.size())
        {
            if(m_verbose)
                cout << "\tFailed" << endl;
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
                    cout << "\tFailed" << endl;
                ASSERTL0(false,"edge not listed twice");
            }
        }
        if(m_verbose)
            cout << "\tPassed" << endl;
    }


}
}
