////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.cpp
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
//  Description: tet meshing methods
//
////////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <NekMeshUtils/BLMeshing/BLMesh.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

#include <ANN/ANN.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

NekDouble BLMesh::Visability(vector<ElementSharedPtr> tris, Array<OneD, NekDouble> N)
{
    NekDouble mn = numeric_limits<double>::max();

    for(int i = 0; i < tris.size(); i++)
    {
        vector<NodeSharedPtr> ns = tris[i]->GetVertexList();

        Array<OneD, NekDouble> tmp(3,0.0);
        tmp[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                 (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y);
        tmp[1] = (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x) -
                 (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z);
        tmp[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                 (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x);

        NekDouble mt = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
        mt = sqrt(mt);
        NekDouble dt = tmp[0]*N[0]/mt + tmp[1]*N[1]/mt + tmp[2]*N[2]/mt;
        mn = min(mn,dt);
    }
    return mn;
}

Array<OneD, NekDouble> BLMesh::GetNormal(vector<ElementSharedPtr> tris)
{
    //compile list of normals
    vector<Array<OneD, NekDouble> > N;
    for(int i = 0; i < tris.size(); i++)
    {
        vector<NodeSharedPtr> ns = tris[i]->GetVertexList();

        Array<OneD, NekDouble> tmp(3,0.0);
        tmp[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                 (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y);
        tmp[1] = (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x) -
                 (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z);
        tmp[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                 (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x);

        NekDouble mt = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
        mt = sqrt(mt);

        tmp[0] /= mt;
        tmp[1] /= mt;
        tmp[2] /= mt;

        N.push_back(tmp);
    }

    vector<NekDouble> w(N.size());
    Array<OneD, NekDouble> Np(3,0.0);

    for(int i = 0; i < N.size(); i++)
    {
        w[i] = 1.0/N.size();
    }

    for(int i = 0; i < N.size(); i++)
    {
        Np[0] += w[i] * N[i][0];
        Np[1] += w[i] * N[i][1];
        Np[2] += w[i] * N[i][2];
    }
    NekDouble mag = sqrt(Np[0]*Np[0] + Np[1]*Np[1] + Np[2]*Np[2]);
    Np[0] /= mag;
    Np[1] /= mag;
    Np[2] /= mag;

    Array<OneD, NekDouble> Ninital = Np;

    NekDouble dot = 0.0;
    int ct = 0;
    while(fabs(dot - 1) > 1e-3)
    {
        ct++;
        vector<NekDouble> a(N.size());
        NekDouble aSum = 0.0;
        for(int i = 0; i < N.size(); i++)
        {
            a[i] = acos(Np[0]*N[i][0] + Np[1]*N[i][1] + Np[2]*N[i][2]);

            aSum += a[i];
        }

        NekDouble wSum = 0.0;
        for(int i = 0; i < N.size(); i++)
        {
            w[i] = w[i] * a[i] / aSum;

            wSum += w[i];
        }

        for(int i = 0; i < N.size(); i++)
        {
            w[i] /= wSum;
        }

        Array<OneD, NekDouble> NpN(3,0.0);
        for(int i = 0; i < N.size(); i++)
        {
            NpN[0] += w[i] * N[i][0];
            NpN[1] += w[i] * N[i][1];
            NpN[2] += w[i] * N[i][2];
        }
        mag = sqrt(NpN[0]*NpN[0] + NpN[1]*NpN[1] + NpN[2]*NpN[2]);
        NpN[0] /= mag;
        NpN[1] /= mag;
        NpN[2] /= mag;

        Np[0] = 0.5* NpN[0] + (1.0-0.5)*Np[0];
        Np[1] = 0.5* NpN[1] + (1.0-0.5)*Np[1];
        Np[2] = 0.5* NpN[2] + (1.0-0.5)*Np[2];

        dot = Np[0] * NpN[0] + Np[1] * NpN[1] + Np[2] * NpN[2];

        if(ct > 100000)
        {
            Np = Ninital;
            break;
        }
    }

    Array<OneD, NekDouble> bestN = Np;

    NekDouble val = -1.0*numeric_limits<double>::max();
    NekDouble dtheta = 3.142/5.0;
    NekDouble dphi = 3.142/5.0;
    while(dtheta > 3.142/300.0)
    {
        NekDouble theta0 = acos(bestN[2]);
        NekDouble phi0   = atan2(bestN[1],bestN[0]);

        //sample grid
        for(int i = -10; i <= 10; i++)
        {
            for(int j = -10; j <= 10; j++)
            {
                Array<OneD, NekDouble> tmp(3);
                NekDouble theta = theta0 + i * dtheta;
                NekDouble phi   = phi0 + j * dphi;
                tmp[0] = sin(theta) * cos(phi);
                tmp[1] = sin(theta) * sin(phi);
                tmp[2] = cos(theta);

                NekDouble valt = Visability(tris,tmp);

                if(valt > val)
                {
                    val = valt;
                    bestN = tmp;
                }
            }
        }

        dtheta /= 2.0;
        dphi /= 2.0;
    }

    Np = bestN;

    return Np;
}

void BLMesh::Mesh()
{
    //At this stage the surface mesh is complete and the elements know their
    //neigbours through element links in the edges,,

    // here elements are made for the boundary layer they will need to know
    // links (maybe facelinks), so that the tetmeshing module can extract the
    // surface upon which it needs to mesh (top of the bl and the rest of the
    // surface).

    //need a map from vertex idx to surface elements
    //but do not care about triangles which are not in the bl
    map<int, vector<ElementSharedPtr> > nIdxToTri;
    map<NodeSharedPtr, blInfo> nodeToBL;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        //orientate the triangle
        if(m_cad->GetSurf(m_mesh->m_element[2][i]->CADSurfId)->IsReversedNormal())
        {
            m_mesh->m_element[2][i]->Flip();
        }

        vector<unsigned int>::iterator f = find(m_blsurfs.begin(), m_blsurfs.end(),
                                                m_mesh->m_element[2][i]->CADSurfId);

        if(f == m_blsurfs.end())
        {
            //if this triangle is not in bl surfs continue
            continue;
        }

        vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            nIdxToTri[ns[j]->m_id].push_back(m_mesh->m_element[2][i]);
        }
    }

    int nlayers = 10;
    NekDouble r = 1.05;
    //NekDouble a = (1.0 - r) / (1.0 - pow(r,nlayers+1));
    NekDouble a = 1.0/nlayers;

    NekDouble blprog[nlayers+1];
    blprog[0] = 0.0;
    for(int i = 1; i <= nlayers; i++)
    {
        blprog[i] = blprog[i-1] + m_bl * a;// * pow(r,i);
    }

    set<int> symSurfs;

    NodeSet::iterator it;
    int ct = 0;
    int failed = 0;
    for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++, ct++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                ct, m_mesh->m_vertexSet.size(), "Build info\t");
        }

        vector<pair<int, CADSurfSharedPtr> > ss = (*it)->GetCADSurfs();
        vector<unsigned int> surfs;
        for(int i = 0; i < ss.size(); i++)
        {
            surfs.push_back(ss[i].first);
        }
        sort(surfs.begin(), surfs.end());
        vector<unsigned int> inter, diff;

        set_intersection(m_blsurfs.begin(), m_blsurfs.end(),
                         surfs.begin(), surfs.end(),
                         back_inserter(inter));
        set_symmetric_difference(inter.begin(), inter.end(),
                         surfs.begin(), surfs.end(),
                         back_inserter(diff));

        // is somewhere on a bl surface
        if (inter.size() > 0)
        {
            //initialise a new bl boudnary node
            blInfo bln;
            bln.bl = blprog[1];
            bln.oNode = (*it);
            bln.symsurf = 0;

            map<int, vector<ElementSharedPtr> >::iterator g = nIdxToTri.find((*it)->m_id);
            ASSERTL0(g != nIdxToTri.end(),"failed to find");

            //calculate mesh normal using normal average as first guess
            bln.N = Array<OneD, NekDouble> (3,0.0);
            for(int i = 0; i < g->second.size(); i++)
            {
                vector<NodeSharedPtr> ns = g->second[i]->GetVertexList();

                Array<OneD, NekDouble> tmp(3,0.0);
                tmp[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                         (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y);
                tmp[1] = (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x) -
                         (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z);
                tmp[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                         (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x);

                NekDouble mt = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                mt = sqrt(mt);

                bln.N[0] += tmp[0] / mt;
                bln.N[1] += tmp[1] / mt;
                bln.N[2] += tmp[2] / mt;
            }
            NekDouble mag = 0.0;
            for(int k = 0; k < 3; k++)
            {
                mag += bln.N[k]*bln.N[k];
            }
            mag = sqrt(mag);
            for(int k = 0; k < 3; k++)
            {
                bln.N[k] /= mag;
            }

            if(Visability(g->second,bln.N) < 0.342)
            {
                bln.N = GetNormal(g->second);
            }

            if(Visability(g->second,bln.N) < 0.0)
            {
                cout << "failed " << (*it)->m_x << " " << (*it)->m_y << " "
                                  << (*it)->m_z << " " << Visability(g->second,bln.N) << endl;
                failed++;
            }

            Array<OneD, NekDouble> loc = (*it)->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bln.N[k] * bln.bl;
            }
            bln.pNode = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                                            loc[0], loc[1], loc[2]));

            nodeToBL[bln.pNode] = bln;

            if(diff.size() > 0)
            {
                //if the diff size is greater than 1 there is a curve that needs remeshing
                ASSERTL0(diff.size() <= 1,"not setup for curve bl refinement");
                symSurfs.insert(diff[0]);
                bln.symsurf = diff[0];
            }

            blData[(*it)] = bln;
        }
    }

    if(m_mesh->m_verbose)
    {
        cout << endl;
    }

    ASSERTL0(failed == 0, "some normals failed to generate");

    map<NodeSharedPtr, blInfo>::iterator bit;

    //make prisms

    map<int,int> nm;
    map<ElementSharedPtr,ElementSharedPtr> priToTri;
    map<ElementSharedPtr,ElementSharedPtr> priToPsd;

    ElmtConfig pconf(LibUtilities::ePrism,1,false,false);
    ElmtConfig tconf(LibUtilities::eTriangle,1,false,false);

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int>::iterator f = find(m_blsurfs.begin(), m_blsurfs.end(),
                                                el->CADSurfId);

        if(f == m_blsurfs.end())
        {
            //if this triangle is not in bl surfs continue
            continue;
        }

        vector<NodeSharedPtr> tn(3); //nodes for pseduo surface
        vector<NodeSharedPtr> pn(6); //all prism nodes
        vector<NodeSharedPtr> n = el->GetVertexList();

        nm[0] = 0;
        nm[1] = 3;
        nm[2] = 4;
        nm[3] = 5;
        nm[4] = 1;
        nm[5] = 2;

        for(int j = 0; j < 3; j++)
        {
            pn[nm[j*2+0]] = n[j];
            pn[nm[j*2+1]] = blData[n[j]].pNode;
            tn[j] = blData[n[j]].pNode;
        }


        vector<int> tags;
        tags.push_back(1); //all prisms are comp 1
        ElementSharedPtr E = GetElementFactory().
                    CreateInstance(LibUtilities::ePrism, pconf, pn, tags);
        E->SetId(i);

        m_mesh->m_element[3].push_back(E);

        //tag of this element doesnt matter so can just be 1
        ElementSharedPtr T = GetElementFactory().
                    CreateInstance(LibUtilities::eTriangle, tconf, tn, tags);
        m_psuedoSurface.push_back(T);

        priToTri[E] = el;
        priToPsd[E] = T;
    }

    vector<ElementSharedPtr> prisms = m_mesh->m_element[3];
    set<int> stopped; //ids of nodes to no longer grow
    //loop over all the prisms, grow them one step, if it becomes invalid or
    //intersects mark it for stepping back. Step back all the marked elements
    //and then remove them from consideration (but keep intersection testing
    //them)

    //at this point the pseduo surface should be a connectivity fixed entity
    //therefore can be processed for data structures
    NodeSet pseduoNodes;
    EdgeSet pseduoEdges;
    for(int i = 0; i < m_psuedoSurface.size(); i++)
    {
        vector<NodeSharedPtr> ns = m_psuedoSurface[i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            pseduoNodes.insert(ns[j]);
        }

        vector<EdgeSharedPtr> es = m_psuedoSurface[i]->GetEdgeList();
        for(int j = 0; j < es.size(); j++)
        {
            pair<EdgeSet::iterator,bool> testIns;
            testIns = pseduoEdges.insert(es[j]);
            if (testIns.second)
            {
                (*testIns.first)->m_elLink.push_back(pair<ElementSharedPtr,int>(m_psuedoSurface[i],j));
            }
            else
            {
                m_psuedoSurface[i]->SetEdge(j, (*testIns.first));
                (*testIns.first)->m_elLink.push_back(pair<ElementSharedPtr,int>(m_psuedoSurface[i],j));
            }
        }
        m_psuedoSurface[i]->SetId(i);
    }

    //need to build a list of nodes to neigbours
    map<ElementSharedPtr, ElementSharedPtr>::iterator eit;
    map<NodeSharedPtr, NodeSet> nodeToNear;
    for(eit = priToTri.begin(); eit != priToTri.end(); eit++)
    {
        vector<EdgeSharedPtr> es = eit->second->GetEdgeList();
        for(int j = 0; j < es.size(); j++)
        {
            nodeToNear[es[j]->m_n1].insert(es[j]->m_n2);
            nodeToNear[es[j]->m_n2].insert(es[j]->m_n1);
        }
    }

    map<int, EdgeSharedPtr> pedges;
    int ect = 0;
    EdgeSet::iterator et;
    for(et = pseduoEdges.begin(); et != pseduoEdges.end(); et++, ect++)
    {
        pedges[ect] = (*et);
    }

    for(int i = 2; i <= nlayers; i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, nlayers, "layers\t");
        }

        cout << endl;
        cout << "stopped " << stopped.size() << endl;
        vector<ElementSharedPtr> revert;
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            set<int>::iterator f = stopped.find(bit->first->m_id);
            if(f != stopped.end())
            {
                continue;
            }

            bit->second.bl = blprog[i];

            Array<OneD, NekDouble> loc = bit->first->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bit->second.N[k] * bit->second.bl;
            }

            bit->second.pNode->m_x = loc[0];
            bit->second.pNode->m_y = loc[1];
            bit->second.pNode->m_z = loc[2];
        }

        set<int> toAdd;
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            set<int>::iterator f = stopped.find(bit->first->m_id);
            if(f != stopped.end())
            {
                continue;
            }

            //if any of the neigbours are previously stopped we need
            //to stop this one for smoothness
            map<NodeSharedPtr, NodeSet>::iterator l = nodeToNear.find(bit->first);
            ASSERTL0(l!=nodeToNear.end(),"not found");
            NodeSet::iterator nit;
            for(nit = l->second.begin(); nit != l->second.end(); nit++)
            {
                f = stopped.find((*nit)->m_id);
                if(f != stopped.end())
                {
                    toAdd.insert(bit->first->m_id);
                }
            }
        }
        set<int>::iterator f;
        for(f = toAdd.begin(); f != toAdd.end(); f++)
        {
            stopped.insert((*f));
        }

        vector<ElementSharedPtr> tmp = prisms;
        prisms.clear();
        for(int j = 0; j < tmp.size(); j++)
        {
            map<ElementSharedPtr, ElementSharedPtr>::iterator f = priToTri.find(tmp[j]);
            ASSERTL0(f != priToTri.end(), "not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            set<int>::iterator f1 = stopped.find(ns[0]->m_id);
            set<int>::iterator f2 = stopped.find(ns[1]->m_id);
            set<int>::iterator f3 = stopped.find(ns[2]->m_id);

            if(!(f1 != stopped.end() && f2 != stopped.end() && f3 != stopped.end()))
            {
                prisms.push_back(tmp[j]);
            }
        }


        for(int j = 0; j < prisms.size(); j++)
        {
            SpatialDomains::GeometrySharedPtr geom =
                                        prisms[j]->GetGeom(m_mesh->m_spaceDim);
            SpatialDomains::GeomFactorsSharedPtr gfac =
                                                geom->GetGeomFactors();

            if(!gfac->IsValid())
            {
                revert.push_back(prisms[j]);
            }
        }

        //this is where proximity goes
        //build a ANN tree from the center point of triangles
        ANNpointArray dataPts;
        ANNpoint queryPt;
        ANNidxArray nnIdx;
        ANNdistArray dists;
        ANNkd_tree* kdTree;
        queryPt = annAllocPt(3);
        dataPts = annAllocPts(pedges.size(), 3);

        NekDouble maxsize = 0.0;

        map<int, EdgeSharedPtr>::iterator etr;
        for(etr = pedges.begin(); etr != pedges.end(); etr++)
        {
            dataPts[etr->first][0] = (etr->second->m_n1->m_x + etr->second->m_n2->m_x) / 2.0;
            dataPts[etr->first][1] = (etr->second->m_n1->m_y + etr->second->m_n2->m_y) / 2.0;
            dataPts[etr->first][2] = (etr->second->m_n1->m_z + etr->second->m_n2->m_z) / 2.0;

            NekDouble size = sqrt((etr->second->m_n1->m_x - etr->second->m_n2->m_x) *
                                  (etr->second->m_n1->m_x - etr->second->m_n2->m_x) +
                                  (etr->second->m_n1->m_y - etr->second->m_n2->m_y) *
                                  (etr->second->m_n1->m_y - etr->second->m_n2->m_y) +
                                  (etr->second->m_n1->m_z - etr->second->m_n2->m_z) *
                                  (etr->second->m_n1->m_z - etr->second->m_n2->m_z));
            maxsize = max(size,maxsize);
        }

        kdTree = new ANNkd_tree(dataPts, pedges.size(), 3);

        cout << "prisms " << prisms.size() << endl;

        for(int j = 0; j < prisms.size(); j++)
        {
            map<ElementSharedPtr,ElementSharedPtr>::iterator f = priToPsd.find(prisms[j]);
            ASSERTL0(f != priToPsd.end(),"not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            queryPt[0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x) / 3.0;
            queryPt[1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y) / 3.0;
            queryPt[2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z) / 3.0;

            int sample = 0;
            sample = kdTree->annkFRSearch(queryPt, maxsize*maxsize*4, sample);
            nnIdx = new ANNidx[sample];
            dists = new ANNdist[sample];
            kdTree->annkSearch(queryPt, sample, nnIdx, dists);

            int tested = 0;
            int hited = 0;
            for(int s = 0; s < sample; s++)
            {
                EdgeSharedPtr e = pedges[nnIdx[s]];

                if(ns[0] == e->m_n1 ||
                   ns[0] == e->m_n2 ||
                   ns[1] == e->m_n1 ||
                   ns[1] == e->m_n2 ||
                   ns[2] == e->m_n1 ||
                   ns[2] == e->m_n2)
                {
                    continue;
                }

                DNekMat A(3,3,0.0);
                NekVector<NekDouble> B(3,0.0);
                A(0,0) = (e->m_n2->m_x - e->m_n1->m_x) * -1.0;
                A(1,0) = (e->m_n2->m_y - e->m_n1->m_y) * -1.0;
                A(2,0) = (e->m_n2->m_z - e->m_n1->m_z) * -1.0;
                A(0,1) = ns[1]->m_x - ns[0]->m_x;
                A(1,1) = ns[1]->m_y - ns[0]->m_y;
                A(2,1) = ns[1]->m_z - ns[0]->m_z;
                A(0,2) = ns[2]->m_x - ns[0]->m_x;
                A(1,2) = ns[2]->m_y - ns[0]->m_y;
                A(2,2) = ns[2]->m_z - ns[0]->m_z;

                NekDouble det = A(0,0) * (A(1,1)*A(2,2) - A(2,1)*A(1,2))
                               -A(0,1) * (A(1,0)*A(2,2) - A(2,0)*A(1,2))
                               +A(0,2) * (A(1,0)*A(2,1) - A(2,0)*A(1,1));

                if(fabs(det) < 1e-15)
                {
                    //no intersecton
                    continue;
                }
                B(0) = e->m_n1->m_x - ns[0]->m_x;
                B(1) = e->m_n1->m_y - ns[0]->m_y;
                B(2) = e->m_n1->m_z - ns[0]->m_z;

                tested++;

                A.Invert();

                NekVector<NekDouble> X = A * B;

                //check triangle intersecton
                if(X(1) > -1e-6 && X(2) > 1e-6 && X(1) + X(2) < 1.000001
                   && X(0) > -1e-6 && X(0) < 1.000001)
                {
                    //hit
                    hited++;
                    revert.push_back(prisms[j]);
                    break;
                }
            }
            //cout << sample << " " << tested << " " << hited << endl;
        }

        tmp = prisms;
        prisms.clear();
        set<int> reverted;
        for(int j = 0; j < revert.size(); j++)
        {
            set<int>::iterator rev = reverted.find(revert[j]->GetId());
            if(rev != reverted.end())
            {
                //prism has already been reverted
                continue;
            }
            reverted.insert(revert[j]->GetId());
            map<ElementSharedPtr, ElementSharedPtr>::iterator f = priToTri.find(revert[j]);
            ASSERTL0(f != priToTri.end(), "not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            int skipped = 0;
            for(int k = 0; k < ns.size(); k++)
            {
                set<int>::iterator s = stopped.find(ns[k]->m_id);
                if(s != stopped.end())
                {
                    skipped++;
                    continue;
                }

                map<NodeSharedPtr, blInfo>::iterator bli = blData.find(ns[k]);
                bli->second.bl = blprog[i-1];

                Array<OneD, NekDouble> loc = bli->first->GetLoc();
                for(int k = 0; k < 3; k++)
                {
                    loc[k] += bli->second.N[k] * bli->second.bl;
                }

                bli->second.pNode->m_x = loc[0];
                bli->second.pNode->m_y = loc[1];
                bli->second.pNode->m_z = loc[2];

                stopped.insert(ns[k]->m_id);
            }
            ASSERTL0(skipped < 3, "had reverted and skipped all three");
        }

        for(int j = 0; j < tmp.size(); j++)
        {
            set<int>::iterator f = reverted.find(tmp[j]->GetId());
            if(f == reverted.end())
            {
                prisms.push_back(tmp[j]);
            }
        }

    }

    /*
    //collision detection and element shrinking

    vector<ElementSharedPtr> intersecting;
    bool repeat = true;
    while(repeat)
    {
        repeat = false;
        NodeSet done;

        //build a ANN tree from the center point of triangles
        ANNpointArray dataPts;
        ANNpoint queryPt;
        ANNidxArray nnIdx;
        ANNdistArray dists;
        ANNkd_tree* kdTree;
        queryPt = annAllocPt(3);
        dataPts = annAllocPts(m_psuedoSurface.size(), 3);

        for(int i = 0; i < m_psuedoSurface.size(); i++)
        {
            vector<NodeSharedPtr> ns = m_psuedoSurface[i]->GetVertexList();
            dataPts[i][0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x ) / 3.0;
            dataPts[i][1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y ) / 3.0;
            dataPts[i][2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z ) / 3.0;
        }

        kdTree = new ANNkd_tree(dataPts, m_psuedoSurface.size(), 3);

        for(int i = 0; i < m_psuedoSurface.size(); i++)
        {
            if (m_mesh->m_verbose)
            {
                LibUtilities::PrintProgressbar(
                    i, m_psuedoSurface.size(), "Proximity sweeping\t");
            }

            vector<NodeSharedPtr> ns = m_psuedoSurface[i]->GetVertexList();

            queryPt[0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x ) / 3.0;
            queryPt[1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y ) / 3.0;
            queryPt[2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z ) / 3.0;
            int sample = 0;
            sample = kdTree->annkFRSearch(queryPt, 0.1*0.1, sample);
            nnIdx = new ANNidx[sample];
            dists = new ANNdist[sample];
            kdTree->annkFRSearch(queryPt, 0.1*0.1, sample, nnIdx, dists);

            EdgeSet toTest;

            for(int j = 0; j < sample; j++)
            {
                ElementSharedPtr el = m_psuedoSurface[nnIdx[j]];
                vector<EdgeSharedPtr> es = el->GetEdgeList();
                for(int k = 0; k < es.size(); k++)
                {
                    if(ns[0] == es[k]->m_n1 ||
                       ns[0] == es[k]->m_n2 ||
                       ns[1] == es[k]->m_n1 ||
                       ns[1] == es[k]->m_n2 ||
                       ns[2] == es[k]->m_n1 ||
                       ns[2] == es[k]->m_n2)
                    {
                        continue;
                    }
                    toTest.insert(es[k]);
                }
            }

            //cout << sample << " " << toTest.size() << endl;

            EdgeSet::iterator eit;
            for(eit = toTest.begin(); eit != toTest.end(); eit++)
            {
                NekDouble A[3][3];
                NekDouble B[3];
                A[0][0] = ((*eit)->m_n2->m_x - (*eit)->m_n1->m_x) * -1.0;
                A[1][0] = ((*eit)->m_n2->m_y - (*eit)->m_n1->m_y) * -1.0;
                A[2][0] = ((*eit)->m_n2->m_z - (*eit)->m_n1->m_z) * -1.0;
                A[0][1] = ns[1]->m_x - ns[0]->m_x;
                A[1][1] = ns[1]->m_y - ns[0]->m_y;
                A[2][1] = ns[1]->m_z - ns[0]->m_z;
                A[0][2] = ns[2]->m_x - ns[0]->m_x;
                A[1][2] = ns[2]->m_y - ns[0]->m_y;
                A[2][2] = ns[2]->m_z - ns[0]->m_z;

                NekDouble det = A[0][0] * (A[1][1]*A[2][2] - A[2][1]*A[1][2])
                               -A[0][1] * (A[1][0]*A[2][2] - A[2][0]*A[1][2])
                               +A[0][2] * (A[1][0]*A[2][1] - A[2][0]*A[1][1]);
                if(fabs(det) < 1e-12)
                {
                    //no intersecton
                    continue;
                }
                B[0] = (*eit)->m_n1->m_x - ns[0]->m_x;
                B[1] = (*eit)->m_n1->m_y - ns[0]->m_y;
                B[2] = (*eit)->m_n1->m_z - ns[0]->m_z;

                NekDouble X[3];

                X[0] = 1.0 / det * (B[0] * (A[1][1]*A[2][2] - A[2][1]*A[1][2]) +
                                    B[1] * (A[2][1]*A[0][2] - A[0][1]*A[2][2]) +
                                    B[2] * (A[0][1]*A[1][2] - A[1][1]*A[0][2]));
                X[1] = 1.0 / det * (B[0] * (A[2][0]*A[1][2] - A[1][0]*A[2][2]) +
                                    B[1] * (A[0][0]*A[2][2] - A[2][0]*A[0][2]) +
                                    B[2] * (A[1][0]*A[0][2] - A[0][0]*A[1][2]));
                X[2] = 1.0 / det * (B[0] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]) +
                                    B[1] * (A[2][0]*A[0][1] - A[0][0]*A[2][1]) +
                                    B[2] * (A[0][0]*A[1][1] - A[1][0]*A[0][1]));

                //check triangle intersecton
                if(X[1] >= 0.0 && X[2] >= 0.0 && X[1] + X[2] <= 1.0
                   && X[0] >= 0.0 && X[0] <= 1.0)
                {
                    //hit
                    intersecting.push_back(m_psuedoSurface[i]);
                    break;
                }
            }
        }

        cout << endl << intersecting.size() << endl << endl;

        if(intersecting.size() > 0)
            repeat = true;

        for(int i = 0; i < intersecting.size(); i++)
        {
            vector<NodeSharedPtr> ns = intersecting[i]->GetVertexList();
            for(int j = 0; j < ns.size(); j++)
            {
                NodeSet::iterator f = done.find(ns[j]);
                if(f == done.end())
                {
                    map<NodeSharedPtr, blInfo>::iterator b = nodeToBL.find(ns[j]);
                    ASSERTL0(b != nodeToBL.end(),"could not find");

                    b->second.bl *= 0.75;

                    Array<OneD, NekDouble> loc = b->second.oNode->GetLoc();
                    for(int k = 0; k < 3; k++)
                    {
                        loc[k] += b->second.N[k] * b->second.bl;
                    }

                    ns[j]->m_x = loc[0];
                    ns[j]->m_y = loc[1];
                    ns[j]->m_z = loc[2];

                    done.insert(ns[j]);
                }
            }
        }

        intersecting.clear();
    }

    //loop over all prisms, if invalid shrink until it is
    //being careful to act on nodes which have already been shrunk
*//*
    //smoothing
    //need to build a list of nodes to neigbours
    map<ElementSharedPtr, ElementSharedPtr>::iterator eit;
    map<NodeSharedPtr, NodeSet> nodeToNear;
    for(eit = priToTri.begin(); eit != priToTri.end(); eit++)
    {
        vector<EdgeSharedPtr> es = eit->second->GetEdgeList();
        for(int j = 0; j < es.size(); j++)
        {
            nodeToNear[es[j]->m_n1].insert(es[j]->m_n2);
            nodeToNear[es[j]->m_n2].insert(es[j]->m_n1);
        }
    }

    bool repeat = true;
    while(repeat)
    {
        repeat = false;
        bool repeat2 = true;
        while(repeat2)
        {
            repeat2 = false;
            for(bit = blData.begin(); bit != blData.end(); bit++)
            {
                map<NodeSharedPtr, NodeSet>::iterator mit = nodeToNear.find(bit->first);
                ASSERTL0(mit != nodeToNear.end(),"not found");
                NodeSet::iterator nit;
                for(nit = mit->second.begin(); nit != mit->second.end(); nit++)
                {
                    map<NodeSharedPtr, blInfo>::iterator bli = blData.find((*nit));
                    ASSERTL0(bli != blData.end(),"not found");
                    if(bli->second.bl < 0.75 * bit->second.bl)
                    {
                        bit->second.bl = bli->second.bl * 1.3;
                        bit->second.pNode->m_x = bit->first->m_x + bit->second.N[0] * bit->second.bl;
                        bit->second.pNode->m_y = bit->first->m_y + bit->second.N[1] * bit->second.bl;
                        bit->second.pNode->m_z = bit->first->m_z + bit->second.N[2] * bit->second.bl;
                        repeat2 = true;
                    }
                }
            }
        }

        for(int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            if (m_mesh->m_verbose)
            {
                LibUtilities::PrintProgressbar(
                    i, m_mesh->m_element[3].size(), "Invalidity sweeping\t");
            }

            ElementSharedPtr el = m_mesh->m_element[3][i];
            SpatialDomains::GeometrySharedPtr geom =
                                                el->GetGeom(m_mesh->m_spaceDim);
            SpatialDomains::GeomFactorsSharedPtr gfac =
                                                geom->GetGeomFactors();

            map<ElementSharedPtr, ElementSharedPtr>::iterator j = priToTri.find(el);
            ASSERTL0(j != priToTri.end(), "not found");
            vector<NodeSharedPtr> ns = j->second->GetVertexList();

            while(!gfac->IsValid())
            {
                NekDouble maxbl = max(blData[ns[0]].bl, blData[ns[1]].bl);
                maxbl = max(maxbl, blData[ns[2]].bl);

                if(maxbl < 1E-6)
                {
                    cout << "shrunk element too far, invalid mesh" << endl;
                    abort();
                    break;
                }

                int ct = 0;
                for(int j = 0; j < 3; j++)
                {
                    map<NodeSharedPtr, blInfo>::iterator bli = blData.find(ns[j]);
                    ASSERTL0(bli != blData.end(), "not found");
                    if(bli->second.bl <= maxbl * 0.75)
                    {
                        ct++;
                        continue;
                    }
                    ASSERTL0(ct < 3,"skipped all 3");

                    bli->second.bl *= 0.75;

                    bli->second.pNode->m_x = ns[j]->m_x + bli->second.N[0] * bli->second.bl;
                    bli->second.pNode->m_y = ns[j]->m_y + bli->second.N[1] * bli->second.bl;
                    bli->second.pNode->m_z = ns[j]->m_z + bli->second.N[2] * bli->second.bl;
                    repeat = true;

                }

                geom = el->GetGeom(m_mesh->m_spaceDim);
                gfac = geom->GetGeomFactors();
            }
        }
    }

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

    m_symSurfs = vector<int>(symSurfs.begin(), symSurfs.end());
    //compile a map of the nodes needed for systemetry surfs
    for(int i = 0; i < m_symSurfs.size(); i++)
    {
        int s = m_symSurfs[i];
        map<NodeSharedPtr, NodeSharedPtr> nmap;
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(bit->second.symsurf == s)
            {
                nmap[bit->first] = bit->second.pNode;
            }
        }
        m_symNodes[s] = nmap;
    }
*/
}
}
}
