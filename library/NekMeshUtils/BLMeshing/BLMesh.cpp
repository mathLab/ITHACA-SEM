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
    NekDouble dtheta = M_PI/5.0;
    NekDouble dphi = M_PI/5.0;
    while(dtheta > M_PI/300.0)
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
            bln.bl = 1;
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
                cerr << "failed " << (*it)->m_x << " " << (*it)->m_y << " "
                                  << (*it)->m_z << " "
                                  << Visability(g->second,bln.N) << endl;
                failed++;
            }

            Array<OneD, NekDouble> loc = (*it)->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bln.N[k] * blprog[bln.bl];
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
                bln.onSym = true;
            }
            else
            {
                bln.onSym = false;
            }

            blData[(*it)] = bln;
        }
    }

    m_symSurfs = vector<int>(symSurfs.begin(), symSurfs.end());

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
    map<NodeSharedPtr, vector<ElementSharedPtr> > nodeToNearPri;
    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[3][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            nodeToNearPri[ns[j]].push_back(el);
        }
    }

    map<int, EdgeSharedPtr> pedges;
    int ect = 0;
    EdgeSet::iterator et;
    for(et = pseduoEdges.begin(); et != pseduoEdges.end(); et++, ect++)
    {
        pedges[ect] = (*et);
    }

    {
        //before iterating over the layers, do an intial intersection test
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

        for(int j = 0; j < m_mesh->m_element[3].size(); j++)
        {
            map<ElementSharedPtr,ElementSharedPtr>::iterator f = priToPsd.find(m_mesh->m_element[3][j]);
            ASSERTL0(f != priToPsd.end(),"not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            ElementSharedPtr el = m_mesh->m_element[3][j];
            SpatialDomains::GeometrySharedPtr geom =
                                                el->GetGeom(m_mesh->m_spaceDim);
            SpatialDomains::GeomFactorsSharedPtr gfac =
                                                geom->GetGeomFactors();

            if(!gfac->IsValid())
            {
                cout << "intial invalid element" << endl;
                cout << ns[0]->m_x << " " << ns[0]->m_y << " " << ns[0]->m_z << endl;
            }

            queryPt[0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x) / 3.0;
            queryPt[1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y) / 3.0;
            queryPt[2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z) / 3.0;

            int sample = 0;
            sample = kdTree->annkFRSearch(queryPt, maxsize*maxsize*2, sample);
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

                NekDouble A0,A1,A2,A3,A4,A5,A6,A7,A8;
                NekDouble B0,B1,B2;
                A0 = (e->m_n2->m_x - e->m_n1->m_x) * -1.0;
                A1 = (e->m_n2->m_y - e->m_n1->m_y) * -1.0;
                A2 = (e->m_n2->m_z - e->m_n1->m_z) * -1.0;
                A3 = ns[1]->m_x - ns[0]->m_x;
                A4 = ns[1]->m_y - ns[0]->m_y;
                A5 = ns[1]->m_z - ns[0]->m_z;
                A6 = ns[2]->m_x - ns[0]->m_x;
                A7 = ns[2]->m_y - ns[0]->m_y;
                A8 = ns[2]->m_z - ns[0]->m_z;

                NekDouble det = A0 * (A4*A8 - A7*A5)
                               -A3 * (A1*A8 - A7*A2)
                               +A6 * (A1*A5 - A4*A2);

                if(fabs(det) < 1e-15)
                {
                    //no intersecton
                    continue;
                }
                B0 = e->m_n1->m_x - ns[0]->m_x;
                B1 = e->m_n1->m_y - ns[0]->m_y;
                B2 = e->m_n1->m_z - ns[0]->m_z;

                tested++;

                NekDouble X0,X1,X2;

                X0 = B0 * (A4*A8 - A7*A5)
                    -A3 * (B1*A8 - A7*B2)
                    +A6 * (B1*A5 - A4*B2);

                X1 = A0 * (B1*A8 - A7*B2)
                    -B0 * (A1*A8 - A7*A2)
                    +A6 * (A1*B2 - B1*A2);

                X2 = A0 * (A4*B2 - B1*A5)
                    -A3 * (A1*B2 - B1*A2)
                    +B0 * (A1*A5 - A4*A2);


                X0 /= det;
                X1 /= det;
                X2 /= det;

                //check triangle intersecton
                if(X1 >= 0.0 && X2 >=0.0 && X1 + X2 <= 1.0
                   && X0 >= 0.0 && X0 <= 1.0)
                {
                    cout << "initialy intersecting element" << endl;
                    cout << ns[0]->m_x << " " << ns[0]->m_y << " " << ns[0]->m_z << endl;
                }
            }
        }
    }

    for(int i = 2; i <= nlayers; i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, nlayers, "layers\t");
        }

        cout << endl;
        vector<ElementSharedPtr> revert;
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            set<int>::iterator f = stopped.find(bit->first->m_id);
            if(f != stopped.end())
            {
                continue;
            }

            bit->second.bl = i;

            Array<OneD, NekDouble> loc = bit->first->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bit->second.N[k] * blprog[bit->second.bl];
            }

            bit->second.pNode->m_x = loc[0];
            bit->second.pNode->m_y = loc[1];
            bit->second.pNode->m_z = loc[2];

            if(bit->second.onSym)
            {
                CADSurfSharedPtr s = m_cad->GetSurf(bit->second.symsurf);
                Array<OneD, NekDouble> uv(2);
                s->ProjectTo(loc,uv);
                bit->second.pNode->m_x = loc[0];
                bit->second.pNode->m_y = loc[1];
                bit->second.pNode->m_z = loc[2];
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

        for(int j = 0; j < prisms.size(); j++)
        {
            ElementSharedPtr el = prisms[j];
            SpatialDomains::GeometrySharedPtr geom =
                                                el->GetGeom(m_mesh->m_spaceDim);
            SpatialDomains::GeomFactorsSharedPtr gfac =
                                                geom->GetGeomFactors();

            if(!gfac->IsValid())
            {
                revert.push_back(prisms[j]);
                continue;
            }

            map<ElementSharedPtr,ElementSharedPtr>::iterator f = priToPsd.find(prisms[j]);
            ASSERTL0(f != priToPsd.end(),"not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            queryPt[0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x) / 3.0;
            queryPt[1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y) / 3.0;
            queryPt[2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z) / 3.0;

            int sample = 0;
            sample = kdTree->annkFRSearch(queryPt, maxsize*maxsize*2, sample);
            nnIdx = new ANNidx[sample];
            dists = new ANNdist[sample];
            kdTree->annkSearch(queryPt, sample, nnIdx, dists);

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

                NekDouble A0,A1,A2,A3,A4,A5,A6,A7,A8;
                NekDouble B0,B1,B2;
                A0 = (e->m_n2->m_x - e->m_n1->m_x) * -1.0;
                A1 = (e->m_n2->m_y - e->m_n1->m_y) * -1.0;
                A2 = (e->m_n2->m_z - e->m_n1->m_z) * -1.0;
                A3 = ns[1]->m_x - ns[0]->m_x;
                A4 = ns[1]->m_y - ns[0]->m_y;
                A5 = ns[1]->m_z - ns[0]->m_z;
                A6 = ns[2]->m_x - ns[0]->m_x;
                A7 = ns[2]->m_y - ns[0]->m_y;
                A8 = ns[2]->m_z - ns[0]->m_z;

                NekDouble det = A0 * (A4*A8 - A7*A5)
                               -A3 * (A1*A8 - A7*A2)
                               +A6 * (A1*A5 - A4*A2);

                if(fabs(det) < 1e-15)
                {
                    //no intersecton
                    continue;
                }
                B0 = e->m_n1->m_x - ns[0]->m_x;
                B1 = e->m_n1->m_y - ns[0]->m_y;
                B2 = e->m_n1->m_z - ns[0]->m_z;

                NekDouble X0,X1,X2;

                X0 = B0 * (A4*A8 - A7*A5)
                    -A3 * (B1*A8 - A7*B2)
                    +A6 * (B1*A5 - A4*B2);

                X1 = A0 * (B1*A8 - A7*B2)
                    -B0 * (A1*A8 - A7*A2)
                    +A6 * (A1*B2 - B1*A2);

                X2 = A0 * (A4*B2 - B1*A5)
                    -A3 * (A1*B2 - B1*A2)
                    +B0 * (A1*A5 - A4*A2);


                X0 /= det;
                X1 /= det;
                X2 /= det;

                //check triangle intersecton
                if(X1 > -1e-6 && X2 > 1e-6 && X1 + X2 < 1.000001
                   && X0 > -1e-6 && X0 < 1.000001)
                {
                    //hit
                    revert.push_back(prisms[j]);
                    break;
                }
            }
        }

        delete [] nnIdx;
        delete [] dists;
        delete kdTree;

        //at this point we have a list of elements which were made invalid
        //by the advancement of the layer
        //now we need to loop over these and revert their heights
        //in doing so connecting elements may also be made invalid therefore
        //this is a recursive process

        cout << "initial reversion " << revert.size() << endl;

        while(revert.size() > 0)
        {
            set<int> reverted;
            vector<ElementSharedPtr> toCheck;
            for(int j = 0; j < revert.size(); j++)
            {
                map<ElementSharedPtr, ElementSharedPtr>::iterator f = priToTri.find(revert[j]);
                ASSERTL0(f != priToTri.end(), "not found");
                vector<NodeSharedPtr> ns = f->second->GetVertexList();

                for(int k = 0; k < ns.size(); k++)
                {
                    map<NodeSharedPtr, blInfo>::iterator bli = blData.find(ns[k]);
                    if(bli->second.bl < i-1)
                    {
                        cout << "not smooth specification " << i << " " << bli->second.bl << endl;
                    }
                    if(bli->second.bl == i-1)
                    {
                        continue;
                    }

                    bli->second.bl = i-1;

                    Array<OneD, NekDouble> loc = bli->first->GetLoc();
                    for(int k = 0; k < 3; k++)
                    {
                        loc[k] += bli->second.N[k] * blprog[bli->second.bl];
                    }

                    bli->second.pNode->m_x = loc[0];
                    bli->second.pNode->m_y = loc[1];
                    bli->second.pNode->m_z = loc[2];

                    if(bit->second.onSym)
                    {
                        CADSurfSharedPtr s = m_cad->GetSurf(bit->second.symsurf);
                        Array<OneD, NekDouble> uv(2);
                        s->ProjectTo(loc,uv);
                        bit->second.pNode->m_x = loc[0];
                        bit->second.pNode->m_y = loc[1];
                        bit->second.pNode->m_z = loc[2];
                    }

                    stopped.insert(ns[k]->m_id);

                    map<NodeSharedPtr, vector<ElementSharedPtr> >::iterator pit =
                                            nodeToNearPri.find(ns[k]);
                    for(int l = 0; l < pit->second.size(); l++)
                    {
                        if(pit->second[l]->GetId() == revert[j]->GetId())
                        {
                            continue;
                        }
                        toCheck.push_back(pit->second[l]);
                    }
                }
                reverted.insert(revert[j]->GetId());
            }

            revert.clear();

            set<int> ids;
            vector<ElementSharedPtr> tmp = toCheck;
            toCheck.clear();
            for(int k = 0; k < tmp.size(); k++)
            {
                set<int>::iterator f2 = ids.find(tmp[k]->GetId());
                set<int>::iterator f3 = reverted.find(tmp[k]->GetId());
                if(f2 == ids.end() && f3 == reverted.end())
                {
                    toCheck.push_back(tmp[k]);
                    ids.insert(tmp[k]->GetId());
                }
            }

            for(int k = 0; k < toCheck.size(); k++)
            {
                SpatialDomains::GeometrySharedPtr geom =
                                                    toCheck[k]->GetGeom(m_mesh->m_spaceDim);
                SpatialDomains::GeomFactorsSharedPtr gfac =
                                                    geom->GetGeomFactors();

                if(!gfac->IsValid())
                {
                    revert.push_back(toCheck[k]);
                }
            }
        }

        vector<ElementSharedPtr> tmp = prisms;
        prisms.clear();

        set<int> toAdd;

        for(int j = 0; j < tmp.size(); j++)
        {
            map<ElementSharedPtr, ElementSharedPtr>::iterator f = priToTri.find(tmp[j]);
            ASSERTL0(f != priToTri.end(), "not found");
            vector<NodeSharedPtr> ns = f->second->GetVertexList();

            int stop = 0;
            for(int k = 0; k < ns.size(); k++)
            {
                set<int>::iterator s = stopped.find(ns[k]->m_id);
                if(s != stopped.end())
                {
                    stop++;
                }
            }

            if(stop > 0 && stop < 3)
            {
                int mn = 1000, mx = 0;
                for(int k = 0; k < ns.size(); k++)
                {
                    map<NodeSharedPtr, blInfo>::iterator bli = blData.find(ns[k]);
                    mn = min(mn, bli->second.bl);
                    mx = min(mx, bli->second.bl);
                }
                if(mx - mn > 1)
                {
                    cout << "error in smoothness" << endl;
                }
                for(int k = 0; k < ns.size(); k++)
                {
                    toAdd.insert(ns[k]->m_id);
                }
            }
            else if(stop == 0)
            {
                prisms.push_back(tmp[j]);
            }
        }

        set<int>::iterator iit;
        for(iit = toAdd.begin(); iit != toAdd.end(); iit++)
        {
            stopped.insert((*iit));
        }
    }

    if(m_mesh->m_verbose)
    {
        cout << endl;
    }

    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        SpatialDomains::GeometrySharedPtr geom =
                                            m_mesh->m_element[3][i]->GetGeom(m_mesh->m_spaceDim);
        SpatialDomains::GeomFactorsSharedPtr gfac =
                                            geom->GetGeomFactors();

        if(!gfac->IsValid())
        {
            cout << "still got an invalid element" << endl;
        }
    }


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

    for(int j = 0; j < m_mesh->m_element[3].size(); j++)
    {
        map<ElementSharedPtr,ElementSharedPtr>::iterator f = priToPsd.find(m_mesh->m_element[3][j]);
        ASSERTL0(f != priToPsd.end(),"not found");
        vector<NodeSharedPtr> ns = f->second->GetVertexList();

        queryPt[0] = (ns[0]->m_x + ns[1]->m_x + ns[2]->m_x) / 3.0;
        queryPt[1] = (ns[0]->m_y + ns[1]->m_y + ns[2]->m_y) / 3.0;
        queryPt[2] = (ns[0]->m_z + ns[1]->m_z + ns[2]->m_z) / 3.0;

        int sample = 0;
        sample = kdTree->annkFRSearch(queryPt, maxsize*maxsize*2, sample);
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
                cout << "still got an intersection" << endl;
            }
        }
        //cout << sample << " " << tested << " " << hited << endl;
    }

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

    //m_mesh->m_element[2] = m_psuedoSurface;

}
}
}
