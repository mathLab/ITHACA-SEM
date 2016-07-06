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
        vector<unsigned int>::iterator f = find(m_blsurfs.begin(), m_blsurfs.end(),
                                                tris[i]->CADSurfId);
        if(f == m_blsurfs.end())
        {
            //if this triangle is not in inter continue
            continue;
        }

        vector<NodeSharedPtr> ns = tris[i]->GetVertexList();
        if(m_cad->GetSurf(tris[i]->CADSurfId)->IsReversedNormal())
        {
            swap(ns[0],ns[1]);
        }

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
        vector<unsigned int>::iterator f = find(m_blsurfs.begin(), m_blsurfs.end(),
                                                tris[i]->CADSurfId);
        if(f == m_blsurfs.end())
        {
            //if this triangle is not in inter continue
            continue;
        }

        vector<NodeSharedPtr> ns = tris[i]->GetVertexList();
        if(m_cad->GetSurf(tris[i]->CADSurfId)->IsReversedNormal())
        {
            swap(ns[0],ns[1]);
        }

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

    NekDouble dot = 0.0;
    int ct = 0;
    while(fabs(dot - 1) > 1e-4)
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

        ASSERTL0(ct < 10000,"failed to find normal");
    }

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
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            nIdxToTri[ns[j]->m_id].push_back(m_mesh->m_element[2][i]);
        }
    }

    set<int> symSurfs;

    NodeSet::iterator it;
    int ct = 0;
    int failed = 0;
    for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++, ct++)
    {
        //if (m_mesh->m_verbose)
        //{
        //    LibUtilities::PrintProgressbar(
        //        ct, m_mesh->m_vertexSet.size(), "Building BL info");
        //}

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
            bln.bl = m_bl;
            bln.symsurf = 0;

            //calculate mesh normal
            bln.N = Array<OneD, NekDouble> (3,0.0);
            map<int, vector<ElementSharedPtr> >::iterator g = nIdxToTri.find((*it)->m_id);
            for(int i = 0; i < g->second.size(); i++)
            {
                vector<unsigned int>::iterator f = find(m_blsurfs.begin(), m_blsurfs.end(),
                                                        g->second[i]->CADSurfId);
                if(f == m_blsurfs.end())
                {
                    //if this triangle is not in inter continue
                    continue;
                }

                vector<NodeSharedPtr> ns = g->second[i]->GetVertexList();
                if(m_cad->GetSurf(g->second[i]->CADSurfId)->IsReversedNormal())
                {
                    swap(ns[0],ns[1]);
                }

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
                Array<OneD, NekDouble> bestN = GetNormal(g->second);

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

                            NekDouble valt = Visability(g->second,tmp);

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

                bln.N = bestN;

                if(Visability(g->second,bln.N) < 0.0)
                {
                    cout << "failed " << (*it)->m_x << " " << (*it)->m_y << " "
                                      << (*it)->m_z << endl;
                    failed++;
                }

            }

            Array<OneD, NekDouble> loc = (*it)->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bln.N[k] * bln.bl;
            }
            bln.pNode = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                                            loc[0], loc[1], loc[2]));


            //if the diff size is greater than 1 there is a curve that needs remeshing
            if(diff.size() > 1)
            {
                cout << diff.size() << endl;
            }
            ASSERTL0(diff.size() <= 1,"not setup for curve bl refinement");
            for(int i = 0; i < diff.size(); i++)
            {
                symSurfs.insert(diff[i]);
                bln.symsurf = diff[i];
            }
            blData[(*it)] = bln;
        }
    }

    if(failed > 0)
    {
        cout << "some normals are not visible" << endl;
        abort();
    }

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

    map<NodeSharedPtr, blInfo>::iterator bit;

    //make prisms

    map<int,int> nm;
    map<ElementSharedPtr,ElementSharedPtr> priToTri;

    ElmtConfig pconf(LibUtilities::ePrism,1,false,false);
    ElmtConfig tconf(LibUtilities::eTriangle,1,false,false);

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_mesh->m_element[2].size(), "Building BL elements\t");
        }

        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int> s;
        s.push_back(el->CADSurfId);

        vector<unsigned int> inter;

        set_intersection(m_blsurfs.begin(), m_blsurfs.end(),
                         s.begin(), s.end(),
                         back_inserter(inter));

        if(inter.size() == 0)
        {
            //triangle is not on a boundary layer surface
            continue;
        }

        vector<NodeSharedPtr> tn(3); //nodes for pseduo surface
        vector<NodeSharedPtr> pn(6); //all prism nodes
        vector<NodeSharedPtr> n = el->GetVertexList();

        CADSurfSharedPtr tmps = m_cad->GetSurf(el->CADSurfId);

        if(tmps->IsReversedNormal())
        {
            nm[0] = 0;
            nm[1] = 3;
            nm[2] = 1;
            nm[3] = 2;
            nm[4] = 4;
            nm[5] = 5;
        }
        else
        {
            nm[0] = 0;
            nm[1] = 3;
            nm[2] = 4;
            nm[3] = 5;
            nm[4] = 1;
            nm[5] = 2;
        }

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

        m_mesh->m_element[3].push_back(E);

        //tag of this element doesnt matter so can just be 1
        ElementSharedPtr T = GetElementFactory().
                    CreateInstance(LibUtilities::eTriangle, tconf, tn, tags);
        m_psuedoSurface.push_back(T);

        priToTri[E] = el;
    }

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

    //this is where it should do some clever collision dectecting and reduce the bl parameter

    //all nodes in the vertex set are unique ordered and in surface elements
    //so these will form the dataset for ANN
    ANNpointArray dataPts;
    ANNpoint queryPt;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    ANNkd_tree* kdTree;
    queryPt = annAllocPt(3);
    dataPts = annAllocPts(m_mesh->m_vertexSet.size(), 3);

    for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++)
    {
        dataPts[(*it)->m_id][0] = (*it)->m_x;
        dataPts[(*it)->m_id][1] = (*it)->m_y;
        dataPts[(*it)->m_id][2] = (*it)->m_z;
    }

    kdTree = new ANNkd_tree(dataPts, m_mesh->m_vertexSet.size(), 3);

    ct = 0;
    for(bit = blData.begin(); bit != blData.end(); bit++, ct++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                ct, blData.size(), "Proximity sweeping\t");
        }

        queryPt[0] = bit->first->m_x;
        queryPt[1] = bit->first->m_y;
        queryPt[2] = bit->first->m_z;

        //do get a decent data set must increase the sample set until the last
        //point is further away than the bl (i.e no possible intersection)
        //set an inital sample size of 200 and keep doubling

        int sample = 100;

        do
        {
            sample *= 2;
            nnIdx = new ANNidx[sample];
            dists = new ANNdist[sample];
            if(sample > m_mesh->m_vertexSet.size())
            {
                sample = m_mesh->m_vertexSet.size();
            }
            kdTree->annkSearch(queryPt, sample, nnIdx, dists);
            if(sample == m_mesh->m_vertexSet.size())
            {
                break;
            }
        }
        while(sqrt(dists[sample-1]) < bit->second.bl * 2.5);

        //now need to build a set of triagnles to test against
        //use set to make sure its unique
        set<ElementSharedPtr> tris;
        for(int i = 0; i < sample; i++)
        {
            NekDouble t = (dataPts[nnIdx[i]][0] - bit->first->m_x) * bit->second.N[0] +
                          (dataPts[nnIdx[i]][1] - bit->first->m_y) * bit->second.N[1] +
                          (dataPts[nnIdx[i]][2] - bit->first->m_z) * bit->second.N[2];
            t /= sqrt((dataPts[nnIdx[i]][0] - bit->first->m_x)*(dataPts[nnIdx[i]][0] - bit->first->m_x)+
                      (dataPts[nnIdx[i]][1] - bit->first->m_y)*(dataPts[nnIdx[i]][1] - bit->first->m_y)+
                      (dataPts[nnIdx[i]][2] - bit->first->m_z)*(dataPts[nnIdx[i]][2] - bit->first->m_z));
            t /= sqrt(bit->second.N[0]*bit->second.N[0] + bit->second.N[1]*bit->second.N[1] +
                      bit->second.N[2]*bit->second.N[2]);
            if(t < 0.5)
            {
                continue;
            }

            if(dists[i] < bit->second.bl * 2.5)
            {
                map<int, vector<ElementSharedPtr> >::iterator s = nIdxToTri.find(nnIdx[i]);
                ASSERTL0(s != nIdxToTri.end(),"not found");

                for(int j = 0; j < s->second.size(); j++)
                {
                    tris.insert(s->second[j]);
                }
            }
        }

        //make it so that it checks a number of normals
        //the averaged Visability and others from the triangles

        vector<Array<OneD, NekDouble> > norms;

        norms.push_back(bit->second.N);

        map<int, vector<ElementSharedPtr> >::iterator g = nIdxToTri.find(bit->first->m_id);
        for(int i = 0; i < g->second.size(); i++)
        {
            vector<NodeSharedPtr> ns = g->second[i]->GetVertexList();
            if(m_cad->GetSurf(g->second[i]->CADSurfId)->IsReversedNormal())
            {
                swap(ns[0],ns[1]);
            }

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
            norms.push_back(tmp);
        }

        NekDouble mind = numeric_limits<double>::max();
        set<ElementSharedPtr>::iterator s;

        for(s = tris.begin(); s != tris.end(); s++)
        {
            ElementSharedPtr el = (*s);
            vector<NodeSharedPtr> ns = el->GetVertexList();

            for(int i = 0 ; i < norms.size(); i++)
            {
                DNekMat A(3,3,0.0);
                DNekMat B(3,1,0.0);
                A(0,0) = norms[i][0] * -1.0;
                A(1,0) = norms[i][1] * -1.0;
                A(2,0) = norms[i][2] * -1.0;
                A(0,1) = ns[1]->m_x - ns[0]->m_x;
                A(1,1) = ns[1]->m_y - ns[0]->m_y;
                A(2,1) = ns[1]->m_z - ns[0]->m_z;
                A(0,2) = ns[2]->m_x - ns[0]->m_x;
                A(1,2) = ns[2]->m_y - ns[0]->m_y;
                A(2,2) = ns[2]->m_z - ns[0]->m_z;

                NekDouble det = A(0,0) * (A(1,1)*A(2,2) - A(2,1)*A(1,2))
                               -A(0,1) * (A(1,0)*A(2,2) - A(2,0)*A(1,2))
                               +A(0,2) * (A(1,0)*A(2,1) - A(2,0)*A(1,1));
                if(fabs(det) < 1e-12)
                {
                    //no intersecton
                    continue;
                }
                B(0,0) = bit->first->m_x - ns[0]->m_x;
                B(1,0) = bit->first->m_y - ns[0]->m_y;
                B(2,0) = bit->first->m_z - ns[0]->m_z;

                A.Invert();

                DNekMat X = A * B; //t u v

                if(X(0,0) < 1e-6 || X(0,0) > bit->second.bl * 2.5)
                {
                    //no plane intersecton possible
                    continue;
                }
                //check triangle intersecton
                if(X(1,0) >= 0.0 && X(2,0) >= 0.0 && X(1,0) + X(2,0) <= 1.0)
                {
                    //hit
                    NekDouble tmp = X(0,0);
                    mind = min(mind, tmp);
                }
            }

        }

        if(mind < bit->second.bl * 2.5 && mind * 0.25 < bit->second.bl)
        {
            bit->second.bl = mind * 0.25;
            bit->second.pNode->m_x = bit->first->m_x + bit->second.N[0] * bit->second.bl;
            bit->second.pNode->m_y = bit->first->m_y + bit->second.N[1] * bit->second.bl;
            bit->second.pNode->m_z = bit->first->m_z + bit->second.N[2] * bit->second.bl;
        }
    }

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

    //loop over all prisms, if invalid shrink until it is
    //being careful to act on nodes which have already been shrunk
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
                cout << ns[0]->m_id << endl;
                cout << ns[1]->m_id << endl;
                cout << ns[2]->m_id << endl;
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

            }

            geom = el->GetGeom(m_mesh->m_spaceDim);
            gfac = geom->GetGeomFactors();
        }
    }

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

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
                    repeat = true;
                }
            }
        }
    }

    /*for(int i = 0; i < m_psuedoSurface.size(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_psuedoSurface.size(), "Intersection sweeping");
        }

        for(int j = i+1; j < m_psuedoSurface.size(); j++)
        {
            vector<NodeSharedPtr> V1 = m_psuedoSurface[i]->GetVertexList();
            vector<NodeSharedPtr> V2 = m_psuedoSurface[j]->GetVertexList();

            Array<OneD, NekDouble> N1(3,0.0);
            N1[0] = (V1[1]->m_y - V1[0]->m_y) * (V1[2]->m_z - V1[0]->m_z) -
                    (V1[1]->m_z - V1[0]->m_z) * (V1[2]->m_y - V1[0]->m_y);
            N1[1] = (V1[1]->m_z - V1[0]->m_z) * (V1[2]->m_x - V1[0]->m_x) -
                    (V1[1]->m_x - V1[0]->m_x) * (V1[2]->m_z - V1[0]->m_z);
            N1[2] = (V1[1]->m_x - V1[0]->m_x) * (V1[2]->m_y - V1[0]->m_y) -
                    (V1[1]->m_y - V1[0]->m_y) * (V1[2]->m_x - V1[0]->m_x);

            Array<OneD, NekDouble> N2(3,0.0);
            N2[0] = (V2[1]->m_y - V2[0]->m_y) * (V2[2]->m_z - V2[0]->m_z) -
                    (V2[1]->m_z - V2[0]->m_z) * (V2[2]->m_y - V2[0]->m_y);
            N2[1] = (V2[1]->m_z - V2[0]->m_z) * (V2[2]->m_x - V2[0]->m_x) -
                    (V2[1]->m_x - V2[0]->m_x) * (V2[2]->m_z - V2[0]->m_z);
            N2[2] = (V2[1]->m_x - V2[0]->m_x) * (V2[2]->m_y - V2[0]->m_y) -
                    (V2[1]->m_y - V2[0]->m_y) * (V2[2]->m_x - V2[0]->m_x);

            NekDouble d2 = N2[0] * V2[0]->m_x + N2[1] * V2[0]->m_y + N2[2] * V2[0]->m_z;
            NekDouble d1 = N1[0] * V1[0]->m_x + N1[1] * V1[0]->m_y + N1[2] * V1[0]->m_z;

            NekDouble dV10 = N2[0] * V1[0]->m_x + N2[1] * V1[0]->m_y + N2[2] * V1[0]->m_z + d2;
            NekDouble dV11 = N2[0] * V1[1]->m_x + N2[1] * V1[1]->m_y + N2[2] * V1[1]->m_z + d2;
            NekDouble dV12 = N2[0] * V1[2]->m_x + N2[1] * V1[2]->m_y + N2[2] * V1[2]->m_z + d2;

            NekDouble dV20 = N1[0] * V2[0]->m_x + N1[1] * V2[0]->m_y + N1[2] * V2[0]->m_z + d1;
            NekDouble dV21 = N1[0] * V2[1]->m_x + N1[1] * V2[1]->m_y + N1[2] * V2[1]->m_z + d1;
            NekDouble dV22 = N1[0] * V2[2]->m_x + N1[1] * V2[2]->m_y + N1[2] * V2[2]->m_z + d1;

            if(dV10 < 0.0 && dV11 < 0.0 && dV11 < 0.0 ||
               dV10 >= 0.0 && dV11 >= 0.0 && dV11 >= 0.0)
            {
                continue;
            }
            if(dV20 < 0.0 && dV21 < 0.0 && dV21 < 0.0 ||
               dV20 >= 0.0 && dV21 >= 0.0 && dV21 >= 0.0)
            {
                continue;
            }

            Array<OneD, NekDouble> D(3,0.0);
            D[0] = N1[1] * N2[2] -
                   N1[2] * N2[1];
            D[1] = N1[2] * N2[0] -
                   N1[0] * N2[2];
            D[2] = N1[0] * N2[1] -
                   N1[1] * N2[0];

            NekDouble pV10 = D[0] * V1[0]->m_x + D[1] * V1[0]->m_y + D[2] * V1[0]->m_z;
            NekDouble pV11 = D[0] * V1[1]->m_x + D[1] * V1[1]->m_y + D[2] * V1[1]->m_z;
            NekDouble pV12 = D[0] * V1[2]->m_x + D[1] * V1[2]->m_y + D[2] * V1[2]->m_z;

            NekDouble pV20 = D[0] * V2[0]->m_x + D[1] * V2[0]->m_y + D[2] * V2[0]->m_z;
            NekDouble pV21 = D[0] * V2[1]->m_x + D[1] * V2[1]->m_y + D[2] * V2[1]->m_z;
            NekDouble pV22 = D[0] * V2[2]->m_x + D[1] * V2[2]->m_y + D[2] * V2[2]->m_z;

            NekDouble t11, t22, t12, t21;

            if(dV10 > 0.0 && dV11 <= 0.0 && dV12 <= 0.0 ||
               dV10 < 0.0 && dV11 >= 0.0 && dV12 >= 0.0)
            {
                t11 = pV10 + (pV12 - pV10) * dV10 / (dV10 - dV12);
                t12 = pV10 + (pV11 - pV10) * dV10 / (dV10 - dV11);
            }
            else if (dV11 > 0.0 && dV10 <= 0.0 && dV12 <= 0.0 ||
                     dV11 < 0.0 && dV10 >= 0.0 && dV12 >= 0.0)
            {
                t11 = pV11 + (pV10 - pV11) * dV11 / (dV11 - dV10);
                t12 = pV11 + (pV12 - pV11) * dV11 / (dV11 - dV12);
            }
            else if (dV12 > 0.0 && dV11 <= 0.0 && dV10 <= 0.0 ||
                     dV12 < 0.0 && dV11 >= 0.0 && dV10 >= 0.0)
            {
                t11 = pV12 + (pV11 - pV12) * dV12 / (dV12 - dV11);
                t12 = pV12 + (pV10 - pV12) * dV12 / (dV12 - dV10);
            }
            else
            {
                ASSERTL0(false,"failed to find tri orient");
            }

            if(dV20 > 0.0 && dV21 <= 0.0 && dV22 <= 0.0 ||
               dV20 < 0.0 && dV21 >= 0.0 && dV22 >= 0.0)
            {
                t21 = pV20 + (pV22 - pV20) * dV20 / (dV20 - dV22);
                t22 = pV20 + (pV21 - pV20) * dV20 / (dV20 - dV21);
            }
            else if (dV21 > 0.0 && dV20 <= 0.0 && dV22 <= 0.0 ||
                     dV21 < 0.0 && dV20 >= 0.0 && dV22 >= 0.0)
            {
                t21 = pV21 + (pV20 - pV21) * dV21 / (dV21 - dV20);
                t22 = pV21 + (pV22 - pV21) * dV21 / (dV21 - dV22);
            }
            else if (dV22 > 0.0 && dV21 <= 0.0 && dV20 <= 0.0 ||
                     dV22 < 0.0 && dV21 >= 0.0 && dV20 >= 0.0)
            {
                t21 = pV22 + (pV21 - pV22) * dV22 / (dV22 - dV21);
                t22 = pV22 + (pV20 - pV22) * dV22 / (dV22 - dV20);
            }
            else
            {
                ASSERTL0(false,"failed to find tri orient");
            }

            if(t11 > t12)
            {
                swap(t11,t12);
            }

            if(t21 > t22)
            {
                swap(t21,t22);
            }

            if(t21 < t12 && t21 > t11)
            {
                cout << endl << "intersection" << endl;
            }
        }
    }

    exit(-1);*/

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
                cout << ns[0]->m_id << endl;
                cout << ns[1]->m_id << endl;
                cout << ns[2]->m_id << endl;
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

            }

            geom = el->GetGeom(m_mesh->m_spaceDim);
            gfac = geom->GetGeomFactors();
        }
    }

    repeat = true;
    while(repeat)
    {
        repeat = false;
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
                    repeat = true;
                }
            }
        }
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

}
}
}
