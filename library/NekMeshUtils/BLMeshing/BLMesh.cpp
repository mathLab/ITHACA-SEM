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

#include <NekMeshUtils/BLMeshing/BLMesh.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

#include <ANN/ANN.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void BLMesh::Mesh()
{
    //At this stage the surface mesh is complete and the elements know their
    //neigbours through element links in the edges,,

    // here elements are made for the boundary layer they will need to know
    // links (maybe facelinks), so that the tetmeshing module can extract the
    // surface upon which it needs to mesh (top of the bl and the rest of the
    // surface).

    set<int> symSurfs;

    NodeSet::iterator it;
    for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++)
    {
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
            bln.surfs = inter;
            bln.bl = m_bl;
            bln.symsurf = 0;

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

    map<int, vector<NodeSharedPtr> > fcnodes;
    vector<ElementSharedPtr> pTri;
    int ct = 0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<unsigned int> s;
        s.push_back(m_mesh->m_element[2][i]->CADSurfId);
        vector<unsigned int> inter;
        set_intersection(m_blsurfs.begin(), m_blsurfs.end(),
                         s.begin(), s.end(),
                         back_inserter(inter));
        if(inter.size() > 0)
        {
            pTri.push_back(m_mesh->m_element[2][i]);
            vector<NodeSharedPtr> fcnode = m_mesh->m_element[2][i]->GetVertexList();
            fcnodes[ct++] = fcnode;
        }
    }

    map<NodeSharedPtr, blInfo>::iterator bit;
    //compile normals
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        bit->second.N = Array<OneD, NekDouble> (3,0.0);
        for(int j = 0; j < bit->second.surfs.size(); j++)
        {
            Array<OneD, NekDouble> uv = bit->first->GetCADSurfInfo(bit->second.surfs[j]);
            Array<OneD, NekDouble> N = m_cad->GetSurf(bit->second.surfs[j])->N(uv);
            for(int k = 0; k < 3; k++)
            {
                bit->second.N[k] += N[k];
            }
        }
        NekDouble mag = 0.0;
        for(int k = 0; k < 3; k++)
        {
            mag += bit->second.N[k]*bit->second.N[k];
        }
        mag = sqrt(mag);
        for(int k = 0; k < 3; k++)
        {
            bit->second.N[k] /= mag;
        }
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

    NodeSet::iterator nit;
    for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
    {
        dataPts[(*nit)->m_id][0] = (*nit)->m_x;
        dataPts[(*nit)->m_id][1] = (*nit)->m_y;
        dataPts[(*nit)->m_id][2] = (*nit)->m_z;
    }

    //need a map from vertex idx to surface elements
    map<int, vector<int> > nIdxToTriIdx;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            nIdxToTriIdx[ns[j]->m_id].push_back(m_mesh->m_element[2][i]->GetId());
        }
    }

    kdTree = new ANNkd_tree(dataPts, m_mesh->m_vertexSet.size(), 3);

    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        queryPt[0] = bit->first->m_x;
        queryPt[1] = bit->first->m_y;
        queryPt[2] = bit->first->m_z;

        //do get a decent data set must increase the sample set until the last
        //point is further away than the bl (i.e no possible intersection)
        //set an inital sample size of 50 and keep doubling

        int sample = 50;

        do
        {
            sample *= 2;
            nnIdx = new ANNidx[sample];
            dists = new ANNdist[sample];
            kdTree->annkSearch(queryPt, sample, nnIdx, dists);
        }
        while(sqrt(dists[sample-1]) < bit->second.bl * 2.5);

        //now need to build a set of triagnles to test against
        //use set to make sure its unique
        map<int, vector<int> >::iterator o = nIdxToTriIdx.find(bit->first->m_id);
        ASSERTL0(o != nIdxToTriIdx.end(),"not found");
        set<int> tris;
        for(int i = 0; i < sample; i++)
        {
            map<int, vector<int> >::iterator s = nIdxToTriIdx.find(nnIdx[i]);
            ASSERTL0(s != nIdxToTriIdx.end(),"not found");
            vector<int> inter;
            set_intersection(o->second.begin(), o->second.end(),
                             s->second.begin(), s->second.end(),
                             back_inserter(inter));
            for(int j = 0; j < s->second.size(); j++)
            {
                bool add = true;
                for(int k = 0; k < inter.size(); k++)
                {
                    if(s->second[j] == inter[k])
                    {
                        add = false;
                        break;
                    }
                }
                if(add)
                {
                    tris.insert(s->second[j]);
                }
            }
        }

        NekDouble mind = bit->second.bl * 2.5;
        set<int>::iterator s;
        for(s = tris.begin(); s != tris.end(); s++)
        {
            ElementSharedPtr el = m_mesh->m_element[2][(*s)];
            vector<NodeSharedPtr> ns = el->GetVertexList();
            Array<OneD, NekDouble> norm(3,0.0);
            norm[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                      (ns[2]->m_y - ns[0]->m_y) * (ns[1]->m_z - ns[0]->m_z);
            norm[1] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z) -
                      (ns[2]->m_x - ns[0]->m_x) * (ns[1]->m_z - ns[0]->m_z);
            norm[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                      (ns[2]->m_x - ns[0]->m_x) * (ns[1]->m_y - ns[0]->m_y);
            NekDouble mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
            norm[0] /= mag;
            norm[1] /= mag;
            norm[2] /= mag;

            NekDouble nu = norm[0] * bit->second.N[0] +
                           norm[1] * bit->second.N[1] +
                           norm[2] * bit->second.N[2];
            if(nu < 1E-15)
            {
                //no intersection
                continue;
            }

            NekDouble t = (bit->first->m_x - ns[0]->m_x)*norm[0] +
                          (bit->first->m_y - ns[0]->m_y)*norm[1] +
                          (bit->first->m_z - ns[0]->m_z)*norm[2];
            t *= -1.0 / nu;
            if(t < 1E-6 || t > bit->second.bl * 2.5)
            {
                //no intersection worth worrying about
                continue;
            }

            //so by this point there is an intersecton with the plane of the triangle
            //within the range of the boudnary layer
            //need to determine if it acutally hits the triangle
            Array<OneD, NekDouble> x(3,0.0);
            x[0] = bit->first->m_x + bit->second.N[0] * t;
            x[1] = bit->first->m_y + bit->second.N[1] * t;
            x[2] = bit->first->m_z + bit->second.N[2] * t;
            Array<OneD, NekDouble> c1(3,0.0), c2(3,0.0), c3(3,0.0);
            c1[0] = (ns[1]->m_y - ns[0]->m_y) * (x[2] - ns[0]->m_z) -
                    (x[1] - ns[0]->m_y) * (ns[1]->m_z - ns[0]->m_z);
            c1[1] = (ns[1]->m_x - ns[0]->m_x) * (x[2] - ns[0]->m_z) -
                    (x[0] - ns[0]->m_x) * (ns[1]->m_z - ns[0]->m_z);
            c1[2] = (ns[1]->m_x - ns[0]->m_x) * (x[1] - ns[0]->m_y) -
                    (x[0] - ns[0]->m_x) * (ns[1]->m_y - ns[0]->m_y);
            c2[0] = (ns[2]->m_y - ns[1]->m_y) * (x[2] - ns[1]->m_z) -
                    (x[1] - ns[1]->m_y) * (ns[2]->m_z - ns[1]->m_z);
            c2[1] = (ns[2]->m_x - ns[1]->m_x) * (x[2] - ns[1]->m_z) -
                    (x[0] - ns[1]->m_x) * (ns[2]->m_z - ns[1]->m_z);
            c2[2] = (ns[2]->m_x - ns[1]->m_x) * (x[1] - ns[1]->m_y) -
                    (x[0] - ns[1]->m_x) * (ns[2]->m_y - ns[1]->m_y);
            c3[0] = (ns[0]->m_y - ns[2]->m_y) * (x[2] - ns[2]->m_z) -
                    (x[1] - ns[2]->m_y) * (ns[0]->m_z - ns[2]->m_z);
            c3[1] = (ns[0]->m_x - ns[2]->m_x) * (x[2] - ns[2]->m_z) -
                    (x[0] - ns[2]->m_x) * (ns[0]->m_z - ns[2]->m_z);
            c3[2] = (ns[0]->m_x - ns[2]->m_x) * (x[1] - ns[2]->m_y) -
                    (x[0] - ns[2]->m_x) * (ns[0]->m_y - ns[2]->m_y);

            bool tst1 = c1[0]*norm[0] + c1[1]*norm[1] + c1[2]* norm[2] >= 0.0;
            bool tst2 = c2[0]*norm[0] + c2[1]*norm[1] + c2[2]* norm[2] >= 0.0;
            bool tst3 = c3[0]*norm[0] + c3[1]*norm[1] + c3[2]* norm[2] >= 0.0;

            if(tst1 && tst2 && tst3)
            {
                //hit ?
                mind = min(mind, t);
            }
        }
        if(mind < bit->second.bl * 2.5)
        {
            bit->second.bl = mind * 0.25;
        }
    }

    //build new nodes
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        Array<OneD, NekDouble> loc = bit->first->GetLoc();
        for(int k = 0; k < 3; k++)
        {
            loc[k] += bit->second.N[k] * bit->second.bl;
        }
        bit->second.pNode = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                                        loc[0], loc[1], loc[2]));
    }

    //make prisms

    map<int,int> nm;

    ElmtConfig pconf(LibUtilities::ePrism,1,false,false);
    ElmtConfig tconf(LibUtilities::eTriangle,1,false,false);

    for(int i = 0; i < pTri.size(); i++)
    {
        vector<NodeSharedPtr> tn(3); //nodes for pseduo surface
        vector<NodeSharedPtr> pn(6); //all prism nodes
        vector<NodeSharedPtr> n = pTri[i]->GetVertexList();

        vector<pair<int, CADSurfSharedPtr> > tmpss = n[0]->GetCADSurfs();
        CADSurfSharedPtr tmps = m_cad->GetSurf(pTri[i]->CADSurfId);

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
    }

    //loop over all prisms, if invalid shrink until it is
    //being careful to act on nodes which have already been shrunk
    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[3][i];
        SpatialDomains::GeometrySharedPtr geom =
                                            el->GetGeom(m_mesh->m_spaceDim);
        SpatialDomains::GeomFactorsSharedPtr gfac =
                                            geom->GetGeomFactors();

        map<int, vector<NodeSharedPtr> >::iterator fit = fcnodes.find(i);
        ASSERTL0(fit != fcnodes.end(), "not found");
        vector<NodeSharedPtr> ns = fit->second;

        while(!gfac->IsValid())
        {
            NekDouble maxbl = max(blData[ns[0]].bl, blData[ns[1]].bl);
            maxbl = max(maxbl, blData[ns[2]].bl);

            if(maxbl < 1E-6)
            {
                cout << "shrunk too far" << endl;
                break;
            }

            int ct = 0;
            for(int j = 0; j < 3; j++)
            {
                map<NodeSharedPtr, blInfo>::iterator bli = blData.find(ns[j]);
                ASSERTL0(bli != blData.end(), "not found");
                if(bli->second.bl <= maxbl * 0.5)
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


    //smoothing
    //need to build a list of nodes to neigbours
    map<NodeSharedPtr, NodeSet> nodeToNear;
    for(int i = 0; i < pTri.size(); i++)
    {
        vector<EdgeSharedPtr> es = pTri[i]->GetEdgeList();
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
