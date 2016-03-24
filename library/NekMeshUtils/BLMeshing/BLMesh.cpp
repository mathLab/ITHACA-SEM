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

    //need a map from vertex idx to surface elements
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
            bln.bl = m_bl;
            bln.symsurf = 0;

            bln.N = Array<OneD, NekDouble> (3,0.0);
            for(int j = 0; j < inter.size(); j++)
            {
                Array<OneD, NekDouble> uv = (*it)->GetCADSurfInfo(inter[j]);
                Array<OneD, NekDouble> N = m_cad->GetSurf(inter[j])->N(uv);
                for(int k = 0; k < 3; k++)
                {
                    bln.N[k] += N[k];
                }
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
            Array<OneD, NekDouble> loc = (*it)->GetLoc();
            for(int k = 0; k < 3; k++)
            {
                loc[k] += bln.N[k] * bln.bl;
            }
            bln.pNode = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                                            loc[0], loc[1], loc[2]));

            //calculate mesh normal
            Array<OneD, NekDouble> mNorm(3,0.0);
            map<int, vector<ElementSharedPtr> >::iterator g = nIdxToTri.find((*it)->m_id);
            for(int i = 0; i < g->second.size(); i++)
            {
                vector<NodeSharedPtr> ns = g->second[i]->GetVertexList();
                if(m_cad->GetSurf(g->second[i]->CADSurfId)->IsReversedNormal())
                {
                    swap(ns[0],ns[1]);
                }
                mNorm[0] += ((ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                             (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y));
                mNorm[1] -= ((ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z) -
                             (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x));
                mNorm[2] += ((ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                             (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x));
            }
            mag = 0.0;
            for(int k = 0; k < 3; k++)
            {
                mag += mNorm[k]*mNorm[k];
            }
            mag = sqrt(mag);
            for(int k = 0; k < 3; k++)
            {
                mNorm[k] /= mag;
            }

            if(mNorm[0] * bln.N[0] + mNorm[1] * bln.N[1] + mNorm[2] * bln.N[2] < 0.9)
            {
                cout << "Norm irregularity ";
                cout << mNorm[0] * bln.N[0] + mNorm[1] * bln.N[1] + mNorm[2] * bln.N[2];
                if(inter.size() ==3)
                {
                    cout << " with 3 normal";
                }
                cout << " " << (*it)->m_id << endl;
            }

            //bln.N = mNorm;

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


    map<NodeSharedPtr, blInfo>::iterator bit;

    //make prisms

    map<int,int> nm;
    map<ElementSharedPtr,ElementSharedPtr> priToTri;

    ElmtConfig pconf(LibUtilities::ePrism,1,false,false);
    ElmtConfig tconf(LibUtilities::eTriangle,1,false,false);

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
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

    //loop over all prisms, if invalid shrink until it is
    //being careful to act on nodes which have already been shrunk
/*    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
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
*/
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
        while(sqrt(dists[sample-1]) < bit->second.bl * 1.1);

        //now need to build a set of triagnles to test against
        //use set to make sure its unique
        set<ElementSharedPtr> tris;
        for(int i = 0; i < sample; i++)
        {
            map<int, vector<ElementSharedPtr> >::iterator s = nIdxToTri.find(nnIdx[i]);
            ASSERTL0(s != nIdxToTri.end(),"not found");

            for(int j = 0; j < s->second.size(); j++)
            {
                tris.insert(s->second[j]);
            }
        }

        NekDouble mind = numeric_limits<double>::max();
        set<ElementSharedPtr>::iterator s;
        for(s = tris.begin(); s != tris.end(); s++)
        {
            ElementSharedPtr el = (*s);
            vector<NodeSharedPtr> ns = el->GetVertexList();
            if(m_cad->GetSurf(el->CADSurfId)->IsReversedNormal())
            {
                swap(ns[0], ns[1]);
            }

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

            if(fabs(nu) < 1E-6)
            {
                //no intersection
                continue;
            }
            else if(nu > 0)
            {
                continue;
            }

            NekDouble d = ns[0]->m_x * norm[0] +
                          ns[0]->m_y * norm[1] +
                          ns[0]->m_z * norm[2];

            NekDouble t =  d - bit->first->m_x*norm[0] -
                               bit->first->m_y*norm[1] -
                               bit->first->m_z*norm[2];
            t /= nu;
            if(t < 1E-6 || t > bit->second.bl * 1.1)
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
            //quick test it is on surface
            NekDouble tst = x[0] * norm[0] + x[1] * norm[1] + x[2] * norm[2];
            ASSERTL0(fabs(tst - d) < 1e-6,"failed planar test");

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

            bool tst1 = (c1[0]*norm[0] + c1[1]*norm[1] + c1[2]* norm[2] >= 0.0);
            bool tst2 = (c2[0]*norm[0] + c2[1]*norm[1] + c2[2]* norm[2] >= 0.0);
            bool tst3 = (c3[0]*norm[0] + c3[1]*norm[1] + c3[2]* norm[2] >= 0.0);

            if(tst1 && tst2 && tst3)
            {
                //hit ?
                mind = min(mind, t);
            }
        }
        if(mind < bit->second.bl * 1.1)
        {
            bit->second.bl = mind * 0.25;
            bit->second.pNode->m_x = bit->first->m_x + bit->second.N[0] * bit->second.bl;
            bit->second.pNode->m_y = bit->first->m_y + bit->second.N[1] * bit->second.bl;
            bit->second.pNode->m_z = bit->first->m_z + bit->second.N[2] * bit->second.bl;
        }
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
