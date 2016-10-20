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

#include <NekMeshUtils/VolumeMeshing/BLMeshing/BLMesh.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

#include <algorithm>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, unsigned int> boxI;

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

inline box GetBox(ElementSharedPtr el)
{
    NekDouble xmin =        numeric_limits<double>::max(),
              xmax = -1.0 * numeric_limits<double>::max(),
              ymin =        numeric_limits<double>::max(),
              ymax = -1.0 * numeric_limits<double>::max(),
              zmin =        numeric_limits<double>::max(),
              zmax = -1.0 * numeric_limits<double>::max();

    vector<NodeSharedPtr> ns = el->GetVertexList();
    for(int i = 0; i < ns.size(); i++)
    {
        xmin = min(xmin,ns[i]->m_x);
        xmax = max(xmax,ns[i]->m_x);
        ymin = min(ymin,ns[i]->m_y);
        ymax = max(ymax,ns[i]->m_y);
        zmin = min(zmin,ns[i]->m_z);
        zmax = max(zmax,ns[i]->m_z);
    }

    return box(point(xmin,ymin,zmin),point(xmax,ymax,zmax));
}

void BLMesh::Mesh()
{
    //At this stage the surface mesh is complete and the elements know their
    //neigbours through element links in the edges,,

    // here elements are made for the boundary layer they will need to know
    // links (maybe facelinks), so that the tetmeshing module can extract the
    // surface upon which it needs to mesh (top of the bl and the rest of the
    // surface).

    NekDouble a = 2.0 * (1.0 - m_prog) / (1.0 - pow(m_prog,m_layer+1));
    Array<OneD, NekDouble> layerT(m_layer);
    layerT[0] = a * m_prog * m_bl;
    for(int i = 1; i < m_layer; i++)
    {
        layerT[i] = layerT[i-1] + a * pow(m_prog,i) * m_bl;
    }

    cout << "First layer height " << layerT[0] << endl;

    //this sets up all the boundary layer normals data holder
    set<int> symSurfs;
    NodeSet::iterator it;
    int ct = 0;
    int failed = 0;

    ofstream file1;
    file1.open("pts.3D");
    file1 << "X Y Z value" << endl;
    for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++, ct++)
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
            blInfoSharedPtr bln = boost::shared_ptr<blInfo>(new blInfo);
            bln->oNode = (*it);
            bln->stopped = false;
            bln->stop = false;

            file1 << (*it)->m_x << " " << (*it)->m_y << " " << (*it)->m_z << " " << ss.size() << endl;

            if(diff.size() > 0)
            {
                //if the diff size is greater than 1 there is a curve that needs remeshing
                ASSERTL0(diff.size() <= 1,"not setup for curve bl refinement");
                symSurfs.insert(diff[0]);
                bln->symsurf = diff[0];
                bln->onSym = true;
            }
            else
            {
                bln->onSym = false;
            }

            blData[(*it)] = bln;
        }
    }
    file1.close();

    //need a map from vertex idx to surface elements
    //but do not care about triangles which are not in the bl
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        //orientate the triangle
        if(m_mesh->m_cad->GetSurf(m_mesh->m_element[2][i]->CADSurfId)
                                                        ->IsReversedNormal())
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
            blData[ns[j]]->els.push_back(m_mesh->m_element[2][i]);
        }
    }

    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        //calculate mesh normal
        bit->second->N = GetNormal(bit->second->els);

        if(bit->second->N[0] != bit->second->N[0])
        {
            cout << "nan in normal" << endl;
            exit(-1);
        }

        if(Visability(bit->second->els,bit->second->N) < 0.0)
        {
            cerr << "failed " << bit->first->m_x << " " << bit->first->m_y << " "
                              << bit->first->m_z << " "
                              << Visability(bit->second->els,bit->second->N) << endl;
            failed++;
        }

        Array<OneD, NekDouble> loc = bit->first->GetLoc();
        for(int k = 0; k < 3; k++)
        {
            loc[k] += bit->second->N[k] * layerT[0];
        }

        bit->second->pNode = boost::shared_ptr<Node>(new Node(
                                        m_mesh->m_numNodes++,
                                        loc[0], loc[1], loc[2]));
        bit->second->bl = 0;
    }

    m_symSurfs = vector<unsigned int>(symSurfs.begin(), symSurfs.end());

    //now need to enforce that all symmetry plane nodes have their normal
    //forced onto the symmetry surface
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        if(!bit->second->onSym)
        {
            continue;
        }

        Array<OneD, NekDouble> uv(2);
        Array<OneD, NekDouble> loc = bit->second->pNode->GetLoc();
        m_mesh->m_cad->GetSurf(bit->second->symsurf)->ProjectTo(loc, uv);

        Array<OneD, NekDouble> nl = m_mesh->m_cad->
                                        GetSurf(bit->second->symsurf)->P(uv);

        Array<OneD, NekDouble> N(3);
        N[0] = nl[0] - bit->first->m_x;
        N[1] = nl[1] - bit->first->m_y;
        N[2] = nl[2] - bit->first->m_z;

        NekDouble mag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
        N[0] /= mag;
        N[1] /= mag;
        N[2] /= mag;

        bit->second->N = N;
        bit->second->AlignNode(layerT[0]);
    }


    //now smooth all the normals by distance weighted average
    //keep normals on curves constant
    map<NodeSharedPtr, vector<blInfoSharedPtr> > nToNInfo; //node to neighbouring information
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        set<int> added;
        added.insert(bit->first->m_id);
        for(int i = 0; i < bit->second->els.size(); i++)
        {
            vector<NodeSharedPtr> ns = bit->second->els[i]->GetVertexList();
            for(int j = 0; j < ns.size(); j++)
            {
                set<int>::iterator t = added.find(ns[j]->m_id);
                if(t == added.end())
                {
                    nToNInfo[bit->first].push_back(blData[ns[j]]);
                }
            }
        }
    }

    for(int l = 0; l < 10; l++)
    {
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(bit->first->GetNumCADSurf() > 1)
            {
                continue;
            }

            Array<OneD, NekDouble> sumV(3,0.0);
            vector<blInfoSharedPtr> data = nToNInfo[bit->first];
            NekDouble Dtotal = 0.0;
            for(int i = 0; i < data.size(); i++)
            {
                NekDouble d = bit->first->Distance(data[i]->oNode);
                Dtotal += d;
                sumV[0] += data[i]->N[0] / d;
                sumV[1] += data[i]->N[1] / d;
                sumV[2] += data[i]->N[2] / d;
            }
            sumV[0] *= Dtotal;
            sumV[1] *= Dtotal;
            sumV[2] *= Dtotal;
            NekDouble mag = sqrt(sumV[0]*sumV[0] + sumV[1]*sumV[1] + sumV[2]*sumV[2]);
            sumV[0] /= mag;
            sumV[1] /= mag;
            sumV[2] /= mag;

            Array<OneD, NekDouble> N(3);

            N[0] = (1.0-0.8) * bit->second->N[0] + 0.8 * sumV[0];
            N[1] = (1.0-0.8) * bit->second->N[1] + 0.8 * sumV[1];
            N[2] = (1.0-0.8) * bit->second->N[2] + 0.8 * sumV[2];

            mag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
            N[0] /= mag;
            N[1] /= mag;
            N[2] /= mag;

            bit->second->N = N;
            bit->second->AlignNode(layerT[0]);
        }
    }

    ofstream file;
    file.open("bl.lines");
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        NekDouble l = 0.05;
        file << bit->first->m_x << ", " << bit->first->m_y << ", " << bit->first->m_z << endl;
        file << bit->first->m_x + bit->second->N[0]*l << ", "
             << bit->first->m_y + bit->second->N[1]*l << ", "
             << bit->first->m_z + bit->second->N[2]*l << endl;
        file << endl;
    }
    file.close();


    ASSERTL0(failed == 0, "some normals failed to generate");

    //make prisms
    map<int,int> nm;
    nm[0] = 0;
    nm[1] = 3;
    nm[2] = 4;
    nm[3] = 5;
    nm[4] = 1;
    nm[5] = 2;

    map<ElementSharedPtr,ElementSharedPtr> priToTri;
    map<ElementSharedPtr,ElementSharedPtr> priToPsd;

    ElmtConfig pconf(LibUtilities::ePrism,1,false,false);
    ElmtConfig tconf(LibUtilities::eTriangle,1,false,false);

    vector<boost::tuple<blInfoSharedPtr,blInfoSharedPtr,blInfoSharedPtr> > prisms;

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int>::iterator f = find(m_blsurfs.begin(),
                                                m_blsurfs.end(),
                                                el->CADSurfId);

        if(f == m_blsurfs.end())
        {
            //if this triangle is not in bl surfs continue
            continue;
        }

        vector<NodeSharedPtr> tn(3); //nodes for pseduo surface
        vector<NodeSharedPtr> pn(6); //all prism nodes
        vector<NodeSharedPtr> n = el->GetVertexList();

        for(int j = 0; j < 3; j++)
        {
            pn[nm[j*2+0]] = n[j];
            pn[nm[j*2+1]] = blData[n[j]]->pNode;
            tn[j] = blData[n[j]]->pNode;
        }

        prisms.push_back(make_tuple(blData[n[0]],
                                    blData[n[1]],
                                    blData[n[2]]));

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

        for(int j = 0; j < 3; j++)
        {
            blData[n[j]]->pEls.push_back(T);
        }
    }

    //need to construct a updating tree of the psuedo surface
    //and the reminaing surface mesh minus the symmetry plane
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        m_mesh->m_element[2][i]->SetId(0);
    }
    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        for(int i = 0; i < bit->second->pEls.size(); i++)
        {
            vector<EdgeSharedPtr> es = bit->second->pEls[i]->GetEdgeList();
            for(int j = 0; j < es.size(); j++)
            {
                //if(es[j]->m_n1 == bit->second->pNode ||
                //   es[j]->m_n2 == bit->second->pNode)
                //{
                    bit->second->edges.insert(es[j]);
                //}
            }
        }
    }
    vector<ElementSharedPtr> elsInRtree;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int>::iterator f = find(m_blsurfs.begin(),
                                                m_blsurfs.end(),
                                                el->CADSurfId);

        vector<unsigned int>::iterator s = find(m_symSurfs.begin(),
                                                m_symSurfs.end(),
                                                el->CADSurfId);

        if(f == m_blsurfs.end() && s == m_symSurfs.end())
        {
            elsInRtree.push_back(m_mesh->m_element[2][i]);
        }
    }
    for(int i = 0; i < m_psuedoSurface.size(); i++)
    {
        elsInRtree.push_back(m_psuedoSurface[i]);
    }
    vector<box> boxes;
    for(int i = 0; i < elsInRtree.size(); i++)
    {
        elsInRtree[i]->SetId(i);
        boxes.push_back(GetBox(elsInRtree[i]));
    }

    vector<boxI> inserts;
    for(int i = 0; i < elsInRtree.size(); i++)
    {
        inserts.push_back(make_pair(boxes[i],i));
    }
    bgi::rtree<boxI, bgi::quadratic<16> > rtree(inserts);

    //ofstream file3;
    //file3.open("hit.3D");
    //file3 << "X Y Z value" << endl;
    for(int l = 1; l < m_layer; l++)
    {
        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(bit->second->stopped)
            {
                continue;
            }

            bit->second->bl = l;
            bit->second->AlignNode(layerT[bit->second->bl]);
        }

        for(int i = 0; i < elsInRtree.size(); i++)
        {
            boxes[i] = GetBox(elsInRtree[i]);
        }
        rtree.clear();
        for(int i = 0; i < elsInRtree.size(); i++)
        {
            rtree.insert(make_pair(boxes[i],i));
        }

        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(bit->second->stopped)
            {
                continue;
            }

            vector<blInfoSharedPtr> infos = nToNInfo[bit->first];
            for(int i = 0; i < infos.size(); i++)
            {
                if(bit->second->bl > infos[i]->bl + 1)
                {
                    bit->second->stop = true;
                }
            }
        }

        /*for(int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            if(prisms[i].get<0>()->stopped &&
               prisms[i].get<1>()->stopped &&
               prisms[i].get<2>()->stopped)
            {
                continue;
            }
            ElementSharedPtr el = m_mesh->m_element[3][i];
            SpatialDomains::GeometrySharedPtr geom = el->GetGeom(3);
            SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
            if(!gfac->IsValid())
            {
                if(!prisms[i].get<0>()->stopped)
                {
                    prisms[i].get<0>()->stop = true;
                }
                if(!prisms[i].get<1>()->stopped)
                {
                    prisms[i].get<1>()->stop = true;
                }
                if(!prisms[i].get<2>()->stopped)
                {
                    prisms[i].get<2>()->stop = true;
                }
            }
        }*/

        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(!bit->second->stop && !bit->second->stopped)
            {
                vector<boxI> intersects;
                for(int i = 0; i < bit->second->pEls.size(); i++)
                {
                    rtree.query(bgi::intersects(boxes[bit->second->pEls[i]->GetId()]), back_inserter(intersects));
                }

                set<int> intersectsIds;

                for(int i = 0; i < intersects.size(); i++)
                {
                    intersectsIds.insert(intersects[i].second);
                }

                EdgeSet testEdges;
                set<int>::iterator idit;
                for(idit = intersectsIds.begin(); idit != intersectsIds.end(); idit++)
                {
                    vector<EdgeSharedPtr> es = elsInRtree[(*idit)]->GetEdgeList();
                    for(int j = 0; j < es.size(); j++)
                    {
                        EdgeSet::iterator it = bit->second->edges.find(es[j]);
                        //if(it == bit->second->edges.end())
                        //{
                            testEdges.insert(es[j]);
                        //}
                    }
                }

                bool hit = false;
                for(int i = 0; i < bit->second->pEls.size(); i++)
                {
                    vector<NodeSharedPtr> ns = bit->second->pEls[i]->GetVertexList();
                    NekDouble A0,A1,A2,A3,A4,A5,A6,A7,A8;
                    NekDouble B0,B1,B2;

                    A3 = ns[1]->m_x - ns[0]->m_x;
                    A4 = ns[1]->m_y - ns[0]->m_y;
                    A5 = ns[1]->m_z - ns[0]->m_z;
                    A6 = ns[2]->m_x - ns[0]->m_x;
                    A7 = ns[2]->m_y - ns[0]->m_y;
                    A8 = ns[2]->m_z - ns[0]->m_z;
                    EdgeSet::iterator it;
                    for(it = testEdges.begin(); it != testEdges.end(); it++)
                    {
                        A0 = ((*it)->m_n2->m_x - (*it)->m_n1->m_x) * -1.0;
                        A1 = ((*it)->m_n2->m_y - (*it)->m_n1->m_y) * -1.0;
                        A2 = ((*it)->m_n2->m_z - (*it)->m_n1->m_z) * -1.0;
                        NekDouble det = A0 * (A4*A8 - A7*A5)
                                       -A3 * (A1*A8 - A7*A2)
                                       +A6 * (A1*A5 - A4*A2);

                        if(fabs(det) < 1e-15)
                        {
                            //no intersecton
                            continue;
                        }
                        B0 = (*it)->m_n1->m_x - ns[0]->m_x;
                        B1 = (*it)->m_n1->m_y - ns[0]->m_y;
                        B2 = (*it)->m_n1->m_z - ns[0]->m_z;

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
                        if(X1 > 1e-6 && X2 > 1e-6 && X1 + X2 < 0.999999
                           && X0 > 1e-6 && X0 < 0.999999)
                        {
                            //file3 << bit->second->oNode->m_x << " " << bit->second->oNode->m_y << " " << bit->second->oNode->m_z << " " << 1 << endl;
                            hit = true;
                            bit->second->stop = true;
                            break;
                        }
                    }
                    if(hit) break;
                }
            }
        }

        for(bit = blData.begin(); bit != blData.end(); bit++)
        {
            if(bit->second->stop)
            {
                bit->second->stopped = true;
                bit->second->stop = false;
                bit->second->bl = l-1;
                bit->second->AlignNode(layerT[bit->second->bl]);
            }
        }
    }

    for(bit = blData.begin(); bit != blData.end(); bit++)
    {
        vector<blInfoSharedPtr> infos = nToNInfo[bit->first];
        for(int i = 0; i < infos.size(); i++)
        {
            if(bit->second->bl > infos[i]->bl + 1)
            {
                cout << "non smooth error " << bit->second->bl << " " << infos[i]->bl << endl;
            }
        }
    }

    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[3][i];
        SpatialDomains::GeometrySharedPtr geom = el->GetGeom(3);
        SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
        if(!gfac->IsValid())
        {
            cout << "validity error" << endl;
        }
    }

}

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
    vector<NekDouble> a(N.size());
    while(fabs(dot - 1) > 1e-6)
    {
        ct++;
        Array<OneD, NekDouble> Nplast(3);
        Nplast[0] = Np[0];
        Nplast[1] = Np[1];
        Nplast[2] = Np[2];

        NekDouble aSum = 0.0;
        for(int i = 0; i < N.size(); i++)
        {
            NekDouble dot2 = Np[0]*N[i][0] + Np[1]*N[i][1] + Np[2]*N[i][2];
            if(fabs(dot2 - 1) < 1e-9)
            {
                a[i] = dot2/fabs(dot2) * 1e-9;
            }
            else
            {
                a[i] = acos(dot2);
            }

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
            w[i] = w[i] / wSum;
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

        Np[0] = 0.8* NpN[0] + (1.0-0.8)*Np[0];
        Np[1] = 0.8* NpN[1] + (1.0-0.8)*Np[1];
        Np[2] = 0.8* NpN[2] + (1.0-0.8)*Np[2];
        mag = sqrt(Np[0]*Np[0] + Np[1]*Np[1] + Np[2]*Np[2]);
        Np[0] /= mag;
        Np[1] /= mag;
        Np[2] /= mag;

        dot = Np[0] * Nplast[0] + Np[1] * Nplast[1] + Np[2] * Nplast[2];

        if(ct > 100000)
        {
            cout << "run out of iterations" << endl;
            Np = Ninital;
            break;
        }
    }

    return Np;
}

}
}
