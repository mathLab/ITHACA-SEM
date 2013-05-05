////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSpherigon.cpp
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
//  Description: Apply Spherigon surface smoothing technique to a 3D mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessSpherigon.h"

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#define TOL_BLEND 1.0e-8

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessSpherigon::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "spherigon"),
                ProcessSpherigon::create);
        
        /**
         * @class ProcessSpherigon
         * 
         * This class implements the spherigon surface smoothing technique which
         * is documented in
         * 
         *   "The SPHERIGON: A Simple Polygon Patch for Smoothing Quickly your
         *   Polygonal Meshes": P. Volino and N. Magnenat Thalmann, Computer
         *   Animation Proceedings (1998).
         * 
         * This implementation works in both a 2D manifold setting (for
         * triangles and quadrilaterals embedded in 3-space) and in a full 3D
         * enviroment (prisms, tetrahedra and hexahedra). 
         * 
         * No additional information needs to be supplied directly to the module
         * in order for it to work. However, 3D elements do rely on the mapping
         * Mesh::spherigonSurfs to be populated by the relevant input modules.
         * 
         * Additionally, since the algorithm assumes normals are supplied which
         * are perpendicular to the true surface defined at each vertex. If
         * these are specified in Mesh::vertexNormals by the input module,
         * better smoothing results can be obtained. Otherwise, normals are
         * estimated by taking the average of all edge/face normals which
         * connect to the vertex.
         */
      
        /**
         * @brief Default constructor.
         */
        ProcessSpherigon::ProcessSpherigon(MeshSharedPtr m) : ProcessModule(m)
        {
            config["N"] = ConfigOption(false, "5",
                "Number of points to add to face edges.");
            config["surf"] = ConfigOption(false, "-1",
                "Tag identifying surface to process.");
            config["BothTriFacesOnPrism"] = ConfigOption(true, "-1",
                "Curve both triangular faces of prism on boundary.");
        }
      
        /**
         * @brief Destructor.
         */
        ProcessSpherigon::~ProcessSpherigon()
        {
            
        }

        /*
         * @brief Calcuate the normalised cross product \f$ \vec{c} =
         * \vec{a}\times\vec{b} / \| \vec{a}\times\vec{b} \| \f$.
         */
        void ProcessSpherigon::UnitCrossProd(Node &a, Node &b, Node &c)
        {
            double inv_mag;

            c.x = a.y*b.z - a.z*b.y;
            c.y = a.z*b.x - a.x*b.z;
            c.z = a.x*b.y - a.y*b.x;

            inv_mag = 1.0/sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
            
            c.x = c.x*inv_mag;
            c.y = c.y*inv_mag;
            c.z = c.z*inv_mag;
        }
        
        /**
         * @brief Calculate the magnitude of the cross product \f$
         * \vec{a}\times\vec{b} \f$.
         */
        double ProcessSpherigon::CrossProdMag(Node &a, Node &b)
        {
            double tmp1 = a.y*b.z - a.z*b.y;
            double tmp2 = a.z*b.x - a.x*b.z;
            double tmp3 = a.x*b.y - a.y*b.x;
            return sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3);
        }

        /**
         * @brief Calculate the \f$ C^0 \f$ blending function for spherigon
         * implementation.
         * 
         * See equation (5) of the paper.
         * 
         * @param r  Barycentric coordinate.
         */
        double ProcessSpherigon::Blend(double r)
        {
            return r*r;
        }

        /**
         * @brief Calculate the \f$ C^1 \f$ blending function for spherigon
         * implementation.
         * 
         * See equation (10) of the paper.
         * 
         * @param r      Generalised barycentric coordinates of the point P.
         * @param Q      Vector of vertices denoting this triangle/quad.
         * @param P      Point in the triangle to apply blending to.
         * @param blend  The resulting blending components for each vertex.
         */
        void ProcessSpherigon::SuperBlend(
            vector<double> &r, vector<Node> &Q, Node &P, vector<double> &blend)
        {
            vector<double> tmp(r.size());
            double         totBlend = 0.0;
            int            i;
            int            nV = r.size();

            for (i = 0; i < nV; ++i)
            {
                blend[i] = 0.0;
                tmp  [i] = (Q[i]-P).abs2();
            }

            for (i = 0; i < nV; ++i)
            {
                int ip = (i+1) % nV, im = (i-1+nV) % nV;
                
                if (r[i] > TOL_BLEND && r[i] < (1-TOL_BLEND))
                {
                    blend[i] = r[i]*r[i]*(
                        r[im]*r[im]*tmp[im]/(tmp[im] + tmp[i]) +
                        r[ip]*r[ip]*tmp[ip]/(tmp[ip] + tmp[i]));
                    totBlend += blend[i];
                }
            }
            
            for (i = 0; i < nV; ++i)
            {
                blend[i] /= totBlend;
                if (r[i] >= (1-TOL_BLEND))
                {
                    blend[i] = 1.0;
                }
                if (r[i] <= TOL_BLEND)
                {
                    blend[i] = 0.0;
                }
            }
        }

        /**
         * @brief Generate a set of approximate vertex normals to a surface
         * represented by line segments in 2D and a hybrid
         * triangular/quadrilateral mesh in 3D.
         * 
         * This routine approximates the true vertex normals to a surface by
         * averaging the normals of all edges/faces which connect to the
         * vertex. It is better to use the exact surface normals which can be
         * set in Mesh::vertexNormals, but where they are not supplied this
         * routine calculates an approximation for the spherigon implementation.
         * 
         * @param el  Vector of elements denoting the surface mesh.
         */
        void ProcessSpherigon::GenerateNormals(
            std::vector<ElementSharedPtr> &el)
        {
            boost::unordered_map<int, Node>::iterator nIt;
            
            for (int i = 0; i < el.size(); ++i)
            {
                ElementSharedPtr e = el[i];
                
                // Ensure that element is a line, triangle or quad.
                ASSERTL0(e->GetConf().e == eLine         ||
                         e->GetConf().e == eTriangle     || 
                         e->GetConf().e == eQuadrilateral,
                         "Spherigon expansions must be lines, triangles or "
                         "quadrilaterals.");
                
                // Calculate normal for this element.
                int nV = e->GetVertexCount();
                vector<NodeSharedPtr> node(nV);

                for (int j = 0; j < nV; ++j)
                {
                    node[j] = e->GetVertex(j);
                }
                
                Node n;
                
                if (m->spaceDim == 3)
                {
                    // Create two tangent vectors and take unit cross product.
                    Node v1 = *(node[1]) - *(node[0]);
                    Node v2 = *(node[2]) - *(node[0]);
                    UnitCrossProd(v1, v2, n);
                }
                else
                {
                    // Calculate gradient vector and invert.
                    Node dx  = *(node[1]) - *(node[0]);
                    dx      /= sqrt(dx.abs2());
                    n.x      = -dx.y;
                    n.y      = dx.x;
                    n.z      = 0;
                }
                
                // Insert face normal into vertex normal list or add to existing
                // value.
                for (int j = 0; j < nV; ++j)
                {
                    nIt = m->vertexNormals.find(e->GetVertex(j)->id);
                    if (nIt == m->vertexNormals.end())
                    {
                        m->vertexNormals[e->GetVertex(j)->id] = n;
                    }
                    else
                    {
                        nIt->second += n;
                    }
                }
            }
            
            // Normalize resulting vectors.
            for (nIt  = m->vertexNormals.begin();
                 nIt != m->vertexNormals.end  (); ++nIt)
            {
                Node &n = m->vertexNormals[nIt->first];
                n /= sqrt(n.abs2());
            }
        }

        /**
         * @brief Perform the spherigon smoothing technique on the mesh.
         */
        void ProcessSpherigon::Process()
        {
            ASSERTL0(m->spaceDim == 3 || m->spaceDim == 2,
                     "Spherigon implementation only valid in 2D/3D.");

            boost::unordered_set<int>::iterator eIt;
            boost::unordered_set<int>           visitedEdges;

            // First construct vector of elements to process.
            vector<ElementSharedPtr> el;

            if (m->verbose)
            {
                cout << "ProcessSpherigon: Smoothing mesh..." << endl;
            }

            if (m->expDim == 2 && m->spaceDim == 3)
            {
                // Manifold case - copy expansion dimension.
                el = m->element[m->expDim];
            }
            else if (m->spaceDim == m->expDim)
            {
                // Full 2D or 3D case - iterate over stored edges/faces and
                // create segments/triangles representing those edges/faces.
                set<pair<int,int> >::iterator it;
                vector<int> t;
                t.push_back(0);
                
                // Construct list of spherigon edges/faces from a tag.
                int surfTag = config["surf"].as<int>();
                bool prismTag = config["BothTriFacesOnPrism"].beenSet;
                if (surfTag != -1)
                {
                    m->spherigonSurfs.clear();
                    for (int i = 0; i < m->element[m->expDim].size(); ++i)
                    {
                        ElementSharedPtr el = m->element[m->expDim][i];
                        int nSurf = m->expDim == 3 ? el->GetFaceCount() : 
                                                     el->GetEdgeCount();
                        
                        for (int j = 0; j < nSurf; ++j)
                        {
                            int bl = el->GetBoundaryLink(j);
                            if (bl == -1)
                            {
                                continue;
                            }

                            ElementSharedPtr bEl  = m->element[m->expDim-1][bl];
                            vector<int>      tags = bEl->GetTagList();

                            if (find(tags.begin(), tags.end(), surfTag) !=
                                tags.end())
                            {
                                m->spherigonSurfs.insert(make_pair(i, j));
                                
                                // Curve other tri face on Prism. Note
                                // could be problem on pyramid when
                                // implemented
                                if((nSurf == 5)&&prismTag)
                                {
                                    // add other end of prism on boundary for smoothing
                                    int triFace = (j == 1)? 3:1;
                                    
                                    m->spherigonSurfs.insert(make_pair(i, triFace));
                                }
                            }

                        }
                    }
                }
                
                if (m->spherigonSurfs.size() == 0)
                {
                    cerr << "WARNING: Spherigon surfaces have not been defined "
                         << "-- ignoring smoothing." << endl;
                    return;
                }

                if (m->expDim == 3)
                {
                    for (it  = m->spherigonSurfs.begin();
                         it != m->spherigonSurfs.end  (); ++it)
                    {
                        FaceSharedPtr f = m->element[m->expDim][it->first]->
                            GetFace(it->second);
                        vector<NodeSharedPtr> nodes = f->vertexList;
                        ElementType eType = (ElementType)(nodes.size()-1);
                        ElmtConfig conf(eType, 1, false, false);
                        
                        // Create 2D element.
                        ElementSharedPtr elmt = GetElementFactory().
                            CreateInstance(eType,conf,nodes,t);
                        
                        // Copy vertices/edges from face.
                        for (int i = 0; i < f->vertexList.size(); ++i)
                        {
                            elmt->SetVertex(i, f->vertexList[i]);
                        }
                        for (int i = 0; i < f->edgeList.size(); ++i)
                        {
                            elmt->SetEdge(i, f->edgeList[i]);
                        }
                        
                        el.push_back(elmt);
                    }
                }
                else
                {
                    for (it  = m->spherigonSurfs.begin();
                         it != m->spherigonSurfs.end  (); ++it)
                    {
                        EdgeSharedPtr edge = m->element[m->expDim][it->first]->
                            GetEdge(it->second);
                        vector<NodeSharedPtr> nodes;
                        ElementType eType = eLine;
                        ElmtConfig conf(eType, 1, false, false);
                        
                        nodes.push_back(edge->n1);
                        nodes.push_back(edge->n2);
                        
                        // Create 2D element.
                        ElementSharedPtr elmt = GetElementFactory().
                            CreateInstance(eType,conf,nodes,t);
                        
                        // Copy vertices/edges from original element.
                        elmt->SetVertex(0, nodes[0]);
                        elmt->SetVertex(1, nodes[1]);
                        elmt->SetEdge(
                            0, m->element[m->expDim][it->first]->
                                GetEdge(it->second));
                        el.push_back(elmt);
                    }
                }
            }
            else
            {
                ASSERTL0(false, "Spherigon expansions must be 2/3 dimensional");
            }
            
            // See if vertex normals have been generated. If they have not,
            // approximate them by summing normals of surrounding elements.
            bool normalsGenerated = false;
            if (m->vertexNormals.size() == 0)
            {
                GenerateNormals(el);
                normalsGenerated = true;
            }

            // Allocate storage for interior points.
            int nq = config["N"].as<int>();
            int nquad = m->spaceDim == 3 ? nq*nq : nq;
            Array<OneD, NekDouble> x(nq*nq);
            Array<OneD, NekDouble> y(nq*nq);
            Array<OneD, NekDouble> z(nq*nq);
            
            ASSERTL0(nq > 2, "Number of points must be greater than 2.");

            LibUtilities::BasisKey B0(
                LibUtilities::eOrtho_A, nq,
                LibUtilities::PointsKey(
                    nq, LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B1(
                LibUtilities::eOrtho_B, nq,
                LibUtilities::PointsKey(
                    nq, LibUtilities::eGaussRadauMAlpha1Beta0));
            StdRegions::StdNodalTriExpSharedPtr stdtri =
                MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
                    B0, B1, LibUtilities::eNodalTriElec);

            Array<OneD, NekDouble> xnodal(nq*(nq+1)/2), ynodal(nq*(nq+1)/2);
            stdtri->GetNodalPoints(xnodal, ynodal);
            
            int edgeMap[3][4][2] = {
                {{0, 1}, {-1,   -1}, {-1,        -1 }, {-1,        -1}}, // seg
                {{0, 1}, {nq-1, nq}, {nq*(nq-1), -nq}, {-1,        -1}}, // tri
                {{0, 1}, {nq-1, nq}, {nq*nq-1,   -1 }, {nq*(nq-1), -nq}} // quad
            };
            
            int vertMap[3][4][2] = {
                {{0, 1}, {0, 0}, {0, 0}, {0, 0}}, // seg
                {{0, 1}, {1, 2}, {2, 3}, {0, 0}}, // tri
                {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, // quad
            };
            
            for (int i = 0; i < el.size(); ++i)
            {
                // Construct a Nektar++ element to obtain coordinate points
                // inside the element. TODO: Add options for various
                // nodal/tensor point distributions + number of points to add.
                ElementSharedPtr e = el[i];

                LibUtilities::BasisKey B2(
                    LibUtilities::eModified_A, nq,
                    LibUtilities::PointsKey(
                        nq, LibUtilities::eGaussLobattoLegendre));
                
                if (e->GetConf().e == eLine)
                {
                    SpatialDomains::SegGeomSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                            e->GetGeom(m->spaceDim));
                    LocalRegions::SegExpSharedPtr seg =
                        MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(
                            B2, geom);
                    seg->GetCoords(x,y,z);
                    nquad = nq;
                }
                else if (e->GetConf().e == eTriangle)
                {
                    SpatialDomains::TriGeomSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                            e->GetGeom(3));
                    LocalRegions::NodalTriExpSharedPtr tri =
                        MemoryManager<LocalRegions::NodalTriExp>
                            ::AllocateSharedPtr(
                                B0, B1, LibUtilities::eNodalTriElec, geom);

                    Array<OneD, NekDouble> coord(2);
                    tri->GetCoords(x,y,z);
                    nquad = nq*(nq+1)/2;
                    Vmath::Vcopy(nq*nq, x, 1, stdtri->UpdatePhys(), 1);

                    for (int j = 0; j < nquad; ++j)
                    {
                        coord[0] = xnodal[j];
                        coord[1] = ynodal[j];
                        x[j] = stdtri->PhysEvaluate(coord);
                    }
                    
                    Vmath::Vcopy(nq*nq, y, 1, stdtri->UpdatePhys(), 1);
                    for (int j = 0; j < nquad; ++j)
                    {
                        coord[0] = xnodal[j];
                        coord[1] = ynodal[j];
                        y[j] = stdtri->PhysEvaluate(coord);
                    }

                    Vmath::Vcopy(nq*nq, z, 1, stdtri->UpdatePhys(), 1);
                    for (int j = 0; j < nquad; ++j)
                    {
                        coord[0] = xnodal[j];
                        coord[1] = ynodal[j];
                        z[j] = stdtri->PhysEvaluate(coord);
                    }
                }
                else if (e->GetConf().e == eQuadrilateral)
                {
                    SpatialDomains::QuadGeomSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                            e->GetGeom(3));
                    LocalRegions::QuadExpSharedPtr quad =
                        MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                            B2, B2, geom);
                    quad->GetCoords(x,y,z);
                    nquad = nq*nq;
                }
                else
                {
                    ASSERTL0(false, "Unknown expansion type.");
                }
                
                // Zero z-coordinate in 2D.
                if (m->spaceDim == 2)
                {
                    Vmath::Zero(nquad, z, 1);
                }
                
                // Find vertex normals.
                int nV = e->GetVertexCount();
                vector<Node> v, vN;
                for (int j = 0; j < nV; ++j)
                {
                    v.push_back(*(e->GetVertex(j)));
                    vN.push_back(m->vertexNormals[v[j].id]);
                }

                vector<Node>   tmp  (nV);
                vector<double> r    (nV);
                vector<Node>   K    (nV);
                vector<Node>   Q    (nV);
                vector<Node>   Qp   (nV);
                vector<double> blend(nV);
                vector<Node>   out  (nquad);

                // Calculate segment length for 2D spherigon routine.
                double segLength = sqrt((v[0] - v[1]).abs2());
                
                // Perform Spherigon method to smooth manifold.
                for (int j = 0; j < nquad; ++j)
                {
                    Node P(0, x[j], y[j], z[j]);
                    Node N(0,0,0,0);

                    // Calculate generalised barycentric coordinates r[] and the
                    // Phong normal N = vN . r for this point of the element.
                    if (m->spaceDim == 2)
                    {
                        // In 2D the coordinates are given by a ratio of the
                        // segment length to the distance from one of the
                        // endpoints.
                        r[0] = sqrt((P - v[0]).abs2()) / segLength;
                        r[0] = max(min(1.0, r[0]), 0.0);
                        r[1] = 1.0 - r[0];
                        
                        // Calculate Phong normal.
                        N = vN[0]*r[0] + vN[1]*r[1];
                    }
                    else if (m->spaceDim == 3)
                    {
                        for (int k = 0; k < nV; ++k)
                        {
                            tmp[k] = P - v[k];
                        }
                        
                        // Calculate generalized barycentric coordinate system
                        // (see equation 6 of paper).
                        double weight = 0.0;
                        for (int k = 0; k < nV; ++k)
                        {
                            r[k] = 1.0;
                            for (int l = 0; l < nV-2; ++l)
                            {
                                r[k] *= CrossProdMag(tmp[(k+l+1) % nV], 
                                                     tmp[(k+l+2) % nV]);
                            }
                            weight += r[k];
                        }
                        
                        // Calculate Phong normal (equation 1).
                        for (int k = 0; k < nV; ++k)
                        {
                            r[k] /= weight;
                            N    += vN[k]*r[k];
                        }
                    }
                    
                    // Normalise Phong normal.
                    N /= sqrt(N.abs2());
                    
                    for (int k = 0; k < nV; ++k)
                    {
                        // Perform steps denoted in equations 2, 3, 8 for C1
                        // smoothing.
                        double tmp1;
                        K[k]  = P+N*((v[k]-P).dot(N));
                        tmp1  = (v[k]-K[k]).dot(vN[k]) / (1.0 + N.dot(vN[k]));
                        Q[k]  = K[k] + N*tmp1;
                        Qp[k] = v[k] - N*((v[k]-P).dot(N));
                    }
                    
                    // Apply C1 blending function to the surface. TODO: Add
                    // option to do (more efficient) C0 blending function.
                    SuperBlend(r, Qp, P, blend);
                    P.x = P.y = P.z = 0.0;
                    
                    // Apply blending (equation 4).
                    for (int k = 0; k < nV; ++k)
                    {
                        P += Q[k]*blend[k];
                    }
                    
                    out[j] = P;
                }
                
                // Push nodes into lines - TODO: face interior nodes. 
                // offset = 0 (seg), 1 (tri) or 2 (quad)
                int offset = (int)e->GetConf().e-1;
                
                for (int edge = 0; edge < e->GetEdgeCount(); ++edge)
                {
                    eIt = visitedEdges.find(e->GetEdge(edge)->id);
                    if (eIt == visitedEdges.end())
                    {
                        bool reverseEdge = !(v[vertMap[offset][edge][0]] ==
                                             *(e->GetEdge(edge)->n1));
                        
                        if (e->GetConf().e == eQuadrilateral)
                        {
                            for (int j = 1; j < nq-1; ++j)
                            {
                                int v = edgeMap[offset][edge][0] + 
                                    j*edgeMap[offset][edge][1];
                                e->GetEdge(edge)->edgeNodes.push_back(
                                    NodeSharedPtr(new Node(out[v])));
                            }
                        }
                        else
                        {
                            for (int j = 0; j < nq-2; ++j)
                            {
                                int v = 3 + edge*(nq-2) + j;
                                e->GetEdge(edge)->edgeNodes.push_back(
                                    NodeSharedPtr(new Node(out[v])));
                            }
                        }
                        
                        if (reverseEdge)
                        {
                            reverse(e->GetEdge(edge)->edgeNodes.begin(),
                                    e->GetEdge(edge)->edgeNodes.end());
                        }

                        e->GetEdge(edge)->curveType =
                            LibUtilities::eGaussLobattoLegendre;

                        visitedEdges.insert(e->GetEdge(edge)->id);
                    }
                }
                
                // Add face nodes in manifold and full 3D case, but not for 2D.
                if (m->spaceDim == 3)
                {
                    vector<NodeSharedPtr> volNodes;
                    
                    if (e->GetConf().e == eQuadrilateral)
                    {
                        volNodes.resize((nq-2)*(nq-2));
                        for (int j = 1; j < nq-1; ++j)
                        {
                            for (int k = 1; k < nq-1; ++k)
                            {
                                int v = j*nq+k;
                                volNodes[(j-1)*(nq-2)+(k-1)] =
                                    NodeSharedPtr(new Node(out[v]));
                            }
                        }
                    }
                    else
                    {
                        for (int j = 3+3*(nq-2); j < nquad; ++j)
                        {
                            volNodes.push_back(NodeSharedPtr(new Node(out[j])));
                        }
                    }
                    
                    e->SetVolumeNodes(volNodes);
                }
            }

            // Copy face nodes back into 3D element faces.
            if (m->expDim == 3)
            {
                set<pair<int,int> >::iterator it;
                int elmt = 0;
                for (it  = m->spherigonSurfs.begin();
                     it != m->spherigonSurfs.end  (); ++it, ++elmt)
                {
                    FaceSharedPtr f = m->element[m->expDim][it->first]->
                        GetFace(it->second);
                    f->faceNodes = el[elmt]->GetVolumeNodes();
                    f->curveType = LibUtilities::eNodalTriElec;
                }
            }

            if (normalsGenerated)
            {
                m->vertexNormals.clear();
            }
        }
    }
}
