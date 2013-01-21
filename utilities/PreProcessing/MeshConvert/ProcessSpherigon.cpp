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

#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
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
         * Mesh::spherigonFaces to be populated by the relevant input modules.
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
            config["N"] = ConfigOption(false, "8",
                "Number of points to add to face edges.");
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
         * represented as a hybrid triangular/quadrilateral mesh.
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
            
            // First loop over elements and construct vertex normals.
            for (int i = 0; i < el.size(); ++i)
            {
                ElementSharedPtr e = el[i];
                
                // Ensure that element is either a triangle or quad.
                ASSERTL0(e->GetConf().e == eTriangle || 
                         e->GetConf().e == eQuadrilateral,
                         "Spherigon expansions must be either triangles or "
                         "quadrilaterals.");
                
                // Calculate normal for this element.
                NodeSharedPtr node[3] = {
                    e->GetVertex(0), e->GetVertex(1), 
                    e->GetVertex(e->GetConf().e == eQuadrilateral ? 3 : 2)
                };
                
                Node v1 = *(node[1]) - *(node[0]);
                Node v2 = *(node[2]) - *(node[0]);
                Node n;
                UnitCrossProd(v1, v2, n);
                
                // Insert face normal into vertex normal list or add to existing
                // value.
                for (int j = 0; j < e->GetVertexCount(); ++j)
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
            for (nIt = m->vertexNormals.begin(); nIt != m->vertexNormals.end(); ++nIt)
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
            ASSERTL0(m->spaceDim == 3,
                     "Spherigon implementation only valid in 3D.");

            boost::unordered_set<int>::iterator eIt;
            boost::unordered_set<int>           visitedEdges;

            // First construct vector of elements to process.
            vector<ElementSharedPtr> el;

            if (m->verbose)
            {
                cout << "ProcessSpherigon: Smoothing mesh..." << endl;
            }

            if (m->expDim == 2)
            {
                // Manifold case - copy expansion dimension.
                el = m->element[m->expDim];
            }
            else if (m->expDim == 3)
            {
                // Full 3D case - iterate over stored faces and create
                // triangular elements representing faces.
                set<pair<int,int> >::iterator it;
                vector<int> t;
                t.push_back(0);
                
                if (m->spherigonFaces.size() == 0)
                {
                    cerr << "WARNING: Spherigon faces have not been defined -- "
                         << "ignoring smoothing." << endl;
                }
                
                for (it  = m->spherigonFaces.begin();
                     it != m->spherigonFaces.end(); ++it)
                {
                    FaceSharedPtr f = m->element[m->expDim][it->first]->GetFace(
                        it->second);
                    vector<NodeSharedPtr> nodes = f->vertexList;
                    ElementType eType = 
                        nodes.size() == 3 ? eTriangle : eQuadrilateral;
                    ElmtConfig conf(eType, 1, false, false);
                    
                    // Create 2D element.
                    el.push_back(
                        GetElementFactory().CreateInstance(eType,conf,nodes,t));
                }
            }
            else
            {
                ASSERTL0(false, "Spherigon expansions must be 2/3 dimensional");
            }
            
            // See if vertex normals have been generated. If they have not,
            // approximate them by summing normals of surrounding elements.
            if (m->vertexNormals.size() == 0)
            {
                GenerateNormals(el);
            }
            
            int nq = config["N"].as<int>();
            Array<OneD, NekDouble> x(nq*nq);
            Array<OneD, NekDouble> y(nq*nq);
            Array<OneD, NekDouble> z(nq*nq);
            
            ASSERTL0(nq > 2, "Number of points must be greater than 2.");

            int edgeMap[2][4][2] = {
                {{0, 1}, {nq-1, nq}, {nq*(nq-1), -nq}, {-1, -1}},        // tri
                {{0, 1}, {nq-1, nq}, {nq*nq-1, -1},    {nq*(nq-1), -nq}} // quad
            };
            
            for (int i = 0; i < el.size(); ++i)
            {
                // Construct a Nektar++ element to obtain coordinate points
                // inside the element. TODO: Add options for various
                // nodal/tensor point distributions + number of points to add.
                ElementSharedPtr e = el[i];
                
                LibUtilities::BasisKey B0(
                    LibUtilities::eModified_A, nq-1,
                    LibUtilities::PointsKey(
                        nq, LibUtilities::ePolyEvenlySpaced));
                
                if (e->GetConf().e == eTriangle)
                {
                    SpatialDomains::TriGeomSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                            e->GetGeom(3));
                    LocalRegions::TriExpSharedPtr tri =
                        MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(
                            B0, B0, geom);
                    tri->GetCoords(x,y,z);
                }
                else
                {
                    SpatialDomains::QuadGeomSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                            e->GetGeom(3));
                    LocalRegions::QuadExpSharedPtr quad =
                        MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                            B0, B0, geom);
                    quad->GetCoords(x,y,z);
                }
                
                // Find vertex normals.
                int nV = e->GetVertexCount();
                vector<Node> v, vN;
                for (int j = 0; j < nV; ++j)
                {
                    v.push_back(*(e->GetVertex(j)));
                    vN.push_back(m->vertexNormals[v[j].id]);
                }
                
                // Calculate area of element.
                Node ta  = v[1]    - v[0];
                Node tb  = v[nV-1] - v[0];
                double A = CrossProdMag(ta, tb);
                
                vector<Node>   tmp  (nV);
                vector<double> r    (nV);
                vector<Node>   K    (nV);
                vector<Node>   Q    (nV);
                vector<Node>   Qp   (nV);
                vector<double> blend(nV);
                
                // Perform Spherigon method to smooth manifold.
                for (int j = 0; j < nq*nq; ++j)
                {
                    Node P(0, x[j], y[j], z[j]);
                    for (int k = 0; k < nV; ++k)
                    {
                        tmp[k] = P - v[k];
                    }
                    
                    // Calculate generalized barycentric coordinate system (see
                    // equation 6 of paper).
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
                    
                    // Calculate normalised Phong normals (equation 1).
                    Node N(0,0,0,0);
                    for (int k = 0; k < nV; ++k)
                    {
                        r[k] /= weight;
                        N    += vN[k]*r[k];
                    }
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
                    
                    x[j] = P.x;
                    y[j] = P.y;
                    z[j] = P.z;
                }
                
                int offset = (int)e->GetConf().e-2;
                
                // Push new nodes into edges (TODO: face nodes).
                for (int edge = 0; edge < e->GetEdgeCount(); ++edge)
                {
                    eIt = visitedEdges.find(e->GetEdge(edge)->id);
                    if (eIt == visitedEdges.end() || m->expDim == 3)
                    {
                        for (int j = 1; j < nq-1; ++j)
                        {
                            int v = edgeMap[offset][edge][0] + 
                                  j*edgeMap[offset][edge][1];
                            NodeSharedPtr tmp(new Node(0, x[v], y[v], z[v]));
                            e->GetEdge(edge)->edgeNodes.push_back(tmp);
                        }
                        visitedEdges.insert(e->GetEdge(edge)->id);
                    }
                }
            }
            
            // Full 3D only: Copy high-order edge nodes back into original
            // elements.
            if (m->expDim == 3)
            {
                set<pair<int,int> >::iterator       it;
                boost::unordered_set<int>::iterator it2;
                int elCount = 0;
                
                visitedEdges.clear();
                
                for (it  = m->spherigonFaces.begin();
                     it != m->spherigonFaces.end(); ++it, ++elCount)
                {
                    FaceSharedPtr f = m->element[m->expDim][it->first]->GetFace(
                        it->second);
                    for (int edge = 0; edge < f->edgeList.size(); ++edge)
                    {
                        EdgeSharedPtr e = el[elCount]->GetEdge(edge);
                        bool reverseEdge = e->n1 == f->edgeList[edge]->n2;
                        
                        if (visitedEdges.count(f->edgeList[edge]->id) != 0)
                        {
                            continue;
                        }
                        
                        f->edgeList[edge]->edgeNodes = e->edgeNodes;
                        
                        if (reverseEdge)
                        {
                            reverse(f->edgeList[edge]->edgeNodes.begin(),
                                    f->edgeList[edge]->edgeNodes.end());
                        }
                        
                        visitedEdges.insert(f->edgeList[edge]->id);
                    }
                }
            }
        }
    }
}
