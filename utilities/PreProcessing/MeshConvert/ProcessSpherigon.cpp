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
//  Description: Process a mesh and 'explode' elements from some central point
//  (e.g. for visualisation purposes).
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessSpherigon.h"

#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/PrismExp.h>
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
                ModuleKey(eProcessModule, "spherigon"), ProcessSpherigon::create);
      
        ProcessSpherigon::ProcessSpherigon(MeshSharedPtr m) : ProcessModule(m)
        {

        }
      
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

        double ProcessSpherigon::Blend(double r)
        {
            return r*r;
        }

        void ProcessSpherigon::SuperBlend(
            Node &r, vector<Node> &Q, Node &P, Node &blend)
        {
            double v0_Qi_P_s, v1_Qi_P_s, v2_Qi_P_s;
            int    i;
            
            blend.x = blend.y = blend.z = 0.0;
            
            v0_Qi_P_s = (Q[0]-P).abs2();
            v1_Qi_P_s = (Q[1]-P).abs2();
            v2_Qi_P_s = (Q[2]-P).abs2();
            
            if (r.x > TOL_BLEND && r.x < (1-TOL_BLEND))
            {
                blend.x = r.x*r.x*(
                    r.z*r.z*v2_Qi_P_s/(v2_Qi_P_s+v0_Qi_P_s) +
                    r.y*r.y*v1_Qi_P_s/(v1_Qi_P_s+v0_Qi_P_s));
            }
            if (r.y > TOL_BLEND && r.y < (1-TOL_BLEND))
            {
                blend.y = r.y*r.y*( 
                    r.x*r.x*v0_Qi_P_s/(v0_Qi_P_s+v1_Qi_P_s) + 
                    r.z*r.z*v2_Qi_P_s/(v2_Qi_P_s+v1_Qi_P_s));
            }
            if (r.z > TOL_BLEND && r.z < (1-TOL_BLEND))
            {
                blend.z = r.z*r.z*(
                    r.y*r.y*v1_Qi_P_s/(v1_Qi_P_s+v2_Qi_P_s) + 
                    r.x*r.x*v0_Qi_P_s/(v0_Qi_P_s+v2_Qi_P_s));
            }
            
            blend /= (blend.x + blend.y + blend.z);
            
            if (r.x >= (1-TOL_BLEND))
            {
                blend.x = 1.0;
            }
            if (r.x <= TOL_BLEND)
            {
                blend.x = 0.0;
            }
            if (r.y >= (1-TOL_BLEND))
            {
                blend.y = 1.0;
            }
            if (r.y <= TOL_BLEND)
            {
                blend.y = 0.0;
            }
            if (r.z >= (1-TOL_BLEND))
            {
                blend.z = 1.0;
            }
            if (r.z <= TOL_BLEND)
            {
                blend.z = 0.0;
            }
        }

        void ProcessSpherigon::Process()
        {
            ASSERTL0(m->spaceDim == 3,
                     "Spherigon implementation only valid in 3D.");
            
            // Manifold case.
            if (m->expDim == 2)
            {
                vector<ElementSharedPtr> &el = m->element[m->expDim];
                
                map<int, Node> vertexNormals;
                map<int, Node>::iterator it;
                
                boost::unordered_set<int>  visitedEdges;
                boost::unordered_set<int>::iterator it2;
                
                // First loop over elements and construct vertex normals.
                for (int i = 0; i < el.size(); ++i)
                {
                    ElementSharedPtr e = el[i];
                    
                    // Calculate normal for this element.
                    NodeSharedPtr node[3] = {
                        e->GetVertex(0), e->GetVertex(1), e->GetVertex(2)};
                    
                    Node v1 = *(node[1]) - *(node[0]);
                    Node v2 = *(node[2]) - *(node[0]);
                    Node n(0,0,0,0);
                    UnitCrossProd(v1, v2, n);

                    for (int j = 0; j < 3; ++j)
                    {
                        it = vertexNormals.find(node[j]->id);
                        if (it == vertexNormals.end())
                        {
                            vertexNormals[node[j]->id] = n;
                        }
                        else
                        {
                            it->second += n;
                        }
                    }
                }
                
                // Normalize vertex normals.
                for (it = vertexNormals.begin(); it != vertexNormals.end(); ++it)
                {
                    Node &n = vertexNormals[it->first];
                    double tmp = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
                    n /= tmp;
                }

                for (int i = 0; i < el.size(); ++i)
                {
                    // Construct a Nektar++ element to obtain coordinate
                    // points. TODO: Add options for various nodal/tensor point
                    // distributions.
                    int nq = 8;
                    ElementSharedPtr e = el[i];
                    
                    LibUtilities::BasisKey B0(
                        LibUtilities::eModified_A, nq-1,
                        LibUtilities::PointsKey(
                            nq, LibUtilities::ePolyEvenlySpaced));
                    
                    SpatialDomains::TriGeomSharedPtr geom = 
                        boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                            e->GetGeom(3));
                    
                    LocalRegions::TriExpSharedPtr tri =
                        MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(
                            B0, B0, geom);
                    
                    Array<OneD, NekDouble> x(nq*nq);
                    Array<OneD, NekDouble> y(nq*nq);
                    Array<OneD, NekDouble> z(nq*nq);
                    
                    tri->GetCoords(x,y,z);
                    
                    Node va(*(e->GetVertex(0)));
                    Node vb(*(e->GetVertex(1)));
                    Node vc(*(e->GetVertex(2)));
                    
                    Node ta = vb - va;
                    Node tb = vc - va;
                    
                    Node &vna = vertexNormals[va.id];
                    Node &vnb = vertexNormals[vb.id];
                    Node &vnc = vertexNormals[vc.id];
                    
                    double A = CrossProdMag(ta, tb);
                    
                    // Perform Spherigon method to smooth manifold.
                    for (int j = 0; j < nq*nq; ++j)
                    {
                        Node npos(0, x[j], y[j], z[j]);
                        Node a = npos - va;
                        Node b = npos - vb;
                        Node c = npos - vc;
                        Node r(0, CrossProdMag(b,c) / A, CrossProdMag(c,a) / A,
                               CrossProdMag(a,b) / A);
                        
                        // Calculate normalised Phong normals
                        Node N(0,vna.x*r.x + vnb.x*r.y + vnc.x*r.z,
                                 vna.y*r.x + vnb.y*r.y + vnc.y*r.z,
                                 vna.z*r.x + vnb.z*r.y + vnc.z*r.z);
                        N /= sqrt(N.abs2());
                        
                        double conta = (va.x-x[j])*N.x + (va.y-y[j])*N.y + (va.z-z[j])*N.z;
                        double contb = (vb.x-x[j])*N.x + (vb.y-y[j])*N.y + (vb.z-z[j])*N.z;
                        double contc = (vc.x-x[j])*N.x + (vc.y-y[j])*N.y + (vc.z-z[j])*N.z;
                        
                        //double ka[3], kb[3], kc[3];
                        
                        Node ka = npos + N*conta;
                        Node kb = npos + N*contb;
                        Node kc = npos + N*contc;
                        
                        conta = (va.x-ka.x)*vna.x + (va.y-ka.y)*vna.y +
                            (va.z-ka.z)*vna.z;
                        contb = (vb.x-kb.x)*vnb.x + (vb.y-kb.y)*vnb.y +
                            (vb.z-kb.z)*vnb.z;
                        contc = (vc.x-kc.x)*vnc.x + (vc.y-kc.y)*vnc.y +
                            (vc.z-kc.z)*vnc.z;
                        
                        double ma = conta/(1 + N.dot(vna));
                        double mb = contb/(1 + N.dot(vnb));
                        double mc = contc/(1 + N.dot(vnc));
                        
                        Node qa = ka + N*ma;
                        Node qb = kb + N*mb;
                        Node qc = kc + N*mc;
                        
                        // Uncomment this region and comment out region below
                        // for C0 spherigon.
                        /*
                        double ba = Blend(r.x);
                        double bb = Blend(r.y);
                        double bc = Blend(r.z);
                        double invBlend = 1.0/(ba+bb+bc);

                        x[j] = (ba*qa.x + bb*qb.x + bc*qc.x)*invBlend;
                        y[j] = (ba*qa.y + bb*qb.y + bc*qc.y)*invBlend;
                        z[j] = (ba*qa.z + bb*qb.z + bc*qc.z)*invBlend;
                        */
                        
                        conta = (va.x-x[j])*N.x + (va.y-y[j])*N.y + (va.z-z[j])*N.z;
                        contb = (vb.x-x[j])*N.x + (vb.y-y[j])*N.y + (vb.z-z[j])*N.z;
                        contc = (vc.x-x[j])*N.x + (vc.y-y[j])*N.y + (vc.z-z[j])*N.z;

                        vector<Node> Q(3);
                        Q[0] = va - N*conta;
                        Q[1] = vb - N*contb;
                        Q[2] = vc - N*contc;
                        
                        Node blend;
                        SuperBlend(r, Q, npos, blend);
                        
                        x[j] = blend.x*qa.x + blend.y*qb.x + blend.z*qc.x;
                        y[j] = blend.x*qa.y + blend.y*qb.y + blend.z*qc.y;
                        z[j] = blend.x*qa.z + blend.y*qb.z + blend.z*qc.z;
                    }
                    
                    // Push new nodes into edges (TODO: face nodes).
                    int triEdgeMap[3][2] = {
                        {0, 1}, {nq-1, nq}, {nq*(nq-1), -nq}
                    };
                    
                    for (int edge = 0; edge < 3; ++edge)
                    {
                        it2 = visitedEdges.find(e->GetEdge(edge)->id);
                        if (it2 == visitedEdges.end())
                        {
                            for (int j = 1; j < nq-1; ++j)
                            {
                                int v = triEdgeMap[edge][0] + j*triEdgeMap[edge][1];
                                NodeSharedPtr tmp(new Node(0, x[v], y[v], z[v]));
                                e->GetEdge(edge)->edgeNodes.push_back(tmp);
                            }
                            visitedEdges.insert(e->GetEdge(edge)->id);
                        }
                    }
                }
            }
        }
    }
}
