////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_ELM_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_ELM_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace MeshUtils {
    
    class MeshNode;
    typedef boost::shared_ptr<MeshNode> MeshNodeSharedPtr;
    
    class MeshNode
    {
    public:
        friend class MemoryManager<MeshNode>;
        
        MeshNode(int i, NekDouble x, NekDouble y, NekDouble z) :
                   nid(i), m_x(x), m_y(y), m_z(z)
        {
            eoc = false;
        };
        
        bool IsEOC(){return eoc;}
        void SetEOC(){eoc = true;}
        
        void SetCurve(int i, NekDouble t)
        {
            CADCurve[i] = t;
        }
        
        void SetSurf(int i, NekDouble u, NekDouble v)
        {
            Array<OneD, NekDouble> uv(2);
            uv[0]=u;
            uv[1]=v;
            CADSurf[i] = uv;
        }

        Array<OneD, NekDouble> GetLoc()
        {
            Array<OneD, NekDouble> out(3);
            out[0]=m_x;
            out[1]=m_y;
            out[2]=m_z;
            return out;
        }
        
        NekDouble GetC(int i)
        {
            std::map<int, NekDouble>::iterator search =
                            CADCurve.find(i);
            ASSERTL0(search != CADCurve.end(), "node not on this curve");
            
            return search->second;
        }
        
        Array<OneD, NekDouble>  GetS(int i)
        {
            //I dont know why I ahev to do this to get it to work
            //this really needs bound checking
            std::map<int, Array<OneD, NekDouble> >::iterator search =
                        CADSurf.find(i);
            ASSERTL0(search->first == i,"surface not found");
        
            return search->second;
        }
        
        NekDouble Distance(const MeshNodeSharedPtr &n)
        {
            Array<OneD,NekDouble> loc = n->GetLoc();
            
            return sqrt((m_x-loc[0])*(m_x-loc[0])+
                        (m_y-loc[1])*(m_y-loc[1])+
                        (m_z-loc[2])*(m_z-loc[2]));
        }
        
        void SetID(int i)
        {
            nid = i;
        }

        void SetEdge(int e)
        {
            Edges.push_back(e);
        }

        void SetTri(int t)
        {
            Tris.push_back(t);
        }

        std::vector<int> GetEdges(){return Edges;}
        std::vector<int> GetTtris(){return Tris;}
		
        
        
    private:
        
        int nid;
        NekDouble m_x, m_y, m_z;
        bool eoc;
        
        std::map<int, NekDouble> CADCurve;
        std::map<int, Array<OneD, NekDouble> > CADSurf;
        
        std::vector<int> Edges;
        std::vector<int> Tris;
        
    };
    
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    
    class MeshEdge;
    typedef boost::shared_ptr<MeshEdge> MeshEdgeSharedPtr;
    
    class MeshEdge
    {
    public:
        friend class MemoryManager<MeshEdge>;
        
        MeshEdge(int i, MeshNodeSharedPtr an, MeshNodeSharedPtr bn)
        {
            eid=i;
            nodes = Array<OneD, MeshNodeSharedPtr>(2);
            nodes[0] =an;
            nodes[1] =bn;
            nodes[0]->SetEdge(eid);
            nodes[1]->SetEdge(eid);
        }
        
        void SetCurve(int i)
        {
            curveedge = i;
        }
        
    private:
        
        int eid;
        int curveedge;
        Array<OneD, MeshNodeSharedPtr> nodes;
        
    };
    
    
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    
    
    class MeshTri;
    typedef boost::shared_ptr<MeshTri> MeshTriSharedPtr;
    
    class MeshTri
    {
    public:
        friend class MemoryManager<MeshTri>;
        
        MeshTri(int t, MeshNodeSharedPtr an,
                MeshNodeSharedPtr bn, MeshNodeSharedPtr cn)
        {
            tid = t;
            firstn = an;
            secondn = bn;
            thirdn = cn;
            firstn->SetTri(tid);
            secondn->SetTri(tid);
            thirdn->SetTri(tid);
        }
        
        void AddEdge(MeshEdgeSharedPtr e)
        {
            edges.push_back(e);
        }
        
        void SetNeigh(Array<OneD,int> n)
        {
            neighbours = n;
        }
        
        
        std::vector<int> GetReqEdge()
        {
            std::vector<int> f,s,t,out;
            f = firstn->GetEdges();
            s = secondn->GetEdges();
            t = thirdn->GetEdges();
            
            for(int i = 0; i < f.size(); i++)
            {
                for(int j = 0; j < s.size(); i++)
                {
                    if(f[i]==s[j])
                    {
                        out.push_back(f[i]);
                        break;
                    }
                }
            }
            
            for(int i = 0; i < f.size(); i++)
            {
                for(int j = 0; j < t.size(); i++)
                {
                    if(f[i]==t[j])
                    {
                        out.push_back(f[i]);
                        break;
                    }
                }
            }
            
            for(int i = 0; i < t.size(); i++)
            {
                for(int j = 0; j < s.size(); i++)
                {
                    if(t[i]==s[j])
                    {
                        out.push_back(t[i]);
                        break;
                    }
                }
            }
            
            return out;
        }
        
    private:
        
        int tid;
        MeshNodeSharedPtr firstn,secondn,thirdn;
        std::vector<MeshEdgeSharedPtr> edges;
        Array<OneD,int> neighbours;
    };
    
}
}

#endif
