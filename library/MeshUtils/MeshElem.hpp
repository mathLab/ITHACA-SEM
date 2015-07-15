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
            ASSERTL0(search->first == i, "node not on this curve");

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

        int EdgeInCommon(const MeshNodeSharedPtr &n)
        {
            std::vector<int> ne = n->GetEdges();

            for(int i = 0; i < ne.size(); i++)
            {
                for(int j = 0; j < Edges.size(); j++)
                {
                    if(ne[i] == Edges[j])
                    {
                        return ne[i];
                    }
                }
            }

            return -1;
        }

        std::vector<int> SurfsInCommon(const MeshNodeSharedPtr &n)
        {
            std::vector<int> out;

            std::map<int, Array<OneD, NekDouble> > map = n->GetSurfMap();

            std::map<int, Array<OneD, NekDouble> >::iterator it1,it2;

            for(it1 = CADSurf.begin(); it1 != CADSurf.end(); it1++)
            {
                for(it2 = map.begin(); it2 != map.end(); it2++)
                {
                    if(it1->first == it2->first)
                    {
                        out.push_back(it1->first);
                    }
                }
            }

            return out;
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

        int GetId(){return nid;}
        std::map<int, Array<OneD, NekDouble> > GetSurfMap(){return CADSurf;}

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
            oncurve = false;
        }

        void SetCurve(int i)
        {
            curveedge = i;
            oncurve = true;
        }

        int GetCurve()
        {
            if(oncurve)
            {
                return curveedge;
            }
            else
            {
                return -1;
            }
        }

        Array<OneD, MeshNodeSharedPtr> GetN(){return nodes;}

        void SetTri(int i )
        {
            tris.push_back(i);
        }

        void SetHONodes(Array<OneD, MeshNodeSharedPtr> n)
        {
            honodes = n;
        }

    private:

        int eid;
        int curveedge;
        bool oncurve;
        Array<OneD, MeshNodeSharedPtr> honodes;
        Array<OneD, MeshNodeSharedPtr> nodes;
        std::vector<int> tris;

    };


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////


    class MeshTri;
    typedef boost::shared_ptr<MeshTri> MeshTriSharedPtr;

    class MeshTri
    {
    public:
        friend class MemoryManager<MeshTri>;

        MeshTri(int t, MeshNodeSharedPtr a,
                       MeshNodeSharedPtr b,
                       MeshNodeSharedPtr c,
                       MeshEdgeSharedPtr e1,
                       MeshEdgeSharedPtr e2,
                       MeshEdgeSharedPtr e3,
                       int i)
        {
            tid = t;
            nodes = Array<OneD, MeshNodeSharedPtr>(3);
            nodes[0] = a; nodes[1] = b; nodes[2] = c;
            nodes[0]->SetTri(tid);
            nodes[1]->SetTri(tid);
            nodes[2]->SetTri(tid);
            edges = Array<OneD, MeshEdgeSharedPtr>(3);
            edges[0] = e1; edges[1] = e2; edges[2] = e3;
            edges[0]->SetTri(tid);
            edges[1]->SetTri(tid);
            edges[2]->SetTri(tid);
            cid=i;
        }

        void SetNeigh(Array<OneD,int> n)
        {
            neighbours = n;
        }

        int Getcid(){return cid;}
        Array<OneD, MeshNodeSharedPtr> GetN(){return nodes;}

    private:

        int tid;
        int cid;
        Array<OneD, MeshNodeSharedPtr> nodes;
        Array<OneD, MeshEdgeSharedPtr> edges;
        Array<OneD,int> neighbours;
    };

}
}

#endif
