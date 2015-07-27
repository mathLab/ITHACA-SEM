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

#include <algorithm>
#include <vector>

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

        bool IsEOC(){return eoc;} //end of curve, speeds up checking
        void SetEOC(){eoc = true;}

        void SetCurve(int i, NekDouble t)
        {
            CADCurve[i] = t;
        }

        void SetSurf(int i, Array<OneD, NekDouble> uv)
        {
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

        bool IsOnSurf(int i)
        {
            std::map<int, Array<OneD, NekDouble> >::iterator s;
            for(s=CADSurf.begin(); s!=CADSurf.end(); s++)
            {
                if(s->first==i)
                {
                    return true;
                }
            }
            return false;
        }

        bool IsOnACurve()
        {
            if(CADCurve.size()>0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        NekDouble Distance(const MeshNodeSharedPtr &n)
        {
            Array<OneD,NekDouble> loc = n->GetLoc();

            return sqrt((m_x-loc[0])*(m_x-loc[0])+
                        (m_y-loc[1])*(m_y-loc[1])+
                        (m_z-loc[2])*(m_z-loc[2]));
        }

        NekDouble Angle(const MeshNodeSharedPtr &a, const MeshNodeSharedPtr &b)
        {
            Array<OneD,NekDouble> la = a->GetLoc();
            Array<OneD,NekDouble> lb = b->GetLoc();
            Array<OneD,NekDouble> va(3),vb(3);
            va[0] = la[0] - m_x;
            va[1] = la[1] - m_y;
            va[2] = la[2] - m_z;
            vb[0] = lb[0] - m_x;
            vb[1] = lb[1] - m_y;
            vb[2] = lb[2] - m_z;

            NekDouble ca = va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
            ca /= sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);
            ca /= sqrt(vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2]);

            if(ca < 0)
            {
                return 3.142 + acos(ca);
            }
            else
            {
                return acos(ca);
            }
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

        void SetID(int i)
        {
            nid = i;
        }

        void SetEdge(int e)
        {
            Edges.push_back(e);
        }

        void RemoveEdge(int e)
        {
            std::vector<int> tmp = Edges;
            int tmps = Edges.size();
            Edges.clear();
            for(int i = 0; i < tmp.size(); i++)
            {
                if(tmp[i] != e)
                {
                    Edges.push_back(tmp[i]);
                }
            }
            ASSERTL0(tmps - Edges.size() == 1,"failed to remove edges");
        }

        void RemoveTri(int t)
        {
            std::vector<int> tmp = Tris;
            int tmps = Tris.size();
            Tris.clear();
            for(int i = 0; i < tmp.size(); i++)
            {
                if(tmp[i] != t)
                {
                    Tris.push_back(tmp[i]);
                }
            }
            ASSERTL0(tmps - Tris.size() == 1,"failed to remove tris");
        }

        void SetTri(int t)
        {
            Tris.push_back(t);
        }

        int GetId(){return nid;}
        std::map<int, Array<OneD, NekDouble> > GetSurfMap(){return CADSurf;}

        std::vector<int> GetEdges(){return Edges;}
        std::vector<int> GetTris(){return Tris;}

        void Move(Array<OneD, NekDouble> l, Array<OneD, NekDouble> uv)
        {
            m_x=l[0];
            m_y=l[1];
            m_z=l[2];
            std::map<int, Array<OneD, NekDouble> >::iterator s =
                        CADSurf.begin();
            CADSurf[s->first]=uv;
        }

        void SetTet(int t)
        {
            Tets.push_back(t);
        }

        std::vector<int> GetTets(){return Tets;}


    private:

        int nid;
        NekDouble m_x, m_y, m_z;
        bool eoc;

        std::map<int, NekDouble> CADCurve;
        std::map<int, Array<OneD, NekDouble> > CADSurf;

        std::vector<int> Edges;
        std::vector<int> Tris;
        std::vector<int> Tets;

    };

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    class MeshEdge;
    typedef boost::shared_ptr<MeshEdge> MeshEdgeSharedPtr;

    class MeshEdge
    {
    public:
        friend class MemoryManager<MeshEdge>;

        MeshEdge(int i, int an, int bn)
        {
            eid=i;
            nodes = Array<OneD, int>(2);
            nodes[0] =an;
            nodes[1] =bn;
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

        Array<OneD, int> GetN(){return nodes;}

        void Swap(int i, int j)
        {
            nodes[0] = i;
            nodes[1] = j;
        }

        void SetTri(int i )
        {
            tris.push_back(i);
        }

        void RemoveTri(int t)
        {
            std::vector<int> tmp = tris;
            int tmps = tris.size();
            tris.clear();
            for(int i = 0; i < tmp.size(); i++)
            {
                if(tmp[i] != t)
                {
                    tris.push_back(tmp[i]);
                }
            }
            ASSERTL0(tmps - tris.size() == 1,"failed to remove tris");
        }

        std::vector<int> GetTri(){return tris;}

        void SetHONodes(std::vector<int> n)
        {
            honodes = n;
        }

        std::vector<int> GetHONodes(int first)
        {
            ASSERTL0(nodes[0] == first || nodes[1] == first,
                        "this node is not in this edge");
            if(nodes[0] != first)
            {
                std::reverse(honodes.begin(),honodes.end());
                int tmp = nodes[0];
                nodes[0] = nodes[1];
                nodes[1] = tmp;
            }

            return honodes;
        }

        int OtherNode(int i)
        {
            ASSERTL0(nodes[0] == i || nodes[1] == i,
                        "this node is not in this edge");
            if(nodes[0]==i)
            {
                return nodes[1];
            }
            else
            {
                return nodes[0];
            }
        }

        void SetSurf(int i){surf=i;}
        int GetSurf(){return surf;}
        int GetId(){return eid;}

    private:

        int eid;
        int curveedge;
        bool oncurve;
        int surf;
        std::vector<int> honodes;
        Array<OneD, int> nodes;
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

        MeshTri(int t, int a,
                       int b,
                       int c,
                       int e1,
                       int e2,
                       int e3,
                       int i)
        {
            tid = t;
            nodes = Array<OneD, int>(3);
            nodes[0] = a; nodes[1] = b; nodes[2] = c;

            edges = Array<OneD, int>(3);
            edges[0] = e1; edges[1] = e2; edges[2] = e3;

            cid=i;
        }

        void SetNeigh(Array<OneD,int> n)
        {
            neighbours = n;
        }

        void Swap(int i, int j, int k)
        {
            nodes[0] = i;
            nodes[1] = j;
            nodes[2] = k;
        }

        void ResetEdges(int i, int j, int k)
        {
            edges[0] = i;
            edges[1] = j;
            edges[2] = k;
        }

        int Getcid(){return cid;}
        Array<OneD, int> GetN(){return nodes;}
        Array<OneD, int> GetE(){return edges;}

    private:

        int tid;
        int cid;
        Array<OneD, int> nodes;
        Array<OneD, int> edges;
        Array<OneD,int> neighbours;
    };

    /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////

    class MeshTet;
    typedef boost::shared_ptr<MeshTet> MeshTetSharedPtr;

    class MeshTet
    {
    public:
        friend class MemoryManager<MeshTet>;

        MeshTet(int t, int a,
                       int b,
                       int c,
                       int d)
        {
            tid = t;
            nodes = Array<OneD, int>(4);
            nodes[0] = a; nodes[1] = b; nodes[2] = c; nodes[3] = d;


        }


        Array<OneD, int> GetN(){return nodes;}

    private:

        int tid;
        Array<OneD, int> nodes;
    };

}
}

#endif
