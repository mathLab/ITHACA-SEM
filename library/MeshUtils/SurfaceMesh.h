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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_SURFACEMESH_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_SURFACEMESH_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/CADSystem/CADSurf.h>
#include <MeshUtils/Octree.h>
#include <MeshUtils/CurveMesh.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
namespace MeshUtils {
            
    class SurfaceMesh
    {
    public:
        friend class MemoryManager<SurfaceMesh>;
        
        SurfaceMesh(const LibUtilities::CADSurfSharedPtr &cad,
                    const OctreeSharedPtr &oct,
                    const std::vector<CurveMeshSharedPtr> &cmeshes,
                    const int &order)
                        : m_cadsurf(cad), m_octree(oct),
                          m_curvemeshes(cmeshes),m_order(order)
        
        {
            m_edges = m_cadsurf->GetEdges();
        };
        
        void Mesh();
        
        void Extract(int &nt,
                     int &tnp,
                     Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &p)
        {
            nt = numtris;
            tnp = TotNumPoints;
            p = HOPoints;
        }
        
        
        
    private:
        
        void HOMesh();
        
        void LinearOptimise();
        
        void EdgeSwap();
        
        void Stretching();
        
        bool Validate(int &np,
                      int &nt,
                      Array<OneD, Array<OneD, NekDouble> > &Points,
                      Array<OneD, Array<OneD, int> > &Connec);
        
        void OrientateCurves();
        void AddNewPoint(NekDouble u, NekDouble v);
        
        LibUtilities::CADSurfSharedPtr m_cadsurf;
        OctreeSharedPtr m_octree;
        std::vector<CurveMeshSharedPtr> m_curvemeshes;

        std::vector<std::vector<std::pair<int,int> > > m_edges;
        std::vector<std::vector<std::vector<NekDouble> > > m_uvloops;
        std::vector<std::vector<NekDouble> > m_centers;
        std::vector<std::vector<NekDouble> > m_extrapoints;
        
        Array<OneD, Array<OneD, NekDouble> > Points;
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > HOPoints;
        Array<OneD, Array<OneD, int> > Connec;
        Array<OneD, Array<OneD, int> > neigbourlist;
        int numpoints, numtris;
        int TotNumPoints;
        int m_order;
        
        NekDouble pasr,asr;
        
    };
    
    typedef boost::shared_ptr<SurfaceMesh> SurfaceMeshSharedPtr;
}
}

#endif
