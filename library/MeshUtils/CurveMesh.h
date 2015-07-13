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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_CURVEMESH_CURVEMESH_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_CURVEMESH_CURVEMESH_H

#include <boost/shared_ptr.hpp>

#include <MeshUtils/MeshElem.hpp>
#include <LibUtilities/CADSystem/CADCurve.h>
#include <MeshUtils/Octree.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
namespace MeshUtils {
            
    class CurveMesh
    {
    public:
        friend class MemoryManager<CurveMesh>;
        
        CurveMesh(bool verbose,
                  int id,
                  const LibUtilities::CADCurveSharedPtr &cad,
                  const OctreeSharedPtr &oct) :
                  m_cadcurve(cad),m_octree(oct),m_verbose(verbose),
                  m_id(id)
        {
        };
        
        void Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                  std::map<int, MeshEdgeSharedPtr> &Edges);
        
        int GetFirstPoint(){return m_meshpoints[0];}
        int GetLastPoint(){return m_meshpoints[Ne];}
        std::vector<int> GetMeshPoints()
                    {return m_meshpoints;}
        int GetNumPoints(){return Ne+1;}
        
    private:
        
        void GetSampleFunction();
        void GetPhiFunction();
        NekDouble EvaluateDS(NekDouble s);
        NekDouble EvaluatePS(NekDouble s);
        
        LibUtilities::CADCurveSharedPtr m_cadcurve;
        OctreeSharedPtr m_octree;
        NekDouble m_curvelength;
        int m_numSamplePoints;
        Array<OneD, NekDouble> m_bounds;
        std::vector<std::vector<NekDouble> > m_dst;
        std::vector<std::vector<NekDouble> > m_ps;
        NekDouble Ae;
        NekDouble ds;
        int Ne;
        std::vector<NekDouble> meshsvalue;
        std::vector<int> m_meshpoints;
        
        bool m_verbose;
        
        int m_id;
        
    };
    
    typedef boost::shared_ptr<CurveMesh> CurveMeshSharedPtr;
}
}

#endif
