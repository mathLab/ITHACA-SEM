////////////////////////////////////////////////////////////////////////////////
//
//  File: Octree.h
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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_OCTREE_OCTREE_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_OCTREE_OCTREE_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>
#include <LibUtilities/MeshUtils/CurvaturePoint.h>
#include <LibUtilities/MeshUtils/Octant.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
    namespace LibUtilities {
        namespace MeshUtils {
        
            class Octree
            {
            public:
                friend class MemoryManager<Octree>;
                
                LIB_UTILITIES_EXPORT Octree(const CADSystemSharedPtr &cad) : m_cad(cad)
                {
                };
                
                LIB_UTILITIES_EXPORT void Build(const NekDouble &min,
                                                const NekDouble &max,
                                                const NekDouble &eps);
                
            private:
                
                void CompileCuravturePointList();
                void subdivide(int parent);
                void SmoothOctants();
                void SubDivideLevel(int parent);
                void SmoothSurfaceOctants();
                void PropagateDomain();
                void SmoothAllOctants();
                int CountElemt();
                NekDouble ddx(int i, int j);
                
                NekDouble m_minDelta;
                NekDouble m_maxDelta;
                NekDouble m_eps;
                
                CADSystemSharedPtr m_cad;
                Array<OneD, NekDouble> BoundingBox;
                std::vector<CurvaturePointSharedPtr> m_cpList;
                std::vector<OctantSharedPtr> OctantList;
                int m_totNotDividing;
               
            };
            
            typedef boost::shared_ptr<Octree> OctreeSharedPtr;
        }
    }
}

#endif
