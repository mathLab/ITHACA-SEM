////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry2D.cpp
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
//  Description:  2D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#include <SpatialDomains/Geometry2D.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        Geometry2D::Geometry2D()
        {
        }

        Geometry2D::Geometry2D(const int coordim):
            Geometry(coordim)
        {
            ASSERTL0(m_coordim > 1,
                     "Coordinate dimension should be at least 2 for a 2D geometry");
        }

        Geometry2D::~Geometry2D()
        {
        }


        StdRegions::StdExpansion2DSharedPtr Geometry2D::GetXmap(const int i)
        {
            return m_xmap[i];
        }

        int Geometry2D::GetFid() const
        {
            return v_GetFid();
        }



        const LibUtilities::BasisSharedPtr Geometry2D::GetEdgeBasis(const int i, const int j)
        {
            return v_GetEdgeBasis(i,j);
        }


        const Geometry2DSharedPtr Geometry2D::GetFace(int i) const
        {
            return v_GetFace(i);
        }

        StdRegions::Orientation Geometry2D::GetFaceOrient(const int i) const
        {
            return v_GetFaceOrient(i);
        }

        const Geometry1DSharedPtr Geometry2D::GetEdge(int i) const
        {
            return v_GetEdge(i);
        }

        const VertexComponentSharedPtr Geometry2D::GetVertex(int i) const
        {
            return v_GetVertex(i);
        }

        StdRegions::Orientation Geometry2D::GetCartesianEorient(const int i) const
        {
            return v_GetCartesianEorient(i);
        }

        int Geometry2D::WhichEdge(SegGeomSharedPtr edge)
        {
            return v_WhichEdge(edge);
        }

        int Geometry2D::WhichFace(Geometry2DSharedPtr face)
        {
            return v_WhichFace(face);
        }

        StdRegions::StdExpansion2DSharedPtr Geometry2D::operator[](const int i) const
        {
            if((i>=0)&& (i<m_coordim))
            {
                return m_xmap[i];
            }

            NEKERROR(ErrorUtil::efatal,
                     "Invalid Index used in [] operator");
            return m_xmap[0]; //should never be reached
        }

        int Geometry2D::v_GetFid() const 
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        const LibUtilities::BasisSharedPtr Geometry2D::v_GetEdgeBasis(const int i, const int j)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            LibUtilities::BasisSharedPtr returnval;
            return returnval;
        }


        int Geometry2D::v_GetEid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        const Geometry1DSharedPtr Geometry2D::v_GetEdge(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            SegGeomSharedPtr returnval;
            return returnval;
        }

        const VertexComponentSharedPtr Geometry2D::v_GetVertex(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            VertexComponentSharedPtr returnval;
            return returnval;
        }

        StdRegions::Orientation Geometry2D::v_GetEorient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return StdRegions::eForwards;
        }

        const Geometry2DSharedPtr Geometry2D::v_GetFace(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            Geometry2DSharedPtr returnval;
            return returnval;
        }

        StdRegions::Orientation Geometry2D::v_GetFaceOrient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
        }

        StdRegions::Orientation Geometry2D::v_GetCartesianEorient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return StdRegions::eForwards;
        }

        int Geometry2D::v_WhichEdge(SegGeomSharedPtr edge)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry2D::v_WhichFace(Geometry2DSharedPtr face)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry2D::v_GetShapeDim() const
        {
            return 2;
        }

        bool Geometry2D::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      NekDouble tol)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return false;
        }


  }; //end of namespace
}; //end of namespace

