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
#include <SpatialDomains/SegGeom.h>
#include <boost/shared_ptr.hpp>

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

        int Geometry2D::GetFid() const
        {
            return v_GetFid();
        }

        const LibUtilities::BasisSharedPtr Geometry2D::GetEdgeBasis(const int i)
        {
            return v_GetEdgeBasis(i);
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

        const PointGeomSharedPtr Geometry2D::GetVertex(int i) const
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

        void Geometry2D::NewtonIterationForLocCoord(
            const Array<OneD, const NekDouble> &coords,
            const Array<OneD, const NekDouble> &ptsx,
            const Array<OneD, const NekDouble> &ptsy,
                  Array<OneD,       NekDouble> &Lcoords)
        {
            NekDouble xmap,ymap, F1,F2;
            NekDouble der1_x, der2_x, der1_y, der2_y ;
            const Array<TwoD, const NekDouble> &gmat
                = m_geomFactors->GetDerivFactors(GetPointsKeys());
            
            // Unfortunately need the points in an Array to interpolate
            Array<OneD, NekDouble> D1Dx(ptsx.num_elements(),&gmat[0][0]);
            Array<OneD, NekDouble> D1Dy(ptsx.num_elements(),&gmat[1][0]);
            Array<OneD, NekDouble> D2Dx(ptsx.num_elements(),&gmat[2][0]);
            Array<OneD, NekDouble> D2Dy(ptsx.num_elements(),&gmat[3][0]);
            
            int cnt=0; 
            int MaxIterations = 40;
            NekDouble Tol = 1e-12; // This is error*error; 
                
            F1 = F2 = 2000; // Starting value of Function
            
            while(cnt++ < MaxIterations)
            {
                //calculate the global point `corresponding to Lcoords
                xmap = m_xmap->PhysEvaluate(Lcoords, ptsx);
                ymap = m_xmap->PhysEvaluate(Lcoords, ptsy);
                
                F1 = coords[0] - xmap;
                F2 = coords[1] - ymap;

                // stopping criterion
                if(F1*F1 + F2*F2 < Tol)
                {
                    break;
                }
                
                //Interpolate derivative metric at Lcoords
                der1_x = m_xmap->PhysEvaluate(Lcoords, D1Dx);
                der2_x = m_xmap->PhysEvaluate(Lcoords, D1Dy);
                der1_y = m_xmap->PhysEvaluate(Lcoords, D2Dx);
                der2_y = m_xmap->PhysEvaluate(Lcoords, D2Dy);                  
                
                Lcoords[0] = Lcoords[0] + der1_x*(coords[0]-xmap) + 
                    der1_y*(coords[1]-ymap);
                Lcoords[1] = Lcoords[1] + der2_x*(coords[0]-xmap) + 
                    der2_y*(coords[1]-ymap);
            }
            
            if(cnt >= 40)
            {
                Lcoords[0] = Lcoords[1] = 2.0;
            }
        }

        int Geometry2D::v_GetFid() const 
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        const LibUtilities::BasisSharedPtr Geometry2D::v_GetEdgeBasis(const int i)
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

        const PointGeomSharedPtr Geometry2D::v_GetVertex(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            PointGeomSharedPtr returnval;
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

