////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.cpp
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
//  Description:  This file contains the base class implementation for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace SpatialDomains
    {

        // static class property
        GeomFactorsVector Geometry::m_regGeomFactorsManager;

        Geometry::Geometry():
            m_coordim(0),
            m_geomFactorsState(eNotFilled),
            m_state(eNotFilled),
            m_shapeType(LibUtilities::eNoShapeType),
            m_globalID(-1)
        {
        }

        Geometry::Geometry(const int coordim):
            m_coordim(coordim),
            m_geomFactorsState(eNotFilled),
            m_state(eNotFilled),
            m_shapeType(LibUtilities::eNoShapeType),
            m_globalID(-1)
        {
        }

        Geometry::~Geometry()
        {
        }

        GeomFactorsSharedPtr Geometry::ValidateRegGeomFactor(
                GeomFactorsSharedPtr geomFactor)
        {
            GeomFactorsSharedPtr returnval = geomFactor;

/// \todo should this '#if 0' statement be removed?
#if 0
            bool found = false;
            if (geomFactor->GetGtype() == eRegular)
            {
                for (GeomFactorsVectorIter iter = m_regGeomFactorsManager.begin();
                    iter != m_regGeomFactorsManager.end();
                    ++iter)
                {
                    if (**iter == *geomFactor)
                    {
                        returnval = *iter;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    m_regGeomFactorsManager.push_back(geomFactor);
                    returnval = geomFactor;
                }
            }
#endif
            return returnval;
        }

        bool SortByGlobalId(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs)
        {
            return lhs->GetGlobalID() < rhs->GetGlobalID();
        }

        bool GlobalIdEquality(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs)
        {
            return lhs->GetGlobalID() == rhs->GetGlobalID();
        }

        void Geometry::v_AddElmtConnected(int gvo_id, int locid)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
        }

        int Geometry::v_NumElmtConnected() const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        bool Geometry::v_IsElmtConnected(int gvo_id, int locid) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return false;
        }

        int Geometry::v_GetVid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry::v_GetEid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry::v_GetFid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return 0;
        }

        int Geometry::v_GetNumVerts() const
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is only valid for shape type geometries");
            return 0;
        }

        StdRegions::Orientation Geometry::v_GetEorient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is not valid for this geometry.");
            return StdRegions::eForwards;
        }

        StdRegions::Orientation Geometry::v_GetPorient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is not valid for this geometry.");
                    return StdRegions::eFwd;
        }

        StdRegions::Orientation Geometry::v_GetForient(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is not valid for this geometry.");
                    return StdRegions::eFwd;
        }

        int Geometry::v_GetNumEdges() const
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry::v_GetNumFaces() const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        int Geometry::v_GetShapeDim() const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        StdRegions::StdExpansionSharedPtr Geometry::v_GetXmap() const
        {
            return m_xmap;
        }

        bool Geometry::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      NekDouble tol)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return false;
        }

        bool Geometry::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                Array<OneD, NekDouble>& locCoord,
                NekDouble tol)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return false;
        }

        bool Geometry::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                Array<OneD, NekDouble>& locCoord,
                NekDouble tol,
                NekDouble &resid)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return false;
        }

        int Geometry::v_GetVertexEdgeMap(const int i, const int j) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return 0;
        }

        int Geometry::v_GetVertexFaceMap(const int i, const int j) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return 0;
        }

        int Geometry::v_GetEdgeFaceMap(const int i, const int j) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function has not been defined for this geometry");
            return 0;
        }

        NekDouble Geometry::v_GetCoord(
                    const int i,
                    const Array<OneD,const NekDouble>& Lcoord)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return 0.0;
        }

        NekDouble Geometry::v_GetLocCoords(
                const Array<OneD,const NekDouble> &coords,
                      Array<OneD,NekDouble> &Lcoords)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return 0.0;
        }

        void Geometry::v_FillGeom()
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
        }

        void Geometry::v_SetOwnData()
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
        }

        const LibUtilities::BasisSharedPtr Geometry::v_GetBasis(const int i)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            LibUtilities::BasisSharedPtr returnval;
            return returnval;
        }

        int Geometry::v_GetCoordim() const
        {
            return m_coordim;
        }

        /**
         * @brief Reset this geometry object: unset the current state and remove
         * allocated GeomFactors.
         */
        void Geometry::v_Reset(CurveMap &curvedEdges,
                               CurveMap &curvedFaces)
        {
            // Reset state
            m_state            = eNotFilled;
            m_geomFactorsState = eNotFilled;

            // Junk geometric factors
            m_geomFactors = GeomFactorsSharedPtr();
        }
    }; //end of namespace
}; //end of namespace

