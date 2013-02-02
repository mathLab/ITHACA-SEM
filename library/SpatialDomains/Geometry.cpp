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
            m_state(eNotFilled),
            m_geomFactorsState(eNotFilled),
            m_geomShapeType(eNoGeomShapeType),
            m_globalID(-1)
        {
        }

        Geometry::Geometry(const int coordim):
            m_coordim(coordim),
            m_state(eNotFilled),
            m_geomFactorsState(eNotFilled),
            m_geomShapeType(eNoGeomShapeType),
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


        void Geometry::AddElmtConnected(int gvo_id, int locid)
        {
            return v_AddElmtConnected(gvo_id, locid);
        }

        int  Geometry::NumElmtConnected() const
        {
            return v_NumElmtConnected();
        }

        bool Geometry::IsElmtConnected(int gvo_id, int locid) const
        {
            return v_IsElmtConnected(gvo_id,locid);
        }


        GeomType Geometry::GetGtype()
        {
            return m_geomFactors->GetGtype();
        }

        const Array<OneD, const NekDouble>& Geometry::GetJac()
        {
            return m_geomFactors->GetJac();
        }

        const Array<TwoD, const NekDouble>& Geometry::GetGmat()
        {
            return m_geomFactors->GetGmat();
        }

        const int Geometry::GetCoordim() const
        {
            return v_GetCoordim();
        }

        GeomFactorsSharedPtr Geometry::GetGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis)
        {
            GenGeomFactors(tbasis);
            return ValidateRegGeomFactor(m_geomFactors);
        }

        GeomFactorsSharedPtr Geometry::GetMetricInfo()
        {
            return m_geomFactors;
        }

        GeomShapeType Geometry::GetGeomShapeType(void)
        {
            return m_geomShapeType;
        }

        int Geometry::GetGlobalID(void)
        {
            return m_globalID;
        }

        void Geometry::SetGlobalID(int globalid)
        {
            m_globalID = globalid;
        }

        int Geometry::GetVid(int i) const
        {
            return v_GetVid(i);
        }

        int Geometry::GetEid(int i) const
        {
            return v_GetEid(i);
        }


        int Geometry::GetNumVerts() const
        {
            return v_GetNumVerts();
        }

        StdRegions::Orientation Geometry::GetEorient(const int i) const
        {
            return v_GetEorient(i);
        }

        StdRegions::Orientation Geometry::GetPorient(const int i) const
        {
            return v_GetPorient(i);
        }

        int Geometry::GetNumEdges() const
        {
            return v_GetNumEdges();
        }

        int Geometry::GetNumFaces() const
        {
            return v_GetNumFaces();
        }

        int Geometry::GetShapeDim() const
        {
            return v_GetShapeDim();
        }

        bool Geometry::ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      NekDouble tol)
        {
            return v_ContainsPoint(gloCoord,tol);
        }

        int Geometry::GetVertexEdgeMap(int i, int j) const
	{
	    return v_GetVertexEdgeMap(i,j);
        }

        int Geometry::GetVertexFaceMap(int i, int j) const
        {
            return v_GetVertexFaceMap(i,j);
        }

        int Geometry::GetEdgeFaceMap(int i, int j) const
        {
            return v_GetEdgeFaceMap(i,j);
        }

        void Geometry::GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis)
        {
            return v_GenGeomFactors(tbasis);
        }


       /** 
        * @brief Put all quadrature information into face/edge structure and
        * backward transform.
        * 
        * @see v_FillGeom()
        */
        void Geometry::FillGeom()
        {
            v_FillGeom();
        }

        void Geometry::GetLocCoords(
            const Array<OneD, const NekDouble> &coords,
                  Array<OneD,       NekDouble> &Lcoords)
        {
            v_GetLocCoords(coords, Lcoords);
        }

        /** 
         * @brief Given local collapsed coordinate Lcoord return the value of
         * physical coordinate in direction i.
         */
        NekDouble Geometry::GetCoord(
            const int i, const Array<OneD, const NekDouble> &Lcoord)
        {
            return v_GetCoord(i, Lcoord);
        }

        void Geometry::SetOwnData()
        {
            v_SetOwnData();
        }

        /**
        * @brief Return a reference to the physical space of co-ordinate
        * dimension i.
        */
        Array<OneD,NekDouble>& Geometry::UpdatePhys(const int i)
        {
            return v_UpdatePhys(i);
        }

        /**
         * @brief Return the j-th basis of the i-th co-ordinate dimension.
         */
        const LibUtilities::BasisSharedPtr Geometry::GetBasis(
            const int i, const int j)
        {
            return v_GetBasis(i, j);
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

        int Geometry::v_GetEid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for shape type geometries");
            return 0;
        }

        void Geometry::v_GenGeomFactors(
                    const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis)
        {
            NEKERROR(ErrorUtil::efatal,
                "This function is only valid for shape type geometries");
        }

        int Geometry::v_GetVid(int i) const
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

        bool Geometry::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      NekDouble tol)
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

        void Geometry::v_GetLocCoords(
                const Array<OneD,const NekDouble> &coords,
                      Array<OneD,NekDouble> &Lcoords)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
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

        Array<OneD,NekDouble>& Geometry::v_UpdatePhys(const int i)
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return NullNekDouble1DArray;
        }

        const LibUtilities::BasisSharedPtr Geometry::v_GetBasis(
                    const int i,
                    const int j)
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

    }; //end of namespace
}; //end of namespace

