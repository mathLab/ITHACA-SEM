////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.h
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
//  Description:  This file contains the base class specification for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY_H

#include <SpatialDomains/GeomFactors.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry; // Forward declaration for typedef.
        typedef boost::shared_ptr<Geometry> GeometrySharedPtr;
        typedef std::vector< GeometrySharedPtr > GeometryVector;
        typedef boost::unordered_set< GeometrySharedPtr > GeometrySet;
        typedef boost::shared_ptr <GeometryVector> GeometryVectorSharedPtr;
        typedef std::vector< GeometrySharedPtr >::iterator GeometryVectorIter;

        class PointGeom;
        typedef boost::shared_ptr< PointGeom >  PointGeomSharedPtr;

        struct Curve;
        typedef boost::shared_ptr<Curve> CurveSharedPtr;
        typedef boost::unordered_map<int, CurveSharedPtr> CurveMap;

        /// \brief Less than operator to sort Geometry objects by global id when sorting 
        /// STL containers.
        SPATIAL_DOMAINS_EXPORT  bool SortByGlobalId(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs);

        SPATIAL_DOMAINS_EXPORT  bool GlobalIdEquality(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs);

        /// Base class for shape geometry information
        class Geometry
        {
            public:
                SPATIAL_DOMAINS_EXPORT Geometry();
                SPATIAL_DOMAINS_EXPORT Geometry(int coordim);

                SPATIAL_DOMAINS_EXPORT virtual ~Geometry();

                //---------------------------------------
                // Element connection functions
                //---------------------------------------
                SPATIAL_DOMAINS_EXPORT inline bool IsElmtConnected(
                            int gvo_id,
                            int locid) const;
                SPATIAL_DOMAINS_EXPORT inline void AddElmtConnected(
                            int gvo_id,
                            int locid);
                SPATIAL_DOMAINS_EXPORT inline int  NumElmtConnected() const;

                //---------------------------------------
                // Helper functions
                //---------------------------------------

                SPATIAL_DOMAINS_EXPORT inline int GetCoordim() const;
                SPATIAL_DOMAINS_EXPORT void SetCoordim(int coordim) 
                {
                    m_coordim = coordim;
                }
                SPATIAL_DOMAINS_EXPORT inline GeomFactorsSharedPtr GetGeomFactors();
                SPATIAL_DOMAINS_EXPORT GeomFactorsSharedPtr GetRefGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis);
                SPATIAL_DOMAINS_EXPORT inline GeomFactorsSharedPtr GetMetricInfo();
                SPATIAL_DOMAINS_EXPORT LibUtilities::ShapeType GetShapeType(void);
                SPATIAL_DOMAINS_EXPORT inline int GetGlobalID(void);
                SPATIAL_DOMAINS_EXPORT inline void SetGlobalID(int globalid);
                SPATIAL_DOMAINS_EXPORT inline int GetVid(int i) const;
                SPATIAL_DOMAINS_EXPORT inline int GetEid(int i) const;
                SPATIAL_DOMAINS_EXPORT inline int GetFid(int i) const;
                SPATIAL_DOMAINS_EXPORT inline int GetTid(int i) const;
                SPATIAL_DOMAINS_EXPORT inline int GetNumVerts() const;
                SPATIAL_DOMAINS_EXPORT inline PointGeomSharedPtr GetVertex(int i) const;
                SPATIAL_DOMAINS_EXPORT inline StdRegions::Orientation
                            GetEorient(const int i) const;
                SPATIAL_DOMAINS_EXPORT inline StdRegions::Orientation
                            GetPorient(const int i) const;
                SPATIAL_DOMAINS_EXPORT inline StdRegions::Orientation
                            GetForient(const int i) const;
                SPATIAL_DOMAINS_EXPORT inline int GetNumEdges() const;
                SPATIAL_DOMAINS_EXPORT inline int GetNumFaces() const;
                SPATIAL_DOMAINS_EXPORT inline int GetShapeDim() const;
                SPATIAL_DOMAINS_EXPORT inline StdRegions::StdExpansionSharedPtr
                            GetXmap() const;
                SPATIAL_DOMAINS_EXPORT inline const Array<OneD, const NekDouble> &
                            GetCoeffs(const int i) const;
                SPATIAL_DOMAINS_EXPORT inline bool ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                              NekDouble tol = 0.0);
                SPATIAL_DOMAINS_EXPORT inline bool ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                              Array<OneD, NekDouble> &locCoord,
                              NekDouble tol);
                SPATIAL_DOMAINS_EXPORT inline bool ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                        Array<OneD, NekDouble> &locCoord,
                        NekDouble tol,
                        NekDouble &resid);
                SPATIAL_DOMAINS_EXPORT inline int GetVertexEdgeMap(int i, int j) const;
                SPATIAL_DOMAINS_EXPORT inline int GetVertexFaceMap(int i, int j) const;
                SPATIAL_DOMAINS_EXPORT inline int GetEdgeFaceMap(int i, int j) const;

                SPATIAL_DOMAINS_EXPORT inline void FillGeom();
                SPATIAL_DOMAINS_EXPORT inline NekDouble GetLocCoords(
                        const Array<OneD, const NekDouble> &coords,
                        Array<OneD,       NekDouble> &Lcoords);
                SPATIAL_DOMAINS_EXPORT inline NekDouble GetCoord(
                        const int i, const Array<OneD, const NekDouble> &Lcoord);

                SPATIAL_DOMAINS_EXPORT inline void SetOwnData();
                SPATIAL_DOMAINS_EXPORT inline const LibUtilities::BasisSharedPtr
                            GetBasis(const int i);
                SPATIAL_DOMAINS_EXPORT inline const LibUtilities::PointsKeyVector
                            GetPointsKeys();
                SPATIAL_DOMAINS_EXPORT inline void Reset(
                    CurveMap &curvedEdges,
                    CurveMap &curvedFaces);

            protected:

                SPATIAL_DOMAINS_EXPORT static GeomFactorsSharedPtr
                            ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor);
                static GeomFactorsVector m_regGeomFactorsManager;

                /// coordinate dimension
                int                  m_coordim;
                GeomFactorsSharedPtr m_geomFactors;
                GeomState            m_geomFactorsState;
                StdRegions::StdExpansionSharedPtr m_xmap;

                /// enum identifier to determine if quad points are filled
                GeomState            m_state;
                GeomType             m_geomType;
                LibUtilities::ShapeType   m_shapeType;
                int                  m_globalID;

                Array<OneD, Array<OneD, NekDouble> > m_coeffs;
            
                void GenGeomFactors();

                //---------------------------------------
                // Element connection functions
                //---------------------------------------
                virtual bool v_IsElmtConnected(
                            int gvo_id,
                            int locid) const;
                virtual void v_AddElmtConnected(
                            int gvo_id,
                            int locid);
                virtual int  v_NumElmtConnected() const;

                //---------------------------------------
                // Helper functions
                //---------------------------------------

                virtual int  v_GetEid(int i) const;
                virtual int  v_GetVid(int i) const;
                virtual int  v_GetFid(int i) const;
                virtual void v_GenGeomFactors() = 0;
                virtual int  v_GetNumVerts() const;
                virtual PointGeomSharedPtr v_GetVertex(int i) const = 0;
                virtual StdRegions::Orientation
                             v_GetEorient(const int i) const;
                virtual StdRegions::Orientation
                             v_GetPorient(const int i) const;
                virtual StdRegions::Orientation
                             v_GetForient(const int i) const;
                virtual int  v_GetNumEdges() const;
                virtual int  v_GetNumFaces() const;
                virtual int  v_GetShapeDim() const;
                virtual StdRegions::StdExpansionSharedPtr
                             v_GetXmap() const;
                virtual int  v_GetCoordim() const;
                virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                              NekDouble tol = 0.0);
                virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                        Array<OneD, NekDouble>& locCoord,
                        NekDouble tol);
                virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                        Array<OneD, NekDouble>& locCoord,
                        NekDouble tol,
                        NekDouble &resid);

                virtual int v_GetVertexEdgeMap(int i,int j) const;
                virtual int v_GetVertexFaceMap(int i,int j) const;
                virtual int v_GetEdgeFaceMap(int i,int j) const;

                virtual void v_FillGeom();
                virtual NekDouble v_GetCoord(
                            const int i,
                            const Array<OneD,const NekDouble>& Lcoord);
                virtual NekDouble v_GetLocCoords(
                            const Array<OneD,const NekDouble>& coords,
                            Array<OneD,NekDouble>& Lcoords);

                virtual void v_SetOwnData();
                virtual const LibUtilities::BasisSharedPtr
                             v_GetBasis(const int i);
                virtual void v_Reset(
                    CurveMap &curvedEdges,
                    CurveMap &curvedFaces);
                inline void SetUpCoeffs(const int nCoeffs);
        }; // class Geometry


        struct GeometryHash : std::unary_function<GeometrySharedPtr, std::size_t>
        {
            std::size_t operator()(GeometrySharedPtr const& p) const
            {
                int i;
                size_t seed  = 0;
                int nVert = p->GetNumVerts();
                std::vector<unsigned int> ids(nVert);
                
                for (i = 0; i < nVert; ++i)
                {
                    ids[i] = p->GetVid(i);
                }
                std::sort(ids.begin(), ids.end());
                boost::hash_range(seed, ids.begin(), ids.end());

                return seed;
            }
        };


        inline void Geometry::AddElmtConnected(int gvo_id, int locid)
        {
            return v_AddElmtConnected(gvo_id, locid);
        }

        inline int  Geometry::NumElmtConnected() const
        {
            return v_NumElmtConnected();
        }

        inline bool Geometry::IsElmtConnected(int gvo_id, int locid) const
        {
            return v_IsElmtConnected(gvo_id,locid);
        }

        inline int Geometry::GetCoordim() const
        {
            return v_GetCoordim();
        }

        inline GeomFactorsSharedPtr Geometry::GetGeomFactors()
        {
            GenGeomFactors();
            return ValidateRegGeomFactor(m_geomFactors);
        }

        inline GeomFactorsSharedPtr Geometry::GetMetricInfo()
        {
            return m_geomFactors;
        }

        inline LibUtilities::ShapeType Geometry::GetShapeType()
        {
            return m_shapeType;
        }

        inline int Geometry::GetGlobalID(void)
        {
            return m_globalID;
        }

        inline void Geometry::SetGlobalID(int globalid)
        {
            m_globalID = globalid;
        }

        inline int Geometry::GetVid(int i) const
        {
            return v_GetVid(i);
        }

        inline int Geometry::GetEid(int i) const
        {
            return v_GetEid(i);
        }

        inline int Geometry::GetFid(int i) const
        {
            return v_GetFid(i);
        }

        inline int Geometry::GetTid(int i) const
        {
            const int nDim = GetShapeDim();
            return
                nDim == 1 ? v_GetVid(i) :
                nDim == 2 ? v_GetEid(i) :
                nDim == 3 ? v_GetFid(i) : 0;
        }

        inline int Geometry::GetNumVerts() const
        {
            return v_GetNumVerts();
        }

        inline PointGeomSharedPtr Geometry::GetVertex(int i) const
        {
            return v_GetVertex(i);
        }

        inline StdRegions::Orientation Geometry::GetEorient(const int i) const
        {
            return v_GetEorient(i);
        }

        inline StdRegions::Orientation Geometry::GetPorient(const int i) const
        {
            return v_GetPorient(i);
        }

        inline StdRegions::Orientation Geometry::GetForient(const int i) const
        {
            return v_GetForient(i);
        }

        inline int Geometry::GetNumEdges() const
        {
            return v_GetNumEdges();
        }

        inline int Geometry::GetNumFaces() const
        {
            return v_GetNumFaces();
        }

        inline int Geometry::GetShapeDim() const
        {
            return v_GetShapeDim();
        }

        inline StdRegions::StdExpansionSharedPtr Geometry::GetXmap() const
        {
            return v_GetXmap();
        }

        inline const Array<OneD, const NekDouble> &Geometry::GetCoeffs(const int i) const
        {
            return m_coeffs[i];
        }

        inline bool Geometry::ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                NekDouble tol)
        {
            return v_ContainsPoint(gloCoord,tol);
        }

        inline bool Geometry::ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      Array<OneD, NekDouble> &locCoord,
                      NekDouble tol)
        {
            return v_ContainsPoint(gloCoord,locCoord,tol);
        }

        inline bool Geometry::ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                      Array<OneD, NekDouble> &locCoord,
                NekDouble tol,
                NekDouble &resid)
        {
            return v_ContainsPoint(gloCoord,locCoord,tol,resid);
        }

        inline int Geometry::GetVertexEdgeMap(int i, int j) const
        {
            return v_GetVertexEdgeMap(i,j);
        }

        /// return the id of the \f$j^{th}\f$ face attached to the \f$ i^{th}\f$ vertex
        inline int Geometry::GetVertexFaceMap(int i, int j) const
        {
            return v_GetVertexFaceMap(i,j);
        }

        inline int Geometry::GetEdgeFaceMap(int i, int j) const
        {
            return v_GetEdgeFaceMap(i,j);
        }

        inline void Geometry::GenGeomFactors()
        {
            return v_GenGeomFactors();
        }


       /**
        * @brief Put all quadrature information into face/edge structure and
        * backward transform.
        *
        * @see v_FillGeom()
        */
        inline void Geometry::FillGeom()
        {
            v_FillGeom();
        }

        inline NekDouble Geometry::GetLocCoords(
            const Array<OneD, const NekDouble> &coords,
            Array<OneD,       NekDouble> &Lcoords)
        {
            return v_GetLocCoords(coords, Lcoords);
        }

        /**
         * @brief Given local collapsed coordinate Lcoord return the value of
         * physical coordinate in direction i.
         */
        inline NekDouble Geometry::GetCoord(
            const int i, const Array<OneD, const NekDouble> &Lcoord)
        {
            return v_GetCoord(i, Lcoord);
        }

        inline void Geometry::SetOwnData()
        {
            v_SetOwnData();
        }

        /**
         * @brief Return the j-th basis of the i-th co-ordinate dimension.
         */
        inline const LibUtilities::BasisSharedPtr Geometry::GetBasis(
            const int i)
        {
            return v_GetBasis(i);
        }

        /**
         * @brief Initialise the m_coeffs array.
         */
        inline void Geometry::SetUpCoeffs(const int nCoeffs)
        {
            m_coeffs = Array<OneD, Array<OneD, NekDouble> >(m_coordim);

            for (int i = 0; i < m_coordim; ++i)
            {
                m_coeffs[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
            }
        }

        inline const LibUtilities::PointsKeyVector Geometry::GetPointsKeys()
        {
            return m_xmap->GetPointsKeys();
        }

        inline void Geometry::Reset(CurveMap &curvedEdges,
                                    CurveMap &curvedFaces)
        {
            v_Reset(curvedEdges, curvedFaces);
        }
    }; //end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY_H
