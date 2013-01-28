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
#include <SpatialDomains/GeometryShapeType.h>

#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {

        const char* const GeomShapeTypeMap[] =
        {
            "NoGeomShapeType",
            "Segment",
            "Point",
            "Triangle",
            "Quadrilateral",
            "Tetrahedron",
            "Pyramid",
            "Prism",
            "Hexahedron"
        };

        class Geometry; // Forward declaration for typedef.
        typedef boost::shared_ptr<Geometry> GeometrySharedPtr;
        typedef std::vector< GeometrySharedPtr > GeometryVector;
        typedef boost::unordered_set< GeometrySharedPtr > GeometrySet;
        typedef boost::shared_ptr <GeometryVector> GeometryVectorSharedPtr;
        typedef std::vector< GeometrySharedPtr >::iterator GeometryVectorIter;



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
                SPATIAL_DOMAINS_EXPORT bool IsElmtConnected(
                            int gvo_id,
                            int locid) const;
                SPATIAL_DOMAINS_EXPORT void AddElmtConnected(
                            int gvo_id,
                            int locid);
                SPATIAL_DOMAINS_EXPORT int  NumElmtConnected() const;

                //---------------------------------------
                // Helper functions
                //---------------------------------------

                SPATIAL_DOMAINS_EXPORT void FillGeom();
                SPATIAL_DOMAINS_EXPORT void GetLocCoords(
                        const Array<OneD, const NekDouble> &coords,
                              Array<OneD,       NekDouble> &Lcoords);
                SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(
                        const int i, const Array<OneD, const NekDouble> &Lcoord);

                SPATIAL_DOMAINS_EXPORT void SetOwnData();
                SPATIAL_DOMAINS_EXPORT Array<OneD,NekDouble>&
                            UpdatePhys(const int i);
                SPATIAL_DOMAINS_EXPORT const LibUtilities::BasisSharedPtr 
                            GetBasis(const int i, const int j);

                SPATIAL_DOMAINS_EXPORT GeomType GetGtype();
                SPATIAL_DOMAINS_EXPORT const Array<OneD, const NekDouble>& GetJac();
                SPATIAL_DOMAINS_EXPORT const Array<TwoD, const NekDouble>& GetGmat();
                SPATIAL_DOMAINS_EXPORT const int GetCoordim() const;
                SPATIAL_DOMAINS_EXPORT GeomFactorsSharedPtr GetGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis);
                SPATIAL_DOMAINS_EXPORT GeomFactorsSharedPtr GetMetricInfo();
                SPATIAL_DOMAINS_EXPORT GeomShapeType GetGeomShapeType(void);
                SPATIAL_DOMAINS_EXPORT int GetGlobalID(void);
                SPATIAL_DOMAINS_EXPORT void SetGlobalID(int globalid);
                SPATIAL_DOMAINS_EXPORT int GetVid(int i) const;
                SPATIAL_DOMAINS_EXPORT int GetEid(int i = 0) const;
                SPATIAL_DOMAINS_EXPORT int GetNumVerts() const;
                SPATIAL_DOMAINS_EXPORT StdRegions::Orientation
                            GetEorient(const int i) const;
                SPATIAL_DOMAINS_EXPORT StdRegions::Orientation
                            GetPorient(const int i) const;
                SPATIAL_DOMAINS_EXPORT int GetNumEdges() const;
                SPATIAL_DOMAINS_EXPORT int GetNumFaces() const;
                SPATIAL_DOMAINS_EXPORT int GetShapeDim() const;
                SPATIAL_DOMAINS_EXPORT bool ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                              NekDouble tol = 0.0);
		SPATIAL_DOMAINS_EXPORT int GetVertexEdgeMap(int i, int j) const;
		SPATIAL_DOMAINS_EXPORT int GetVertexFaceMap(int i, int j) const;
		SPATIAL_DOMAINS_EXPORT int GetEdgeFaceMap(int i, int j) const;


            protected:

                SPATIAL_DOMAINS_EXPORT static GeomFactorsSharedPtr
                            ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor);

                /// coordinate dimension
                int                  m_coordim;

                GeomFactorsSharedPtr m_geomFactors;
                GeomState            m_geomFactorsState;

                /// enum identifier to determine if quad points are filled
                GeomState            m_state;

                static GeomFactorsVector m_regGeomFactorsManager;

                GeomShapeType m_geomShapeType;
                int           m_globalID;

                void GenGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis);

        private:
                GeomType m_geomType;


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
                virtual void v_GenGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis);
                virtual int  v_GetVid(int i) const;
                virtual int  v_GetNumVerts() const;
                virtual StdRegions::Orientation
                             v_GetEorient(const int i) const;
                virtual StdRegions::Orientation
                             v_GetPorient(const int i) const;
                virtual int  v_GetNumEdges() const;
                virtual int  v_GetNumFaces() const;
                virtual int  v_GetShapeDim() const;
                virtual int  v_GetCoordim() const;
                virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                              NekDouble tol = 0.0);
		virtual int v_GetVertexEdgeMap(int i,int j) const;
		virtual int v_GetVertexFaceMap(int i,int j) const;
		virtual int v_GetEdgeFaceMap(int i,int j) const;
                virtual void v_FillGeom();
                virtual NekDouble v_GetCoord(
                            const int i,
                            const Array<OneD,const NekDouble>& Lcoord);
                virtual void v_GetLocCoords(
                            const Array<OneD,const NekDouble>& coords,
                                  Array<OneD,NekDouble>& Lcoords);

                virtual void v_SetOwnData();
                virtual Array<OneD,NekDouble>& v_UpdatePhys(const int i);
                virtual const LibUtilities::BasisSharedPtr
                             v_GetBasis(const int i, const int j);

        }; // class Geometry


        struct GeometryHash : std::unary_function<GeometrySharedPtr, std::size_t>
        {
            std::size_t operator()(GeometrySharedPtr const& p) const
            {
                int i;
                size_t seed  = 0;
                int nVert = p->GetNumVerts();
                int nEdge = p->GetNumEdges();
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

    }; //end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY_H
