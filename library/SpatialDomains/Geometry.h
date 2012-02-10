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

#include "pchSpatialDomains.h"

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/GeomFactors.h>

#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {

        class Geometry1D;

        // Types of geometry types.
        enum GeomShapeType
        {
            eNoGeomShapeType,
            eSegment,
            ePoint,
            eTriangle,
            eQuadrilateral,
            eTetrahedron,
            ePyramid,
            ePrism,
            eHexahedron,
            SIZE_GeomShapeType
        };


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

        class Geometry
        {
            public:
                SPATIAL_DOMAINS_EXPORT Geometry();
                SPATIAL_DOMAINS_EXPORT Geometry(int coordim);

                SPATIAL_DOMAINS_EXPORT virtual ~Geometry();


                inline GeomType GetGtype()
                {
                    return m_geomFactors->GetGtype();
                }

                inline const Array<OneD, const NekDouble> &GetJac()
                {
                    return m_geomFactors->GetJac();
                }

                inline const Array<TwoD, const NekDouble>& GetGmat()
                {
                    return m_geomFactors->GetGmat();
                }

                inline const int GetCoordim() const
                {
                    return m_coordim;
                }

                inline GeomFactorsSharedPtr GetGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
                {
                    GenGeomFactors(tbasis);
                    return ValidateRegGeomFactor(m_geomFactors);
                }

                inline GeomFactorsSharedPtr GetMetricInfo()
                {
                    return m_geomFactors;
                }

                inline GeomShapeType GetGeomShapeType(void)
                {
                    return m_geomShapeType;
                }

                inline int GetGlobalID(void)
                {
                    return m_globalID;
                }

                void SetGlobalID(int globalid)
                {
                    m_globalID = globalid;
                }

                // Wrappers around virtual Functions
                inline int GetVid(int i) const
                {
                    return v_GetVid(i);
                }

                inline int GetNumVerts() const
                {
                    return v_GetNumVerts();
                }

                inline StdRegions::EdgeOrientation GetEorient(const int i) const
                {
                    return v_GetEorient(i);
                }
			
				inline StdRegions::PointOrientation GetPorient(const int i) const
				{
					return v_GetPorient(i);
				}

                inline int GetNumEdges() const
                {
                    return v_GetNumEdges();
                }

                inline int GetShapeDim() const
                {
                    return v_GetShapeDim();
                }

                inline bool ContainsPoint(
                                          const Array<OneD, const NekDouble> &gloCoord, NekDouble tol = 0.0)
                {
                    return v_ContainsPoint(gloCoord,tol);
                }

            protected:

                SPATIAL_DOMAINS_EXPORT static GeomFactorsSharedPtr ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor);

                int                  m_coordim;     // coordinate dimension
                GeomFactorsSharedPtr m_geomFactors;
                GeomState            m_state;       // enum identifier to determine if quad points are filled
                static GeomFactorsVector m_regGeomFactorsManager;

                GeomShapeType m_geomShapeType;
                int           m_globalID;

                void GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
                {
                    return v_GenGeomFactors(tbasis);
                }

        private:
                GeomType m_geomType;

                virtual int v_GetEid(int i) const
                {
                    NEKERROR(ErrorUtil::efatal,
                             "This function is only valid for shape type geometries");
                    return 0;
                }

                virtual void v_GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
                {
                    NEKERROR(ErrorUtil::efatal,
                        "This function is only valid for shape type geometries");
                }

                virtual int v_GetVid(int i) const
                {
                    NEKERROR(ErrorUtil::efatal,
                             "This function is only valid for expansion type geometries");
                    return 0;
                }


                virtual int v_GetNumVerts() const
                {
                    NEKERROR(ErrorUtil::efatal,
                        "This function is only valid for shape type geometries");
                    return 0;
                }

                virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
                {
                    NEKERROR(ErrorUtil::efatal,
                        "This function is not valid for this geometry.");
                    return StdRegions::eForwards;
                }
			
				virtual StdRegions::PointOrientation v_GetPorient(const int i) const
				{
					NEKERROR(ErrorUtil::efatal,
							 "This function is not valid for this geometry.");
					return StdRegions::eFwd;
				}

                virtual int v_GetNumEdges() const
                {
                    NEKERROR(ErrorUtil::efatal,
                        "This function is only valid for shape type geometries");
                    return 0;
                }


                virtual int v_GetShapeDim() const
                {
                    NEKERROR(ErrorUtil::efatal,
                             "This function is only valid for shape type geometries");
                    return 0;
                }

                virtual bool v_ContainsPoint(
                                             const Array<OneD, const NekDouble> &gloCoord, NekDouble tol = 0.0)
                {
                    NEKERROR(ErrorUtil::efatal,
                             "This function has not been defined for this geometry");
                    return false;
                }
        };

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

        /// \brief Less than operator to sort Geometry objects by global id when sorting 
        /// STL containers.
       SPATIAL_DOMAINS_EXPORT  bool SortByGlobalId(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs);

       SPATIAL_DOMAINS_EXPORT  bool GlobalIdEquality(const boost::shared_ptr<Geometry>& lhs, 
            const boost::shared_ptr<Geometry>& rhs);

    }; //end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY_H
