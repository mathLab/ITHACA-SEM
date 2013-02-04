////////////////////////////////////////////////////////////////////////////////
//
//  File:  MeshGraph2D.h
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH2D_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH2D_H

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>

#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        typedef boost::shared_ptr<QuadGeom> SharedQuadGeomPtr;
        typedef std::vector< SharedQuadGeomPtr >      QuadGeomVector;

        class MeshGraph2D:
            public MeshGraph
        {

        public:

            SPATIAL_DOMAINS_EXPORT MeshGraph2D();
            SPATIAL_DOMAINS_EXPORT MeshGraph2D(const LibUtilities::SessionReaderSharedPtr &pSession);
            SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph2D();

            SPATIAL_DOMAINS_EXPORT void ReadGeometry(const std::string &infilename);
            SPATIAL_DOMAINS_EXPORT void ReadGeometry(TiXmlDocument &doc);

            SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GetSegGeom(int eID);

            inline int GetCoordim(void){
                return GetSpaceDimension();
            }

            inline const TriGeomMap &GetTrigeoms(void) const
            {
                return m_triGeoms;
            }

            inline const QuadGeomMap &GetQuadgeoms(void) const
            {
                return m_quadGeoms;
            }

#ifdef OLD
            inline bool GetGeofac_defined(void)
            {
                return m_geoFacDefined;
            }
#endif

            void GenXGeoFac();

            inline int GetNseggeoms() const
            {
                return int(m_segGeoms.size());
            }

            inline int GetVidFromElmt(StdRegions::ExpansionType expansion,
                const int vert, const int elmt) const
            {
                if(expansion == StdRegions::eTriangle)
                {
                    ASSERTL2(m_triGeoms.find(elmt) != m_triGeoms.end(),
                        "eid is out of range");

                    return m_triGeoms.find(elmt)->second->GetVid(vert);
                }
                else
                {
                    ASSERTL2(m_quadGeoms.find(elmt) != m_quadGeoms.end(),
                        "eid is out of range");

                    return m_quadGeoms.find(elmt)->second->GetVid(vert);
                }
            }

            inline int GetEidFromElmt(StdRegions::ExpansionType expansion,
                const int edge, const int elmt) const
            {
                if(expansion == StdRegions::eTriangle)
                {
                    ASSERTL2(m_triGeoms.find(elmt) != m_triGeoms.end(),
                        "eid is out of range");

                    return m_triGeoms.find(elmt)->second->GetEid(edge);
                }
                else
                {
                    ASSERTL2(m_quadGeoms.find(elmt) != m_quadGeoms.end(),
                        "eid is out of range");

                    return m_quadGeoms.find(elmt)->second->GetEid(edge);
                }
            }

            inline StdRegions::Orientation GetEorientFromElmt(StdRegions::ExpansionType expansion,const int edge, const int elmt) const
            {
                if(expansion == StdRegions::eTriangle)
                {
                    ASSERTL2(m_triGeoms.find(elmt) != m_triGeoms.end(),
                        "eid is out of range");

                    return m_triGeoms.find(elmt)->second->GetEorient(edge);
                }
                else
                {
                    ASSERTL2(m_quadGeoms.find(elmt) != m_quadGeoms.end(),
                        "eid is out of range");

                    return m_quadGeoms.find(elmt)->second->GetEorient(edge);
                }
            }


            inline StdRegions::Orientation GetCartesianEorientFromElmt(StdRegions::ExpansionType expansion,const int edge, const int elmt) const
            {
                StdRegions::Orientation returnval;

                if(expansion == StdRegions::eTriangle)
                {
                    ASSERTL2(m_triGeoms.find(elmt) != m_triGeoms.end(),
                        "eid is out of range");

                    returnval = m_triGeoms.find(elmt)->second->GetEorient(edge);
                }
                else
                {
                    ASSERTL2(m_quadGeoms.find(elmt) != m_quadGeoms.end(),
                        "eid is out of range");

                    returnval =  m_quadGeoms.find(elmt)->second->GetEorient(edge);
                }

                // swap orientation if on edge 2 & 3 (if quad)
                if(edge >= 2)
                {
                    if(returnval == StdRegions::eForwards)
                    {
                        returnval = StdRegions::eBackwards;
                    }
                    else
                    {
                        returnval = StdRegions::eForwards;
                    }
                }
                return returnval;
            }

            int GetNumComposites(void)
            {
                return int(m_meshComposites.size());
            }

            int GetNumCompositeItems(int whichComposite)
            {
                int returnval = -1;

                try
                {
                    returnval = int(m_meshComposites[whichComposite]->size());
                }
                catch(...)
                {
                    std::ostringstream errStream;
                    errStream << "Unable to access composite item [" << whichComposite << "].";
                    NEKERROR(ErrorUtil::efatal, errStream.str());
                }

                return returnval;
            }

            /// \brief Return the elements (shared ptrs) that have this edge.

            //ElementEdgeVectorSharedPtr GetElementsFromEdge(SegGeomSharedPtr edge);
            SPATIAL_DOMAINS_EXPORT ElementEdgeVectorSharedPtr GetElementsFromEdge(SegGeomSharedPtr edge);

            SPATIAL_DOMAINS_EXPORT ElementEdgeVectorSharedPtr GetElementsFromEdge(Geometry1DSharedPtr edge);

            /** \brief Return the BasisKey corresponding to an edge of an element
	     *  If the expansion is a triangle the Modified_B direction is modified to be one-dimensional Modified_A,GaussLobattoLegendre.
	     **/
            SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKey GetEdgeBasisKey(SegGeomSharedPtr edge, const std::string variable = "DefaultVar");

        protected:
            void ReadEdges    (TiXmlDocument &doc);
            void ReadElements (TiXmlDocument &doc);
            void ReadComposites(TiXmlDocument &doc);
            void ResolveGeomRef(const std::string &prevToken, const std::string &token,
                    Composite& composite);

#ifdef OLD
            bool   m_geoFacDefined;
#endif

        private:

        };

        typedef boost::shared_ptr<MeshGraph2D> MeshGraph2DSharedPtr;

    };
};

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH2D_H

