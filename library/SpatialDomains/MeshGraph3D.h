////////////////////////////////////////////////////////////////////////////////
//
//  File:  MeshGraph3D.h
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
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <unordered_map>

namespace Nektar
{
    namespace SpatialDomains
    {
        class SegGeom;
        class QuadGeom;
        class TetGeom;
        class PyrGeom;
        class PrismGeom;
        class HexGeom;

        class MeshGraph3D: public MeshGraph
        {
        public:

            SPATIAL_DOMAINS_EXPORT MeshGraph3D();
            SPATIAL_DOMAINS_EXPORT MeshGraph3D(const LibUtilities::SessionReaderSharedPtr &pSession, const DomainRangeShPtr &rng = NullDomainRangeShPtr);
            SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph3D();

            SPATIAL_DOMAINS_EXPORT void ReadGeometry(const std::string &infilename);
            SPATIAL_DOMAINS_EXPORT void ReadGeometry(TiXmlDocument &doc);

            SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GetSegGeom(int eID);
            SPATIAL_DOMAINS_EXPORT Geometry2DSharedPtr GetGeometry2D(int gID);

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

            SPATIAL_DOMAINS_EXPORT void GenXGeoFac();

            inline int GetNseggeoms() const
            {
                return int(m_segGeoms.size());
            }

            inline int GetVidFromElmt(LibUtilities::ShapeType shape,
                const int vert, const int elmt) const
            {
                if(shape == LibUtilities::eTriangle)
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

            inline int GetEidFromElmt(LibUtilities::ShapeType shape,
                const int edge, const int elmt) const
            {
                if(shape == LibUtilities::eTriangle)
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

            inline StdRegions::Orientation GetEorientFromElmt(LibUtilities::ShapeType shape,const int edge, const int elmt) const
            {
                if(shape == LibUtilities::eTriangle)
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


            inline StdRegions::Orientation GetCartesianEorientFromElmt(LibUtilities::ShapeType shape,const int edge, const int elmt) const
            {
                StdRegions::Orientation returnval;

                if(shape == LibUtilities::eTriangle)
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

            /// \brief Return the elements (shared ptrs) that have this face.
            SPATIAL_DOMAINS_EXPORT ElementFaceVectorSharedPtr GetElementsFromFace(Geometry2DSharedPtr face);

            /// \brief Return the BasisKey corresponding to a face of an element
            SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKey GetFaceBasisKey(Geometry2DSharedPtr face, const int flag, const std::string variable = "DefaultVar");

        protected:
            void ReadEdges    (TiXmlDocument &doc);
            void ReadFaces    (TiXmlDocument &doc);
            void ReadElements (TiXmlDocument &doc);
            void ReadComposites(TiXmlDocument &doc);
            void ResolveGeomRef(const std::string &prevToken, const std::string &token,
                    Composite& composite);

        private:
            void PopulateFaceToElMap(Geometry3DSharedPtr element, int kNfaces);
            std::unordered_map<int, ElementFaceVectorSharedPtr> m_faceToElMap;

        };

        typedef std::shared_ptr<MeshGraph3D> MeshGraph3DSharedPtr;

    }; // end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H
