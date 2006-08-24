////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshGraph2D.h,v $
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
#include <SpatialDomains/EdgeComponent.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>

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

            MeshGraph2D();
            virtual ~MeshGraph2D();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);
            void Write(std::string &outfilename);

            EdgeComponentSharedPtr GetEdgeComponent(int eID);

            inline const int GetCoordim(void){
                return GetSpaceDimension();
            }

            inline const TriGeomVector &GetTrigeoms(void) const
            {
                return m_trigeoms;
            }

            inline const QuadGeomVector &GetQuadgeoms(void) const
            {
                return m_quadgeoms;
            }

            inline bool GetGeofac_defined(void)
            {
                return m_geofac_defined;
            }

            void GenXGeoFac();

            inline int GetNecomps()
            {
                return int(m_ecomps.size());
            }

            inline int GetVidFromElmt(StdRegions::ShapeType shape, 
                const int vert, const int elmt)
            {
                if(shape == StdRegions::eTriangle)
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    return m_trigeoms[elmt]->GetVid(vert);
                }
                else
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    return m_quadgeoms[elmt]->GetVid(vert);
                }
            }

            inline int GetEidFromElmt(StdRegions::ShapeType shape, 
                const int edge, const int elmt)
            {
                if(shape == StdRegions::eTriangle)
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    return m_trigeoms[elmt]->GetEid(edge);
                }
                else
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_quadgeoms.size()),
                        "eid is out of range");

                    return m_quadgeoms[elmt]->GetEid(edge);
                }
            }

            inline StdRegions::EdgeOrientation GetEorientFromElmt(StdRegions::ShapeType shape,const int edge, const int elmt)
            {
                if(shape == StdRegions::eTriangle)
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    return m_trigeoms[elmt]->GetEorient(edge);
                }
                else
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_quadgeoms.size()),
                        "eid is out of range");

                    return m_quadgeoms[elmt]->GetEorient(edge);
                }
            }


            inline StdRegions::EdgeOrientation GetCartesianEorientFromElmt(StdRegions::ShapeType shape,const int edge, const int elmt)
            {
                StdRegions::EdgeOrientation returnval;

                if(shape == StdRegions::eTriangle)
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    returnval = m_trigeoms[elmt]->GetEorient(edge);
                }
                else
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_quadgeoms.size()),
                        "eid is out of range");

                    returnval =  m_quadgeoms[elmt]->GetEorient(edge);
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
                return int(m_Mesh2DCompositeVector.size());
            }

            int GetNumCompositeItems(int whichComposite)
            {
                int returnval = -1;

                try
                {
                    returnval = int(m_Mesh2DCompositeVector[whichComposite].size());
                }
                catch(...)
                {
                    std::ostrstream errStream;
                    errStream << "Unable to access composite item [" << whichComposite << "].";
                    NEKERROR(ErrorUtil::efatal, errStream.str());
                }

                return returnval;
            }

            GeometrySharedPtr GetCompositeItem(int whichComposite, int whichItem);

        protected:
            void ReadEdges    (TiXmlDocument &doc);
            void ReadElements (TiXmlDocument &doc);
            void ReadComposites(TiXmlDocument &doc);
            GeometrySharedPtr ResolveGeomRef(const std::string &prevToken, const std::string &token);

            bool   m_geofac_defined;

        private:
            EdgeComponentVector m_ecomps;
            TriGeomVector       m_trigeoms;
            QuadGeomVector      m_quadgeoms;

            std::vector< std::vector<GeometrySharedPtr> > m_Mesh2DCompositeVector;
        };
    };
};

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH2D_H

//
// $Log: MeshGraph2D.h,v $
// Revision 1.6  2006/08/18 19:45:29  jfrazier
// Completed composites.
//
// Revision 1.5  2006/08/17 22:55:00  jfrazier
// Continued adding code to process composites in the mesh2d.
//
// Revision 1.4  2006/08/16 23:34:42  jfrazier
// *** empty log message ***
//
// Revision 1.3  2006/07/02 17:16:17  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.2  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.1  2006/05/04 18:59:01  kirby
// *** empty log message ***
//
// Revision 1.11  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.10  2006/04/04 23:12:38  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.9  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.8  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.7  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.6  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.5  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
