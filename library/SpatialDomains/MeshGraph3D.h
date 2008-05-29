////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph3D.h,v $
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
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>

#include <list>

namespace Nektar
{
    namespace SpatialDomains
    {
        class SegGeom;
        class TriGeom;
        class QuadGeom;
        class TetGeom;
        class PyrGeom;
        class PrismGeom;
        class HexGeom;

        class MeshGraph3D: public MeshGraph
        {
        public:

            MeshGraph3D();
            virtual ~MeshGraph3D();

            void ReadGeometry(std::string &infilename);
            void ReadGeometry(TiXmlDocument &doc);

            SegGeomSharedPtr GetSegGeom(int eID);
            TriGeomSharedPtr GetTriGeom(int tID);

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

            void GenXGeoFac();

            inline const int GetNseggeoms() const 
            {
                return int(m_seggeoms.size());
            }

            inline const int GetVidFromElmt(StdRegions::ShapeType shape, 
                const int vert, const int elmt) const 
            {
                if(shape == StdRegions::eTriangle)
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_trigeoms.size()),
                        "eid is out of range");

                    return m_trigeoms[elmt]->GetVid(vert);
                }
                else
                {
                    ASSERTL2((elmt >=0)&&(elmt < m_quadgeoms.size()),
                        "eid is out of range");

                    return m_quadgeoms[elmt]->GetVid(vert);
                }
            }

            inline const int GetEidFromElmt(StdRegions::ShapeType shape, 
                const int edge, const int elmt) const
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

            inline const StdRegions::EdgeOrientation GetEorientFromElmt(StdRegions::ShapeType shape,const int edge, const int elmt) const 
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


            inline const StdRegions::EdgeOrientation GetCartesianEorientFromElmt(StdRegions::ShapeType shape,const int edge, const int elmt) const
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
                return int(m_MeshCompositeVector.size());
            }

            int GetNumCompositeItems(int whichComposite)
            {
                int returnval = -1;

                try
                {
                    returnval = int(m_MeshCompositeVector[whichComposite]->size());
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
            ElementEdgeVectorSharedPtr GetElementsFromEdge(SegGeomSharedPtr edge);

        protected:
            void ReadEdges    (TiXmlDocument &doc);
            void ReadFaces    (TiXmlDocument &doc);
            void ReadElements (TiXmlDocument &doc);
            void ReadComposites(TiXmlDocument &doc);
            void ResolveGeomRef(const std::string &prevToken, const std::string &token);

        private:
            SegGeomVector           m_seggeoms;
            TriGeomVector           m_trigeoms;
            QuadGeomVector          m_quadgeoms;
            TetGeomVector           m_tetgeoms;
            PyrGeomVector           m_pyrgeoms;
            PrismGeomVector         m_prismgeoms;
            HexGeomVector           m_hexgeoms;
        };
    }; // end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H

//
// $Log: MeshGraph3D.h,v $
// Revision 1.4  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.3  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.2  2006/07/02 17:16:17  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.6  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.5  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.4  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.3  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
