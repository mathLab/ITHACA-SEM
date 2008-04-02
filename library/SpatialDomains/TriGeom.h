////////////////////////////////////////////////////////////////////////////////
//
//  File:  TriGeom.h
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
#ifndef NEKTAR_SPATIALDOMAINS_TRIGEOM_H
#define NEKTAR_SPATIALDOMAINS_TRIGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/TriFaceComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class TriGeom;
        typedef boost::shared_ptr<TriGeom> TriGeomSharedPtr;
        typedef std::vector< TriGeomSharedPtr > TriGeomVector;
        typedef std::vector< TriGeomSharedPtr >::iterator TriGeomVectorIter;

        class TriGeom: public TriFaceComponent
        {
            public:
                TriGeom();
                TriGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
                TriGeom(const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
                ~TriGeom();

                inline void SetOwnData()
                {
                    m_owndata = true; 
                }

                void FillGeom();
    
                void GetLocCoords(const ConstArray<OneD,NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);

                inline int GetEid(int i) const
                {
                    ASSERTL2((i >=0) && (i <= 2),"Edge id must be between 0 and 2");
                    return m_edges[i]->GetEid();
                }

                inline int GetVid(const int i) const
                {
                    ASSERTL2((i >=0) && (i <= 2),"Vertex id must be between 0 and 2");
                    return m_verts[i]->GetVid();
                }

                inline StdRegions::EdgeOrientation GetEorient(const int i) const
                {
                    ASSERTL2((i >=0) && (i <= 2),"Edge id must be between 0 and 2");
                    return m_eorient[i];
                }

                /// \brief Return the edge number of the given edge, or -1, if 
                /// not an edge of this element.
                int WhichEdge(SegGeomSharedPtr edge)
                {
                    int returnval = -1;

                    SegGeomVector::iterator edgeIter;
                    int i;

                    for (i=0,edgeIter = m_edges.begin(); edgeIter != m_edges.end(); ++edgeIter,++i)
                    {
                        if (*edgeIter == edge)
                        {
                            returnval = i;
                            break;
                        }
                    }

                    return returnval;
                }

                static const int kNedges = 3;
                static const int kNverts = 3;

            protected:
                VertexComponentVector           m_verts;
                SegGeomVector                   m_edges;
                StdRegions::EdgeOrientation     m_eorient [kNedges];

                void GenGeomFactors(void);
            private:
                bool m_owndata;

                virtual void v_GenGeomFactors(void)
                {
                    GenGeomFactors();
                }

                virtual void v_SetOwnData()
                {
                    SetOwnData();
                }

                virtual void v_FillGeom()
                {
                    FillGeom();
                }            
                
                virtual void v_GetLocCoords(const ConstArray<OneD,NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
                {
                    v_GetLocCoords(coords,Lcoords);
                }
                
                virtual int v_GetEid(int i) const
                {
                    return GetEid(i);
                }
                
                virtual int v_GetVid(int i) const
                {
                    return GetVid(i);
                }
                
                virtual const SegGeomSharedPtr v_GetEdge(int i) const
                {
                    return GetEdge(i);
                }
                
                virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
                {
                    return GetEorient(i);
                }
                
                virtual int v_WhichEdge(SegGeomSharedPtr edge)
                {
                    return WhichEdge(edge);
                }
        };
    }; //end of namespace SpatialDomains
}; //end of namespace Nektar

#endif //NEKTAR_SPATIALDOMAINS_TRIGEOM_H

//
// $Log: TriGeom.h,v $
// Revision 1.15  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.14  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.13  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.12  2007/12/11 21:51:53  jfrazier
// Updated 2d components so elements could be retrieved from edges.
//
// Revision 1.11  2007/07/22 23:04:24  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.10  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.9  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.8  2007/06/01 17:08:07  pvos
// Modification to make LocalRegions/Project2D run correctly (PART1)
//
// Revision 1.7  2007/05/28 21:48:43  sherwin
// Update for 2D functionality
//
// Revision 1.6  2006/07/02 17:16:18  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.5  2006/06/01 14:15:31  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.3  2006/05/29 17:05:17  sherwin
// Updated to use shared_ptr around Geom types - added typedef
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:05  kirby
// *** empty log message ***
//
// Revision 1.23  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.22  2006/04/09 02:08:36  jfrazier
// Added precompiled header.
//
// Revision 1.21  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.20  2006/03/12 14:20:44  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.19  2006/03/12 11:06:41  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.18  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.17  2006/02/19 01:37:35  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
