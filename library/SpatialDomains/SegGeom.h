////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.h
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
#ifndef NEKTAR_SPATIALDOMAINS_SEGGEOM_H
#define NEKTAR_SPATIALDOMAINS_SEGGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/EdgeComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class SegGeom: public EdgeComponent
        {
            public:
                SegGeom();

                SegGeom(int id, const int coordim):
                    EdgeComponent(id,coordim)
                    {
                    }

                SegGeom(int id, const int coordim,
                        const VertexComponentSharedPtr vertex[]): 
                    EdgeComponent(id,coordim,vertex)
                    {
                    }

                SegGeom(const int id, const VertexComponentSharedPtr vert1, const VertexComponentSharedPtr  vert2);

                SegGeom(const SegGeom &in);

                ~SegGeom();

                inline int GetVid(int i) const
                {
                    ASSERTL2((i >=0) && (i <= 1),"Verted id must be between 0 and 1");
                    return m_verts[i]->GetVid();
                }

                void    FillGeom ();

                StdRegions::ShapeType DetShapeType() const
                {
                    return StdRegions::eSegment;
                }

                void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);


                inline void SetOwnData()
                {
                    m_owndata = true; 
                }

                void WriteToFile(std::ofstream &outfile, const int dumpVar);

            protected:

                void GenGeomFactors(void);

            private:
                bool m_owndata;   ///< Boolean indicating whether object owns the data

                virtual void v_GenGeomFactors(void)
                {
                    GenGeomFactors();
                }

                virtual int v_GetVid(int i) const
                {
                    return GetVid(i);
                }
                
                virtual void v_FillGeom()
                {
                    FillGeom();
                }
                
                virtual StdRegions::ShapeType v_DetShapeType() const
                {
                    return DetShapeType();
                }
                
                virtual void v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
                {
                    GetLocCoords(coords,Lcoords);
                }            
                
                virtual  void v_SetOwnData()
                {
                    SetOwnData();
                }
                
                virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
                {
                    WriteToFile(outfile, dumpVar);
                }
        };
        
        // shorthand for boost pointer
        typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< SegGeomSharedPtr >::iterator SegGeomVectorIter;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_SEGGEOM_H

//
// $Log: SegGeom.h,v $
// Revision 1.15  2008/04/02 22:19:04  pvos
// Update for 2D local to global mapping
//
// Revision 1.14  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.13  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.12  2007/07/22 23:04:24  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.11  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.10  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.9  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.8  2007/03/20 09:17:40  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.7  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
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
// Revision 1.2  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.25  2006/04/04 23:12:38  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.24  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.23  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.22  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.21  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.20  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
