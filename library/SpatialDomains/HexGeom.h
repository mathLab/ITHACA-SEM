////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/HexGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_HEXGEOM_H
#define NEKTAR_SPATIALDOMAINS_HEXGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdHexExp.h>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/QuadGeom.h>


namespace Nektar
{
    namespace SpatialDomains
    {
	    class HexGeom;

        typedef boost::shared_ptr<HexGeom> HexGeomSharedPtr;
        typedef std::vector< HexGeomSharedPtr > HexGeomVector;
        typedef std::vector< HexGeomSharedPtr >::iterator HexGeomVectorIter;

        class HexGeom: public LibUtilities::GraphVertexObject, public Geometry3D
        {
        public:
            HexGeom();
            HexGeom(const QuadGeomSharedPtr faces[],  const StdRegions::FaceOrientation forient[]);
            HexGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const QuadGeomSharedPtr faces[],
                    const StdRegions::EdgeOrientation eorient[], const StdRegions::FaceOrientation forient[]);
            ~HexGeom();

            void AddElmtConnected(int gvo_id, int locid);
            int  NumElmtConnected() const;
            bool IsElmtConnected(int gvo_id, int locid) const;
			void FillGeom();
            void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);

            inline void SetOwnData()
            {
               m_owndata = true;
            }

            inline int GetEid() const 
            {
                return m_eid;
            }

            inline int GetVid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 7),"Verted id must be between 0 and 7");
                return m_verts[i]->GetVid();
            }

            inline const SegGeomSharedPtr GetEdge(int i) const
            {
                ASSERTL2((i >=0) && (i <= 11),"Edge id must be between 0 and 11");
                return m_edges[i];
            }

            inline StdRegions::EdgeOrientation GetEorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 11),"Edge id must be between 0 and 11");
                return m_eorient[i];
            }

            inline int GetCoordDim() const 
            {
                return m_coordim;
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
            {
                return m_xmap[i]->GetBasis(j);
            }

            inline Array<OneD,NekDouble> &UpdatePhys(const int i)
            {
                return m_xmap[i]->UpdatePhys();
            }

            NekDouble GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord);

            
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


            static const int kNverts = 8;
            static const int kNedges = 12;
            static const int kNqfaces = 6;
            static const int kNtfaces = 0;
            static const int kNfaces = kNqfaces + kNtfaces;

       protected:

            VertexComponentVector           m_verts;
            SegGeomVector                   m_edges;
            QuadGeomVector                  m_qfaces;
            StdRegions::EdgeOrientation     m_eorient[kNedges];
            StdRegions::FaceOrientation     m_forient[kNfaces];

            int m_eid;
            bool m_ownverts;
            std::list<CompToElmt> m_elmtmap;

            Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;
			void GenGeomFactors(void);
        
        private:

            bool m_owndata;

            virtual void v_GenGeomFactors(void)
            {
                GenGeomFactors( );
            }

            virtual void v_SetOwnData()
            {
                SetOwnData();
            }

            virtual void v_FillGeom()
            {
                FillGeom();
            }

            virtual void v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                GetLocCoords(coords,Lcoords);
            }
            
            virtual void v_AddElmtConnected(int gvo_id, int locid)
            {
                AddElmtConnected(gvo_id,locid);
            }

            virtual int  v_NumElmtConnected() const
            {
                return NumElmtConnected();
            }

            virtual bool v_IsElmtConnected(int gvo_id, int locid) const
            {
                return IsElmtConnected(gvo_id,locid);
            }
            
            virtual int v_GetEid() const 
            {
                return GetEid();
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

            virtual int v_GetCoordDim() const 
            {
                return GetCoordDim();
            }

            virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
            {
                return GetBasis(i,j);
            }

            virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
            {
                return UpdatePhys(i);
            }

            virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                return GetCoord(i,Lcoord);
            }

            virtual int v_WhichEdge(SegGeomSharedPtr edge)
            {
                return WhichEdge(edge);
            }

        };

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_HEXGEOM_H

//
// $Log: HexGeom.h,v $
// Revision 1.12  2008/06/12 21:22:55  delisi
// Added method stubs for GenGeomFactors, FillGeom, and GetLocCoords.
//
// Revision 1.11  2008/06/11 21:34:41  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.10  2008/05/12 17:28:16  ehan
// Added virtual functions
//
// Revision 1.9  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.8  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.7  2008/02/12 01:26:00  ehan
// Included stdExpansion3D to prevent undefined error  of "StdExpansion3DSharedPtr".
//
// Revision 1.6  2008/02/10 01:05:57  jfrazier
// Changed  include order.
//
// Revision 1.5  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.4  2008/02/03 05:05:06  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.3  2008/01/31 11:00:56  ehan
// Added boost pointer  and included MeshGraph.h
//
// Revision 1.2  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.16  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.15  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.14  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.13  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
