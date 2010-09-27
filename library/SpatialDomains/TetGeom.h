////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/TetGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_TETGEOM
#define NEKTAR_SPATIALDOMAINS_TETGEOM

 #include <StdRegions/StdRegions.hpp>
 #include <StdRegions/StdTetExp.h>

 #include <SpatialDomains/Geometry3D.h>
 #include <SpatialDomains/TriGeom.h>
 #include <SpatialDomains/GeomFactors3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class TetGeom;

        typedef boost::shared_ptr<TetGeom> TetGeomSharedPtr;
        typedef std::vector< TetGeomSharedPtr > TetGeomVector;
        typedef std::vector< TetGeomSharedPtr >::iterator TetGeomVectorIter;

        class TetGeom: public LibUtilities::GraphVertexObject, public Geometry3D
        {
        public:
			TetGeom ();

            TetGeom(const TriGeomSharedPtr faces[]);
			TetGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const TriGeomSharedPtr faces[],
					const StdRegions::EdgeOrientation eorient[], const StdRegions::FaceOrientation forient[]);

			TetGeom(const VertexComponentSharedPtr verts[], const TriGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
			TetGeom(const TriGeomSharedPtr faces[], const StdRegions::FaceOrientation forient[]);
			TetGeom(const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
			~TetGeom();

            void AddElmtConnected(int gvo_id, int locid);
            int  NumElmtConnected() const;
            bool IsElmtConnected(int gvo_id, int locid) const;
            void FillGeom();
            void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);

            inline int GetFid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_faces[i]->GetFid();
            }

            inline StdRegions::FaceOrientation GetFaceorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_forient[i];
            }

            inline const TriGeomSharedPtr GetFace(int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_faces[i];
            }

            inline int GetEid() const
            {
                return m_eid;
            }

            inline int GetEid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 5),"Edge id must be between 0 and 5");
                return m_edges[i]->GetEid();
            }

            inline StdRegions::EdgeOrientation GetEorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 5),"Edge id must be between 0 and 5");
                return m_eorient[i];
            }

            inline const SegGeomSharedPtr GetEdge(int i) const
            {
                ASSERTL2((i >=0) && (i <= 5),"Edge id must be between 0 and 5");
                return m_edges[i];
            }

            inline int GetVid(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Vertex id must be between 0 and 3");
                return m_verts[i]->GetVid();
            }


            /// \brief Return the face number of the given face, or -1, if
            /// not an face of this element.
            int WhichFace(Geometry2DSharedPtr face)
            {
                int returnval = -1;

                TriGeomVector::iterator faceIter;
                int i;

                for (i=0,faceIter = m_faces.begin(); faceIter != m_faces.end(); ++faceIter,++i)
                {
                    if (*faceIter == face)
                    {
                        returnval = i;
                        break;
                    }
                }

                return returnval;
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


            inline int GetCoordDim() const
            {
                return m_coordim;
            }

            inline void SetOwnData()
            {
                m_ownData = true;
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

            static const int                    kNverts = 4;
            static const int                    kNedges = 6;
            static const int                    kNqfaces = 0;
            static const int                    kNtfaces = 4;
            static const int                    kNfaces = kNqfaces + kNtfaces;

        protected:
            VertexComponentVector               m_verts;
            SegGeomVector                       m_edges;
            TriGeomVector                       m_faces;
            StdRegions::EdgeOrientation         m_eorient[kNedges];
            StdRegions::FaceOrientation         m_forient[kNfaces];
            int                                 m_eid;
            bool                                m_ownVerts;
            std::list<CompToElmt>               m_elmtMap;

            Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;
            void GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

        private:
            bool                                m_ownData;

            void SetUpLocalEdges();
            void SetUpLocalVertices();
            void SetUpEdgeOrientation();
            void SetUpFaceOrientation();

            virtual void v_GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
            {
                GenGeomFactors(tbasis);
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

            virtual int v_GetFid(int i) const
            {
               return GetFid(i);
            }

            virtual const Geometry2DSharedPtr v_GetFace(int i) const
            {
               return GetFace(i);
            }

            virtual int v_GetEid(int i) const
            {
                return GetEid(i);
            }

            virtual const Geometry1DSharedPtr v_GetEdge(int i) const
            {
                return GetEdge(i);
            }

            virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
            {
                return GetEorient(i);
            }

           virtual StdRegions::FaceOrientation v_GetFaceorient(const int i) const
            {
                return GetFaceorient(i);
            }

            virtual int v_GetVid(int i) const
            {
                return GetVid(i);
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

            virtual int v_WhichFace(Geometry2DSharedPtr face)
            {
                return WhichFace(face);
            }

            virtual int v_GetNumVerts() const
            {
                return kNverts;
            }

            virtual int v_GetNumEdges() const
            {
                return kNedges;
            }

            virtual bool v_ContainsPoint(
                    const Array<OneD, const NekDouble> &gloCoord);

        };

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_TETGEOM

//
// $Log: TetGeom.h,v $
// Revision 1.17  2009/01/21 16:59:04  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.16  2008/11/17 08:58:32  ehan
// Added necessary mapping routines for Tet
//
// Revision 1.15  2008/08/26 02:18:17  ehan
// Changed shared pointer (TriGeomSharedPtr to the Geometry2DSharedPtr) in the face function.
//
// Revision 1.14  2008/06/30 19:35:52  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.13  2008/06/16 22:44:23  ehan
// Added virtual function GetFaceorientation(..)
//
// Revision 1.12  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.11  2008/05/12 17:30:09  ehan
// Added virtual functions
//
// Revision 1.10  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.9  2008/04/02 22:19:04  pvos
// Update for 2D local to global mapping
//
// Revision 1.8  2008/03/25 08:43:16  ehan
// Added constructor, GetLocCoords(), and FillGeom().
//
// Revision 1.7  2008/02/27 22:31:29  ehan
// Added inline function "SetOwnData()".
//
// Revision 1.6  2008/02/12 01:28:58  ehan
// Included TriGeom.h to prevent undefined TriGeom error.
//
// Revision 1.5  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.4  2008/02/05 00:43:22  ehan
// Included geometry3D and meshgraphics inorder to prevent compile error.
//
// Revision 1.3  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.16  2006/04/09 02:08:36  jfrazier
// Added precompiled header.
//
// Revision 1.15  2006/03/12 14:20:44  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.14  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.13  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.12  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
