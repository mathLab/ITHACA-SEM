////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/PrismGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
#define NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/GeomFactors3D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        class PrismGeom;
        typedef boost::shared_ptr<PrismGeom> PrismGeomSharedPtr;
        typedef std::vector< PrismGeomSharedPtr > PrismGeomVector;
        typedef std::vector< PrismGeomSharedPtr >::iterator PrismGeomVectorIter;

        class PrismGeom: public LibUtilities::GraphVertexObject, public Geometry3D
        {
        public:
            SPATIAL_DOMAINS_EXPORT PrismGeom ();
            SPATIAL_DOMAINS_EXPORT PrismGeom(const Geometry2DSharedPtr faces[]);
            SPATIAL_DOMAINS_EXPORT PrismGeom(const TriGeomSharedPtr tfaces[], const QuadGeomSharedPtr qfaces[], const StdRegions::FaceOrientation forient[]);
            SPATIAL_DOMAINS_EXPORT PrismGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const TriGeomSharedPtr tfaces[],
                      const QuadGeomSharedPtr qfaces[],const StdRegions::EdgeOrientation eorient[],
                      const StdRegions::FaceOrientation forient[]);
            SPATIAL_DOMAINS_EXPORT PrismGeom(const Geometry2DSharedPtr faces[], const StdRegions::FaceOrientation forient[]);
            SPATIAL_DOMAINS_EXPORT ~PrismGeom();
            
            SPATIAL_DOMAINS_EXPORT void AddElmtConnected(int gvo_id, int locid);
            SPATIAL_DOMAINS_EXPORT int  NumElmtConnected() const;
            SPATIAL_DOMAINS_EXPORT bool IsElmtConnected(int gvo_id, int locid) const;
			SPATIAL_DOMAINS_EXPORT void FillGeom();
			SPATIAL_DOMAINS_EXPORT void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);

            inline int GetFid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 4),"Edge id must be between 0 and 4");
                return m_faces[i]->GetFid();
            }

            inline StdRegions::FaceOrientation GetFaceorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 4),"Edge id must be between 0 and 4");
                return m_forient[i];
            }

            inline const Geometry2DSharedPtr GetFace(int i) const
            {
                ASSERTL2((i >=0) && (i <= 4),"Edge id must be between 0 and 4");
                return m_faces[i];
            }

            inline int GetEid() const 
            {
                return m_eid;
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

            inline void SetOwnData()
            {
                m_owndata = true; 
            }

            /// \brief Return the face number of the given face, or N, if not a face of this element.
            int WhichFace(Geometry2DSharedPtr face)
            {
                int i = 0;

                Geometry2DVector::iterator f;
                for (i = 0, f = m_faces.begin(); f != m_faces.end(); ++f,++i)
                {
                    if (*f == face)
                    {
                        break;
                    }
                }
                return i;
            }

            SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord);

            SPATIAL_DOMAINS_EXPORT static const int kNverts = 6;
            SPATIAL_DOMAINS_EXPORT static const int kNedges = 9;
            SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 3;
            SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 2;
            SPATIAL_DOMAINS_EXPORT static const int kNfaces = kNqfaces + kNtfaces;
            SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;
        protected:
            bool m_owndata;

            void SetUpLocalEdges();
            void SetUpLocalVertices();
            void SetUpEdgeOrientation();
            void SetUpFaceOrientation();
            
            VertexComponentVector           m_verts;
            SegGeomVector                   m_edges;
            StdRegions::EdgeOrientation     m_eorient [kNedges];
            StdRegions::FaceOrientation     m_forient[kNfaces];
            Geometry2DVector                m_faces;
            TriGeomVector                   m_tfaces;
            QuadGeomVector                  m_qfaces;

            int m_eid;

            bool m_ownverts;
            std::list<CompToElmt> m_elmtmap;

            Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;
			void GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

        private:
 

            virtual void v_GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
            {
                GenGeomFactors(tbasis);
            }

            virtual void v_FillGeom()
            {
                FillGeom();
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

            virtual void v_SetOwnData()
            {
                SetOwnData();   
            }
                        
            virtual int v_GetFid(int i) const
            {
               return GetFid(i);
            }

            virtual const Geometry2DSharedPtr v_GetFace(int i) const
            {
               return GetFace(i);
            }

            virtual StdRegions::FaceOrientation v_GetFaceorient(const int i) const
            {
                return GetFaceorient(i);
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

        };
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

//
// $Log: PrismGeom.h,v $
// Revision 1.13  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.12  2008/11/17 08:59:29  ehan
// Added necessary mapping routines for Tet
//
// Revision 1.11  2008/06/30 19:35:03  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.10  2008/06/16 22:42:18  ehan
// Added inline function GetFace(..), whichFace(..), and GetFaceorient(..).
//
// Revision 1.9  2008/06/14 01:22:38  ehan
// Implemented constructor and FillGeom().
//
// Revision 1.8  2008/06/12 21:22:55  delisi
// Added method stubs for GenGeomFactors, FillGeom, and GetLocCoords.
//
// Revision 1.7  2008/06/11 21:34:42  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.6  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.5  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.4  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.3  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.2  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.15  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.14  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 11:06:39  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.12  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.11  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

