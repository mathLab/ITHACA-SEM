////////////////////////////////////////////////////////////////////////////////
//
//  File: QuadGeom.h
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
#ifndef NEKTAR_SPATIALDOMAINS_QUADGEOM_H
#define NEKTAR_SPATIALDOMAINS_QUADGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdQuadExp.h>

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/MeshComponents.h>


namespace Nektar
{
    namespace SpatialDomains
    {
        class QuadGeom;
        typedef boost::shared_ptr<QuadGeom> QuadGeomSharedPtr;
        typedef std::vector< QuadGeomSharedPtr > QuadGeomVector;
        typedef std::vector< QuadGeomSharedPtr >::iterator QuadGeomVectorIter;
    
        class QuadGeom: public Geometry2D
        {        
        public:
            QuadGeom();
            QuadGeom(int id, const int coordim);
            QuadGeom(const int id, const VertexComponentSharedPtr verts[],  const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
            QuadGeom(const int id, const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[]);
            QuadGeom(const QuadGeom &in);
			~QuadGeom();

			void AddElmtConnected(int gvo_id, int locid);
			int  NumElmtConnected() const;
			bool IsElmtConnected(int gvo_id, int locid) const;

			inline int GetFid() const 
			{
				return m_fid;
			}

			inline int GetCoordDim() const 
			{
				return m_coordim;
			}

			inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
			{
				return m_xmap[i]->GetBasis(j);
			}
                        
                        inline const LibUtilities::BasisSharedPtr GetEdgeBasis(const int i, const int j)
                        {
                            ASSERTL1(j <= 3,"edge is out of range");
                            if((j == 0)||(j == 2))
                            {
                                return m_xmap[i]->GetBasis(0);
                            }
                            else
                            {
                                return m_xmap[i]->GetBasis(1);
                            }
                        }

			inline Array<OneD,NekDouble> &UpdatePhys(const int i)
			{
				return m_xmap[i]->UpdatePhys();
			}

			NekDouble GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord);

            inline void SetOwnData()
            {
                m_owndata = true; 
            }

            void FillGeom();

            void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);

            inline int GetEid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_edges[i]->GetEid();
            }

            inline int GetVid(int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Verted id must be between 0 and 3");
                return m_verts[i]->GetVid();
            }
            
            inline const VertexComponentSharedPtr GetVertex(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Vertex id must be between 0 and 3");
                return m_verts[i];
            }
            
            inline const Geometry1DSharedPtr GetEdge(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_edges[i];
            }

            inline StdRegions::EdgeOrientation GetEorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                return m_eorient[i];            
            }

            /// \brief Get the orientation of face1.
            ///
            static StdRegions::FaceOrientation GetFaceOrientation(const QuadGeom &face1,
                                                                  const QuadGeom &face2);
                                                                      

            inline StdRegions::EdgeOrientation GetCartesianEorient(const int i) const
            {
                ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
                if(i < 2)
                {
                    return m_eorient[i];
                }
                else
                {
                    if(m_eorient[i] == StdRegions::eForwards)
                    {
                        return StdRegions::eBackwards;
                    }
                    else
                    {
                        return StdRegions::eForwards; 
                    }
                }
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


            static const int kNverts = 4;
            static const int kNedges = 4;

        protected:
            VertexComponentVector           m_verts;
            SegGeomVector                   m_edges;
            StdRegions::EdgeOrientation     m_eorient[kNedges];
            int                             m_fid;
            bool                            m_ownverts;
            std::list<CompToElmt>           m_elmtmap;

            void GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

        private:
            bool m_owndata;   ///< Boolean indicating whether object owns the data

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
            
			virtual int v_GetFid() const 
			{
				return GetFid();
			}

			virtual int v_GetCoordDim() const 
			{
				return GetCoordDim();
			}

			virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
			{
				return GetBasis(i,j);
			}

			virtual const LibUtilities::BasisSharedPtr v_GetEdgeBasis(const int i, const int j)
			{
				return GetEdgeBasis(i,j);
			}

			virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
			{
				return UpdatePhys(i);
			}

			virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
			{
				return GetCoord(i,Lcoord);
			}

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
            
            virtual int v_GetEid(int i) const
            {
                return GetEid(i);
            }
            
            virtual int v_GetVid(int i) const
            {
                return GetVid(i);
            }
            
            virtual const VertexComponentSharedPtr v_GetVertex(int i) const
            {
                return GetVertex(i);
            }
            
            virtual const Geometry1DSharedPtr v_GetEdge(int i) const
            {
                return GetEdge(i);
            }
            
            virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
            {
                return GetEorient(i);
            }

            virtual StdRegions::EdgeOrientation v_GetCartesianEorient(const int i) const
            {
                return GetCartesianEorient(i);
            }
            
            virtual int v_WhichEdge(SegGeomSharedPtr edge)
            {
                return WhichEdge(edge);
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

#endif //NEKTAR_SPATIALDOMAINS_QUADGEOM_H

//
// $Log: QuadGeom.h,v $
// Revision 1.25  2008/11/17 08:59:01  ehan
// Added GetNumVerts and GetNumEdges
//
// Revision 1.24  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.23  2008/07/29 22:23:36  sherwin
// various mods for DG advection solver in Multiregions. Added virtual calls to Geometry, Geometry1D, 2D and 3D
//
// Revision 1.22  2008/06/30 19:35:36  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.21  2008/06/14 01:23:17  ehan
// Implemented constructor and FillGeom().
//
// Revision 1.20  2008/06/11 21:34:42  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.19  2008/05/10 18:27:33  sherwin
// Modifications necessary for QuadExp Unified DG Solver
//
// Revision 1.18  2008/05/09 20:33:16  ehan
// Fixed infinity loop of recursive function.
//
// Revision 1.17  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.16  2008/04/02 22:19:04  pvos
// Update for 2D local to global mapping
//
// Revision 1.15  2008/01/26 20:17:28  sherwin
// EdgecComponent to SegGeom
//
// Revision 1.14  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.13  2007/12/11 21:51:53  jfrazier
// Updated 2d components so elements could be retrieved from edges.
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
// Revision 1.9  2007/06/01 17:08:07  pvos
// Modification to make LocalRegions/Project2D run correctly (PART1)
//
// Revision 1.8  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.7  2006/07/02 17:16:18  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.6  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.5  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.4  2006/05/29 17:05:17  sherwin
// Updated to use shared_ptr around Geom types - added typedef
//
// Revision 1.3  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.2  2006/05/06 20:36:16  sherwin
// Modifications to get LocalRegions/Project1D working
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.24  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.23  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.22  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.21  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.20  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.19  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.18  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
