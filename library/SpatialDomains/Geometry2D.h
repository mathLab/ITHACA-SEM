////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry2D.h
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry2D;
        
        // shorthand for boost pointer
        typedef boost::shared_ptr<Geometry2D> Geometry2DSharedPtr;
        typedef std::vector< Geometry2DSharedPtr > Geometry2DVector;
        typedef std::vector< Geometry2DSharedPtr >::iterator Geometry2DVectorIter;
        
        class Geometry2D: public Geometry
        {       
        public:
            SPATIAL_DOMAINS_EXPORT Geometry2D();
            SPATIAL_DOMAINS_EXPORT Geometry2D(const int coordim);
            SPATIAL_DOMAINS_EXPORT ~Geometry2D();


            StdRegions::StdExpansion2DSharedPtr GetXmap(const int i)
            {
                return m_xmap[i];
            }

            // Wrappers around virtual Functions (for the FaceComponent classes)
            void AddElmtConnected(int gvo_id, int locid)
            {
                v_AddElmtConnected(gvo_id,locid);
            }

            int  NumElmtConnected() const
            {
                return v_NumElmtConnected();
            }

            bool IsElmtConnected(int gvo_id, int locid) const
            {
                return v_IsElmtConnected(gvo_id,locid);
            }
            
            inline int GetFid() const
            {
                return v_GetFid();
            }

            inline int GetCoordDim() const 
            {
                return v_GetCoordDim();
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
            {
                return v_GetBasis(i,j);
            }

            inline const LibUtilities::BasisSharedPtr GetEdgeBasis(const int i, const int j)
            {
                return v_GetEdgeBasis(i,j);
            }

            inline Array<OneD,NekDouble> &UpdatePhys(const int i)
            {
                return v_UpdatePhys(i);
            }

            NekDouble GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                return v_GetCoord(i,Lcoord);
            }

            // Wrappers around virtual Functions (for the QuadGeom and TriGeom classes)
            inline void SetOwnData()
            {
                v_SetOwnData();
            }

            void FillGeom()
            {
                v_FillGeom();
            }            
             
            void GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                v_GetLocCoords(coords,Lcoords);
            }

            inline int GetEid(int i) const
            {
                return v_GetEid(i);
            }

            inline int GetVid(int i) const
            {
                return v_GetVid(i);
            }

            inline const Geometry2DSharedPtr GetFace(int i) const
            {
                return v_GetFace(i);
            }

            inline StdRegions::FaceOrientation GetFaceorient(const int i) const
            {
                return v_GetFaceorient(i);
            }
                
            inline const Geometry1DSharedPtr GetEdge(int i) const
            {
                return v_GetEdge(i);
            }
                
            inline const VertexComponentSharedPtr GetVertex(int i) const
            {
                return v_GetVertex(i);
            }

            inline StdRegions::EdgeOrientation GetCartesianEorient(const int i) const
            {
                return v_GetCartesianEorient(i);
            }

            int WhichEdge(SegGeomSharedPtr edge)
            {
                return v_WhichEdge(edge);
            }

            int WhichFace(Geometry2DSharedPtr face)
            {
                return v_WhichFace(face);
            }

            StdRegions::StdExpansion2DSharedPtr operator[](const int i) const
            {
                if((i>=0)&& (i<m_coordim))
                {
                    return m_xmap[i];
                }
                
                NEKERROR(ErrorUtil::efatal,
                         "Invalid Index used in [] operator");
                return m_xmap[0]; //should never be reached
            }
            

        protected:

            Array<OneD, StdRegions::StdExpansion2DSharedPtr> m_xmap;

        private:

            virtual void v_AddElmtConnected(int gvo_id, int locid)
            {  
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");            
            }

            virtual int v_NumElmtConnected() const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual bool v_IsElmtConnected(int gvo_id, int locid) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return false;
            }
            
            virtual int v_GetFid() const 
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual int v_GetCoordDim() const 
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                LibUtilities::BasisSharedPtr returnval;        
                return returnval;    
            }


            virtual const LibUtilities::BasisSharedPtr v_GetEdgeBasis(const int i, const int j)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                LibUtilities::BasisSharedPtr returnval;        
                return returnval;    
            }

            virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return NullNekDouble1DArray;
            }

            virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0.0;
            }

            virtual void v_SetOwnData()
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
            }

            virtual void v_FillGeom()
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
            }            
             
            virtual void v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
            }

            virtual int v_GetEid(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual int v_GetVid(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }
                
            virtual const Geometry1DSharedPtr v_GetEdge(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                SegGeomSharedPtr returnval;
                return returnval;
            }

            virtual const VertexComponentSharedPtr v_GetVertex(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                VertexComponentSharedPtr returnval;
                return returnval;
            }

            virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return StdRegions::eForwards;
            }

            virtual const Geometry2DSharedPtr v_GetFace(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                Geometry2DSharedPtr returnval;
                return returnval;                
            }

            virtual StdRegions::FaceOrientation v_GetFaceorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
            }

            virtual StdRegions::EdgeOrientation v_GetCartesianEorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return StdRegions::eForwards;
            }

            virtual int v_WhichEdge(SegGeomSharedPtr edge)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual int v_WhichFace(Geometry2DSharedPtr face)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }


            virtual int v_GetShapeDim() const
            {
                return 2;
            }
        };


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H

//
// $Log: Geometry2D.h,v $
// Revision 1.9  2008/11/17 09:00:11  ehan
// Added GetNumVerts and GetNumEdges
//
// Revision 1.8  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.7  2008/07/29 22:23:36  sherwin
// various mods for DG advection solver in Multiregions. Added virtual calls to Geometry, Geometry1D, 2D and 3D
//
// Revision 1.6  2008/06/16 22:37:32  ehan
// Added inline function GetFace(..) and GetFaceorient(..).
//
// Revision 1.5  2008/05/10 18:27:33  sherwin
// Modifications necessary for QuadExp Unified DG Solver
//
// Revision 1.4  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.3  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.2  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.16  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.15  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.14  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.13  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//


