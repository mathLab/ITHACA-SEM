////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry3D.h
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry2D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry3D;
        
        // shorthand for boost pointer
        typedef boost::shared_ptr<Geometry3D> Geometry3DSharedPtr;
        typedef std::vector< Geometry3DSharedPtr > Geometry3DVector;
        typedef std::vector< Geometry3DSharedPtr >::iterator Geometry3DVectorIter;

    class Geometry3D: public Geometry
        {
        public:
            Geometry3D();
            Geometry3D(const int coordim);
            ~Geometry3D();

            StdRegions::StdExpansion3DSharedPtr GetXmap(const int i)
            {
                return m_xmap[i];
            }

            Geometry1DSharedPtr GetEdge(const int j)
            {
                return v_GetEdge(j);
            }

            inline const Geometry1DSharedPtr GetEdge(int i) const
            {
                return v_GetEdge(i);
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

            inline int GetVid(int i) const
            {
                return v_GetVid(i);
            }

            inline int GetEid(int i) const
            {
                return v_GetEid(i);
            }
            
            inline int GetFid(int i) const
            {
                return v_GetFid(i);
            }

            void FillGeom()
            {
                v_FillGeom();
            }            
                         
            void GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                v_GetLocCoords(coords,Lcoords);
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

            inline void SetOwnData()
            {
                v_SetOwnData();
            }

            inline const Geometry2DSharedPtr GetFace(int i) const
            {
                return v_GetFace(i);
            }

            inline StdRegions::FaceOrientation GetFaceorient(const int i) const
            {
                return v_GetFaceorient(i);
            }

            inline StdRegions::EdgeOrientation GetEorient(const int i) const
            {
                return v_GetEorient(i);
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
            

        protected:

           Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;

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

            virtual int v_GetFid(int i) const 
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

            virtual int v_GetVid(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual const Geometry1DSharedPtr v_GetEdge(int i) const
            {
                NEKERROR(ErrorUtil::efatal,"This function is only valid for shape type geometries");
                SegGeomSharedPtr returnval;
                return returnval;
            }

            virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for shape type geometries");
                return StdRegions::eForwards;
            }

            virtual const Geometry2DSharedPtr v_GetFace(int i) const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for shape type geometries");
                Geometry2DSharedPtr returnval;
                return returnval;                
            }

            virtual StdRegions::FaceOrientation v_GetFaceorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for shape type geometries");
                return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
            }

            virtual StdRegions::EdgeOrientation v_GetCartesianEorient(const int i) const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for shape type geometries");
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
            
        };

    }; //end of namespace
}; //end of namespace


#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H

//
// $Log: Geometry3D.h,v $
// Revision 1.7  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.6  2008/08/26 02:27:19  ehan
// Added various virtual functions.
//
// Revision 1.5  2008/07/29 22:23:36  sherwin
// various mods for DG advection solver in Multiregions. Added virtual calls to Geometry, Geometry1D, 2D and 3D
//
// Revision 1.4  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.3  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.2  2008/01/31 11:02:55  ehan
// Added constructor and destructor.
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.14  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.12  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.11  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

