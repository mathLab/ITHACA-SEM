////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry3D.h
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
//  Description: 3D geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/GeomFactors3D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry3D : public Geometry
        {
        public:
            SPATIAL_DOMAINS_EXPORT Geometry3D();
            SPATIAL_DOMAINS_EXPORT Geometry3D(const int coordim);
            SPATIAL_DOMAINS_EXPORT virtual ~Geometry3D();
            
            //---------------------------------------
            // 3D Geometry Methods
            //---------------------------------------
            SPATIAL_DOMAINS_EXPORT void FillGeom();
            SPATIAL_DOMAINS_EXPORT void GetLocCoords(
                const Array<OneD, const NekDouble> &coords,
                      Array<OneD,       NekDouble> &Lcoords);
            SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(
                const int i, const Array<OneD, const NekDouble> &Lcoord);

            //---------------------------------------
            // Helper functions
            //---------------------------------------
            SPATIAL_DOMAINS_EXPORT int GetEid(int i) const;
            SPATIAL_DOMAINS_EXPORT Geometry2DSharedPtr 
                GetFace(int i);
            SPATIAL_DOMAINS_EXPORT StdRegions::FaceOrientation 
                GetFaceorient(const int i) const;
            SPATIAL_DOMAINS_EXPORT int GetFid(int i) const;
            SPATIAL_DOMAINS_EXPORT Array<OneD,NekDouble> 
                &UpdatePhys(const int i);
            SPATIAL_DOMAINS_EXPORT const LibUtilities::BasisSharedPtr 
                GetBasis(const int i, const int j);

            //---------------------------------------
            // Element connection functions
            //---------------------------------------
            SPATIAL_DOMAINS_EXPORT void SetOwnData();
            
        protected:
            VertexComponentVector                            m_verts;
            SegGeomVector                                    m_edges;
            Geometry2DVector                                 m_faces;
            std::vector<StdRegions::EdgeOrientation>         m_eorient;
            std::vector<StdRegions::FaceOrientation>         m_forient;
            Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;
            std::list<CompToElmt>                            m_elmtmap;
            bool                                             m_owndata;
            int                                              m_eid;
            bool                                             m_ownverts;

            //---------------------------------------
            // 3D Geometry Methods
            //---------------------------------------
            virtual void      v_FillGeom();
            virtual void      v_GetLocCoords(
                const Array<OneD, const NekDouble> &coords,
                      Array<OneD,       NekDouble> &Lcoords);
            virtual void      v_GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);
            virtual NekDouble v_GetCoord(
                const int i, const Array<OneD, const NekDouble> &Lcoord);

            //---------------------------------------
            // Helper functions
            //---------------------------------------
            StdRegions::StdExpansion3DSharedPtr v_GetXmap(const int i);
            virtual int                         v_GetShapeDim() const;
            virtual int                         v_GetVid(int i) const;
            virtual const SegGeomSharedPtr      v_GetEdge(int i) const;
            virtual StdRegions::EdgeOrientation v_GetEorient(const int i) const;
            virtual int                         v_GetEid(int i) const;
            virtual const Geometry2DSharedPtr   v_GetFace(int i) const;
            virtual StdRegions::FaceOrientation v_GetFaceorient(const int i) const;
            virtual int                         v_GetFid(int i) const;
            virtual int                         v_GetEid() const;
            virtual Array<OneD,NekDouble>      &v_UpdatePhys(const int i);
            virtual int                         v_WhichEdge(SegGeomSharedPtr edge);
            virtual int                         v_WhichFace(Geometry2DSharedPtr face);
            const LibUtilities::BasisSharedPtr  v_GetBasis(const int i, 
                                                           const int j);

            //---------------------------------------
            // Element connection functions
            //---------------------------------------
            virtual void v_AddElmtConnected(int gvo_id, int locid);
            virtual bool v_IsElmtConnected (int gvo_id, int locid) const;
            virtual int  v_NumElmtConnected() const;
            virtual void v_SetOwnData();
        };

        typedef boost::shared_ptr<Geometry3D> Geometry3DSharedPtr;
        typedef std::vector<Geometry3DSharedPtr> Geometry3DVector;
        typedef std::vector<Geometry3DSharedPtr>::iterator Geometry3DVectorIter;
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
