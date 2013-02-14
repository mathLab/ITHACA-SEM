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
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/SegGeom.h>

#include <SpatialDomains/GeomFactors2D.h>
#include <StdRegions/StdTriExp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class  TriGeom;
        class  SegGeom;
        struct Curve;

        typedef boost::shared_ptr<Curve> CurveSharedPtr;
        typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef boost::shared_ptr<TriGeom> TriGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< TriGeomSharedPtr > TriGeomVector;
        typedef std::vector< TriGeomSharedPtr >::iterator TriGeomVectorIter;
        typedef std::map<int, TriGeomSharedPtr> TriGeomMap;
        typedef std::map<int, TriGeomSharedPtr>::iterator TriGeomMapIter;

        class TriGeom: public Geometry2D
        {
            public:
                SPATIAL_DOMAINS_EXPORT TriGeom();

                SPATIAL_DOMAINS_EXPORT TriGeom(int id, const int coordim);

                SPATIAL_DOMAINS_EXPORT TriGeom(
                        const int id,
                        const VertexComponentSharedPtr verts[],
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[]);

                SPATIAL_DOMAINS_EXPORT TriGeom(
                        const int id,
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[]);

                SPATIAL_DOMAINS_EXPORT TriGeom(
                        const int id,
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[],
                        const CurveSharedPtr &curve);

                SPATIAL_DOMAINS_EXPORT TriGeom(const TriGeom &in);

                SPATIAL_DOMAINS_EXPORT ~TriGeom();

                SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(
                        const int i,
                        const Array<OneD, const NekDouble> &Lcoord);

                /// Get the orientation of face1.
                SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation
                            GetFaceOrientation(const TriGeom &face1,
                                               const TriGeom &face2);

                SPATIAL_DOMAINS_EXPORT static const int kNedges = 3;
                SPATIAL_DOMAINS_EXPORT static const int kNverts = 3;

            protected:
                VertexComponentVector           m_verts;
                SegGeomVector                   m_edges;
                StdRegions::Orientation         m_eorient [kNedges];
                int				                m_fid;
                bool				            m_ownVerts;
                std::list<CompToElmt>		    m_elmtMap;

                SPATIAL_DOMAINS_EXPORT virtual void v_AddElmtConnected(
                        int gvo_id,
                        int locid);

                SPATIAL_DOMAINS_EXPORT virtual int v_NumElmtConnected() const;

                SPATIAL_DOMAINS_EXPORT virtual bool v_IsElmtConnected(
                        int gvo_id,
                        int locid) const;

                SPATIAL_DOMAINS_EXPORT virtual int v_GetFid() const;

                SPATIAL_DOMAINS_EXPORT virtual int v_GetCoordim() const;

                SPATIAL_DOMAINS_EXPORT virtual const LibUtilities::BasisSharedPtr
                            v_GetBasis(const int i, const int j);

                SPATIAL_DOMAINS_EXPORT virtual const LibUtilities::BasisSharedPtr
                            v_GetEdgeBasis(const int i, const int j);

                SPATIAL_DOMAINS_EXPORT virtual Array<OneD,NekDouble> &
                            v_UpdatePhys(const int i);

                SPATIAL_DOMAINS_EXPORT virtual NekDouble v_GetCoord(
                        const int i,
                        const Array<OneD,const NekDouble> &Lcoord);

                SPATIAL_DOMAINS_EXPORT virtual void v_GenGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

                SPATIAL_DOMAINS_EXPORT virtual void v_SetOwnData();

                /// Put all quadrature information into edge structure
                SPATIAL_DOMAINS_EXPORT virtual void v_FillGeom();

                SPATIAL_DOMAINS_EXPORT virtual void v_GetLocCoords(
                        const Array<OneD,const NekDouble> &coords,
                              Array<OneD,      NekDouble> &Lcoords);

                SPATIAL_DOMAINS_EXPORT virtual int v_GetEid(int i) const;

                SPATIAL_DOMAINS_EXPORT virtual int v_GetVid(int i) const;

                SPATIAL_DOMAINS_EXPORT virtual const VertexComponentSharedPtr
                            v_GetVertex(int i) const;

                SPATIAL_DOMAINS_EXPORT virtual const Geometry1DSharedPtr
                            v_GetEdge(int i) const;

                SPATIAL_DOMAINS_EXPORT virtual StdRegions::Orientation
                            v_GetEorient(const int i) const;

                SPATIAL_DOMAINS_EXPORT virtual StdRegions::Orientation
                            v_GetCartesianEorient(const int i) const;

                /// Return the edge number of the given edge
                SPATIAL_DOMAINS_EXPORT virtual int v_WhichEdge(
                        SegGeomSharedPtr edge);

                SPATIAL_DOMAINS_EXPORT virtual int v_GetNumVerts() const;

                SPATIAL_DOMAINS_EXPORT virtual int v_GetNumEdges() const;

                SPATIAL_DOMAINS_EXPORT virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble> &gloCoord,
                              NekDouble tol = 0.0);

            private:
                bool                            m_ownData;
        };

    }; //end of namespace SpatialDomains
}; //end of namespace Nektar

#endif
