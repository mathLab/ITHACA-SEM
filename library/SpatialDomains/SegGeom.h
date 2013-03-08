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
//  Description: Segment geometry information
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_SEGGEOM_H
#define NEKTAR_SPATIALDOMAINS_SEGGEOM_H


#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class SegGeom;
        typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< SegGeomSharedPtr >::iterator SegGeomVectorIter;
        typedef std::map<int, SegGeomSharedPtr> SegGeomMap;


        class SegGeom: public Geometry1D
        {
            public:
                SPATIAL_DOMAINS_EXPORT SegGeom();

                SPATIAL_DOMAINS_EXPORT SegGeom(int id, const int coordim);

                SPATIAL_DOMAINS_EXPORT SegGeom(
                        int id,
                        const int coordim,
                        const VertexComponentSharedPtr vertex[]);

                SPATIAL_DOMAINS_EXPORT SegGeom(
                        int id,
                        const int coordim,
                        const VertexComponentSharedPtr vertex[],
                        const CurveSharedPtr &curve);

                SPATIAL_DOMAINS_EXPORT SegGeom(
                        const int id,
                        const VertexComponentSharedPtr& vert1,
                        const VertexComponentSharedPtr& vert2);

                SPATIAL_DOMAINS_EXPORT SegGeom(const SegGeom &in);

                SPATIAL_DOMAINS_EXPORT ~SegGeom();

                SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(
                        const int i,
                        const Array<OneD, const NekDouble> &Lcoord);

                SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation
                        GetEdgeOrientation(const SegGeom& edge1,
                                           const SegGeom& edge2);

                SPATIAL_DOMAINS_EXPORT StdRegions::StdExpansion1DSharedPtr
                        operator[](const int i) const;


                SPATIAL_DOMAINS_EXPORT static const int           kNverts = 2;
                SPATIAL_DOMAINS_EXPORT static const int           kNedges = 1;

            protected:
                int                                               m_eid;
                std::list<CompToElmt>                             m_elmtMap;
                Array<OneD, StdRegions::StdExpansion1DSharedPtr>  m_xmap;
                SpatialDomains::VertexComponentSharedPtr          m_verts[kNverts];
                StdRegions::Orientation                           m_porient[kNverts];


                SPATIAL_DOMAINS_EXPORT virtual int v_GetVid(int i) const;

                SPATIAL_DOMAINS_EXPORT virtual VertexComponentSharedPtr
                        v_GetVertex(const int i) const;

                SPATIAL_DOMAINS_EXPORT virtual int v_GetEid() const;

                SPATIAL_DOMAINS_EXPORT virtual const LibUtilities::BasisSharedPtr
                        v_GetBasis(const int i, const int j);

                SPATIAL_DOMAINS_EXPORT virtual const
                        StdRegions::StdExpansion1DSharedPtr& v_GetXmap(const int i);

                SPATIAL_DOMAINS_EXPORT virtual void v_SetOwnData();

                SPATIAL_DOMAINS_EXPORT virtual void
                        v_AddElmtConnected(int gvo_id, int locid);

                SPATIAL_DOMAINS_EXPORT virtual int v_NumElmtConnected() const;

                SPATIAL_DOMAINS_EXPORT virtual bool
                        v_IsElmtConnected(int gvoId, int locId) const;

                SPATIAL_DOMAINS_EXPORT virtual Array<OneD, NekDouble>&
                        v_UpdatePhys(const int i);

                SPATIAL_DOMAINS_EXPORT virtual LibUtilities::ShapeType
                        v_DetShapeType() const;

                SPATIAL_DOMAINS_EXPORT virtual void v_GetLocCoords(
                        const Array<OneD, const NekDouble>& coords,
                              Array<OneD,NekDouble>& Lcoords);

                SPATIAL_DOMAINS_EXPORT virtual void v_GenGeomFactors(
                        const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

                SPATIAL_DOMAINS_EXPORT virtual StdRegions::Orientation
                        v_GetPorient(const int i) const;

                SPATIAL_DOMAINS_EXPORT virtual void v_FillGeom ();

                SPATIAL_DOMAINS_EXPORT virtual NekDouble v_GetCoord(
                        const int i,
                        const Array<OneD,const NekDouble> &Lcoord);

                SPATIAL_DOMAINS_EXPORT virtual void v_WriteToFile(
                              std::ofstream &outfile,
                        const int dumpVar);

                SPATIAL_DOMAINS_EXPORT virtual int  v_GetNumVerts() const;

                SPATIAL_DOMAINS_EXPORT virtual int  v_GetNumEdges() const;

                SPATIAL_DOMAINS_EXPORT virtual bool v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                        NekDouble tol = 0.0);

            private:
                /// Boolean indicating whether object owns the data
                bool                            m_ownData;
        };


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_SEGGEOM_H

