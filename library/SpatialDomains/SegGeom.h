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

#include <SpatialDomains/GeomFactors1D.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/MeshComponents.h>

#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/Basis.h>

#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {

        class SegGeom: public Geometry1D
        {
            public:
                SPATIAL_DOMAINS_EXPORT SegGeom();

                SPATIAL_DOMAINS_EXPORT SegGeom(int id, const int coordim);

                SPATIAL_DOMAINS_EXPORT SegGeom(int id, const int coordim, const VertexComponentSharedPtr vertex[]);

                SPATIAL_DOMAINS_EXPORT SegGeom(int id, const int coordim, const VertexComponentSharedPtr vertex[], const CurveSharedPtr &curve);

                SPATIAL_DOMAINS_EXPORT SegGeom(const int id, const VertexComponentSharedPtr vert1, const VertexComponentSharedPtr  vert2);

                SPATIAL_DOMAINS_EXPORT SegGeom(const SegGeom &in);

                SPATIAL_DOMAINS_EXPORT ~SegGeom();

                SPATIAL_DOMAINS_EXPORT void AddElmtConnected(int gvoId, int locId);
                SPATIAL_DOMAINS_EXPORT int NumElmtConnected() const;
                SPATIAL_DOMAINS_EXPORT bool IsElmtConnected(int gvoId, int locId) const;

                inline int GetEid() const
                {
                    return m_eid;
                }

                inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
                {
                    return m_xmap[i]->GetBasis(j);
                }

                inline const StdRegions::StdExpansion1DSharedPtr &GetXmap(const int i)
                {
                    return m_xmap[i];
                }

                inline Array<OneD, NekDouble> &UpdatePhys(const int i)
                {
                    return m_xmap[i]->UpdatePhys();
                }

                inline VertexComponentSharedPtr GetVertex(const int i) const
                {
                    VertexComponentSharedPtr returnval;

                    if (i >= 0 && i < kNverts)
                    {
                        returnval = m_verts[i];
                    }

                    return returnval;
                }


                StdRegions::StdExpansion1DSharedPtr operator[](const int i) const
                {
                    if((i>=0)&& (i<m_coordim))
                    {
                        return m_xmap[i];
                    }

                    NEKERROR(ErrorUtil::efatal,
                        "Invalid Index used in [] operator");
                    return m_xmap[0]; //should never be reached
                }

                SPATIAL_DOMAINS_EXPORT NekDouble GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord);

                /// \brief Get the orientation of edge1.
                ///
                /// Since both edges are passed, it does
                /// not need any information from the EdgeComponent instance.
                SPATIAL_DOMAINS_EXPORT static StdRegions::EdgeOrientation GetEdgeOrientation(const SegGeom &edge1,
                    const SegGeom &edge2);

                inline int GetVid(int i) const
                {
                    ASSERTL2((i >=0) && (i <= 1),"Verted id must be between 0 and 1");
                    return m_verts[i]->GetVid();
                }

                inline void SetOwnData()
                {
                    m_ownData = true;
                }

                SPATIAL_DOMAINS_EXPORT void    FillGeom ();

                StdRegions::ExpansionType DetExpansionType() const
                {
                    return StdRegions::eSegment;
                }

                SPATIAL_DOMAINS_EXPORT void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords);


                SPATIAL_DOMAINS_EXPORT void WriteToFile(std::ofstream &outfile, const int dumpVar);

            protected:
                int                             m_eid;
                std::list<CompToElmt>           m_elmtMap;
                Array<OneD, StdRegions::StdExpansion1DSharedPtr> m_xmap;

                static const int                kNverts = 2;
                static const int                kNedges = 1;
                SpatialDomains::VertexComponentSharedPtr m_verts[kNverts];

                void GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);

            private:
                /// Boolean indicating whether object owns the data
                bool                            m_ownData;

                virtual void v_AddElmtConnected(int gvo_id, int locid)
                {
                    AddElmtConnected(gvo_id, locid);
                }

                virtual int v_NumElmtConnected() const
                {
                    return NumElmtConnected();
                }

                virtual bool v_IsElmtConnected(int gvo_id, int locid) const
                {
                    return IsElmtConnected(gvo_id, locid);
                }

                virtual int v_GetEid() const
                {
                    return GetEid();
                }

                virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
                {
                    return GetBasis(i,j);
                }

                virtual const StdRegions::StdExpansion1DSharedPtr &v_GetXmap(const int i)
                {
                    return GetXmap(i);
                }

                virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
                {
                    return UpdatePhys(i);
                }

                virtual VertexComponentSharedPtr v_GetVertex(const int i) const
                {
                    return GetVertex(i);
                }

                virtual void v_GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
                {
                    GenGeomFactors(tbasis);
                }

                virtual void v_SetOwnData()
                {
                    SetOwnData();
                }

                virtual int v_GetVid(int i) const
                {
                    return GetVid(i);
                }

                virtual void v_FillGeom()
                {
                    FillGeom();
                }

                virtual StdRegions::ExpansionType v_DetExpansionType() const
                {
                    return DetExpansionType();
                }

                virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
                {
                    return GetCoord(i,Lcoord);
                }

                virtual void v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
                {
                    GetLocCoords(coords,Lcoords);
                }


                virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
                {
                    WriteToFile(outfile, dumpVar);
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
                                             const Array<OneD, const NekDouble> &gloCoord, NekDouble tol = 0.0);
        };

        // shorthand for boost pointer
        typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< SegGeomSharedPtr >::iterator SegGeomVectorIter;
        typedef std::map<int, SegGeomSharedPtr> SegGeomMap;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_SEGGEOM_H

//
// $Log: SegGeom.h,v $
// Revision 1.25  2009/12/15 18:09:02  cantwell
// Split GeomFactors into 1D, 2D and 3D
// Added generation of tangential basis into GeomFactors
// Updated ADR2DManifold solver to use GeomFactors for tangents
// Added <GEOMINFO> XML session section support in MeshGraph
// Fixed const-correctness in VmathArray
// Cleaned up LocalRegions code to generate GeomFactors
// Removed GenSegExp
// Temporary fix to SubStructuredGraph
// Documentation for GlobalLinSys and GlobalMatrix classes
//
// Revision 1.24  2009/03/04 14:17:38  pvos
// Removed all methods that take and Expansion as argument
//
// Revision 1.23  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.22  2008/12/17 12:29:56  pvos
// Fixed bug
//
// Revision 1.21  2008/11/17 08:58:53  ehan
// Added GetNumVerts and GetNumEdges
//
// Revision 1.20  2008/09/09 14:22:39  sherwin
// Added curved segment constructor methods
//
// Revision 1.19  2008/06/12 23:27:57  delisi
// Removed MeshGraph.h include from SegGeom.h, to get rid of circular includes. Now can use typedefs from SegGeom.h instead of repeating it in MeshGraph.h.
//
// Revision 1.18  2008/06/11 21:34:42  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.17  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.16  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
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
