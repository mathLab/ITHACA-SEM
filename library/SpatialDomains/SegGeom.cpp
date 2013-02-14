////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.cpp
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

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors1D.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

namespace Nektar
{
    namespace SpatialDomains
    {
        SegGeom::SegGeom()
        {
            m_geomShapeType = eSegment;
        }

        SegGeom::SegGeom(int id, const int coordim):
            Geometry1D(coordim),
            m_xmap(coordim)
        {
            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,
                                              LibUtilities::eGaussLobattoLegendre
                                           )
                                          );
            m_geomShapeType = eSegment;
            m_eid = id;

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }
        }

        SegGeom::SegGeom(
                int id,
                const int coordim,
                const VertexComponentSharedPtr vertex[]):
            Geometry1D(coordim)
        {
            m_geomShapeType = eSegment;
            m_eid   = id;
            m_state = eNotFilled;

            if (coordim > 0)
            {
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                               LibUtilities::PointsKey(3,
                                                  LibUtilities::eGaussLobattoLegendre
                                               )
                                              );

                m_xmap = Array<OneD, StdRegions::StdExpansion1DSharedPtr>(m_coordim);

                for(int i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                }
            }

            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }

        SegGeom::SegGeom(
                int id,
                const int coordim,
                const VertexComponentSharedPtr vertex[],
                const CurveSharedPtr& curve):
            Geometry1D(coordim)
        {
            m_geomShapeType = eSegment;
            m_eid = id;
            m_state = eNotFilled;

            if (coordim > 0)
            {
                int npts = curve->m_points.size();
                LibUtilities::PointsKey pkey(npts+1,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, npts, pkey);

                m_xmap = Array<OneD, StdRegions::StdExpansion1DSharedPtr>(m_coordim);

                Array<OneD,NekDouble> tmp(npts);


                if(vertex[0]->dist(*(curve->m_points[0])) > NekConstants::kVertexTheSameDouble)
                { 
                    cout<<"edge="<<id<<endl;                    
                    std::string err = "Vertex 0 is separated from first point by more than ";
                    std::stringstream strstrm;
                    strstrm << NekConstants::kVertexTheSameDouble;
                    err += strstrm.str();
                    NEKERROR(ErrorUtil::ewarning, err.c_str());
                }

                if(vertex[1]->dist(*(curve->m_points[npts-1])) > NekConstants::kVertexTheSameDouble)
                {
                    cout<<"edge="<<id<<endl;      
                    std::string err = "Vertex 1 is separated from last point by more than ";
                    std::stringstream strstrm;
                    strstrm << NekConstants::kVertexTheSameDouble;
                    err += strstrm.str();
                    NEKERROR(ErrorUtil::ewarning, err.c_str());
                }


                for(int i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);

                    // Load up coordinate values into tmp
                    for(int j = 0; j < npts; ++j)
                    {
                        tmp[j] = (curve->m_points[j]->GetPtr())[i];
                    }

                    // Interpolate to GLL points
                    DNekMatSharedPtr I0;

                    LibUtilities::PointsKey fkey(npts,curve->m_ptype);
                    I0 = LibUtilities::PointsManager()[fkey]->GetI(pkey);

                    NekVector<NekDouble> in(npts,tmp,eWrapper);
                    NekVector<NekDouble>       out(npts+1,m_xmap[i]->UpdatePhys(),eWrapper);
                    out  = (*I0)*in;

                    m_xmap[i]->FwdTrans(m_xmap[i]->GetPhys(),m_xmap[i]->UpdateCoeffs());
                }
            }

            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }


        SegGeom::SegGeom(
                const int id,
                const VertexComponentSharedPtr& vert1,
                const VertexComponentSharedPtr& vert2):
            Geometry1D(vert1->GetCoordim()), m_xmap(vert1->GetCoordim())
        {
            m_geomShapeType = eSegment;

            m_verts[0] = vert1;
            m_verts[1] = vert2;

            m_state = eNotFilled;

            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,
                                              LibUtilities::eGaussLobattoLegendre
                                           )
                                          );
            m_eid = id;

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }
        }

        SegGeom::SegGeom(const SegGeom &in)
        {
            // From Geometry class
            m_geomShapeType = in.m_geomShapeType;

            // info from EdgeComponent class
            m_eid     = in.m_eid;
            std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtMap.begin(); def != in.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }
            m_xmap = in.m_xmap;

            // info from SegGeom class
            m_coordim  = in.m_coordim;
            m_verts[0] = in.m_verts[0];
            m_verts[1] = in.m_verts[1];

            m_state = in.m_state;
        }

        SegGeom::~SegGeom()
        {
        }




        /** given local collapsed coordinate Lcoord return the value of
            physical coordinate in direction i **/
        NekDouble SegGeom::v_GetCoord(
                const int i,
                const Array<OneD, const NekDouble>& Lcoord)
        {

            ASSERTL1(m_state == ePtsFilled, "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        void SegGeom::v_AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtMap.push_back(ee);
        }

        int SegGeom::v_NumElmtConnected() const
        {
            return int(m_elmtMap.size());
        }

        bool SegGeom::v_IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtMap.begin(),m_elmtMap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtMap.end())
            {
                return(true);
            }

            return(false);
        }



        ///  \brief Get the orientation of edge1.
        ///
        ///  If edge1 is connected to edge2 in the same direction as
        ///  the points comprising edge1 then it is forward, otherwise
        ///  it is backward.
        ///
        ///  For example, assume edge1 is comprised of points 1 and 2,
        ///  and edge2 is comprised of points 2 and 3, then edge1 is
        ///  forward.
        ///
        ///  If edge1 is comprised of points 2 and 1 and edge2 is
        ///  comprised of points 3 and 2, then edge1 is backward.
        ///
        ///  Since both edges are passed, it does
        ///  not need any information from the EdgeComponent instance.
        StdRegions::Orientation SegGeom::GetEdgeOrientation(
                const SegGeom& edge1,
                const SegGeom& edge2)
        {
            StdRegions::Orientation returnval = StdRegions::eForwards;

            /// Backward direction.  Vertex 0 is connected to edge 2.
            if ((*edge1.GetVertex(0) == *edge2.GetVertex(0)) ||
                (*edge1.GetVertex(0) == *edge2.GetVertex(1)))
            {
                returnval = StdRegions::eBackwards;
            }
            // Not forward either, then we have a problem.
            else if ((*edge1.GetVertex(1) != *edge2.GetVertex(0)) &&
                (*edge1.GetVertex(1) != *edge2.GetVertex(1)))
            {
                std::ostringstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << edge1.GetEid() << ", " << edge2.GetEid();
                ASSERTL0(false, errstrm.str());
            }

            return returnval;
        }

        ///  Set up GeoFac for this geometry using Coord quadrature distribution
        void SegGeom::v_GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr>& tbasis)
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                SpatialDomains::GeomType gType = eRegular;

                SegGeom::v_FillGeom();

                if(m_xmap[0]->GetBasisNumModes(0)!=2)
                {
                    gType = eDeformed;
                }

                m_geomFactors = MemoryManager<GeomFactors1D>::AllocateSharedPtr(gType,
                                                             m_coordim, m_xmap, tbasis);

                m_geomFactorsState = ePtsFilled;
            }
        }


        ///  \brief put all quadrature information into edge structure and
        ///    backward transform
        void SegGeom::v_FillGeom()
        {
            if(m_state != ePtsFilled)
            {
                int i;

                for(i = 0; i < m_coordim; ++i){
                    m_xmap[i]->SetCoeff(0,(*m_verts[0])[i]);
                    m_xmap[i]->SetCoeff(1,(*m_verts[1])[i]);
                    m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(),m_xmap[i]->UpdatePhys());
                }
                m_state = ePtsFilled;
            }
        }

        void SegGeom::v_GetLocCoords(
                const Array<OneD, const NekDouble>& coords,
                      Array<OneD,NekDouble>& Lcoords)
        {
            int i;

            SegGeom::v_FillGeom();

            // calculate local coordinate for coord
            if(GetGtype() == eRegular)
            {
                Array<OneD, const NekDouble> pts;
                NekDouble len = 0.0;
                NekDouble xi  = 0.0;
                int nq;

                // get points;
                //find end points
                for(i = 0; i < m_coordim; ++i)
                {
                    nq   = m_xmap[i]->GetNumPoints(0);
                    pts  = m_xmap[i]->GetPhys();
                    len  += (pts[nq-1]-pts[0])*(pts[nq-1]-pts[0]);
                    xi   += (coords[i]-pts[0])*(coords[i]-pts[0]);
                }

                len = sqrt(len);
                xi  = sqrt(xi);

                Lcoords[0] =  2*xi/len-1.0;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "inverse mapping must be set up to use this call");
            }
        }

        void SegGeom::v_WriteToFile(std::ofstream &outfile, const int dumpVar)
        {

            int i,j;
            int  nquad = m_xmap[0]->GetNumPoints(0);
            NekDouble *coords[3];

            SegGeom::v_FillGeom();

            for(i = 0; i < m_coordim; ++i)
            {
                coords[i] = &(m_xmap[i]->UpdatePhys())[0];
            }

            if(dumpVar)
            {
                outfile << "Variables = x";
                if(m_coordim == 2)
                {
                    outfile << ", y";
                }
                else if (m_coordim == 3)
                {
                    outfile << ", y, z";
                }
                outfile << std::endl;
            }

            outfile << "Zone, I=" << nquad << ", F=Point\n";
            for(i = 0; i < nquad; ++i)
            {
                for(j = 0; j < m_coordim; ++j)
                {
                    outfile << coords[j][i] << " ";
                }
                outfile << std::endl;
            }
        }

        bool SegGeom::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                NekDouble tol)
        {
            Array<OneD,NekDouble> stdCoord(GetCoordim(),0.0);
            GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -(1+tol) && stdCoord[0] <= 1+tol)
            {
                return true;
            }
            return false;
        }

        int SegGeom::v_GetVid(int i) const
        {
            ASSERTL2((i >=0) && (i <= 1),"Verted id must be between 0 and 1");
            return m_verts[i]->GetVid();
        }

        VertexComponentSharedPtr SegGeom::v_GetVertex(const int i) const
        {
            VertexComponentSharedPtr returnval;

            if (i >= 0 && i < kNverts)
            {
                returnval = m_verts[i];
            }

            return returnval;
        }

        int SegGeom::v_GetEid() const
        {
            return m_eid;
        }

        const LibUtilities::BasisSharedPtr SegGeom::v_GetBasis(
                const int i,
                const int j)
        {
            return m_xmap[i]->GetBasis(j);
        }

        StdRegions::StdExpansion1DSharedPtr SegGeom::operator[](const int i) const
        {
            if((i>=0)&& (i<m_coordim))
            {
                return m_xmap[i];
            }

            NEKERROR(ErrorUtil::efatal, "Invalid Index used in [] operator");
            return m_xmap[0]; //should never be reached
        }

        const StdRegions::StdExpansion1DSharedPtr& SegGeom::v_GetXmap(const int i)
        {
            return m_xmap[i];
        }

        void SegGeom::v_SetOwnData()
        {
            m_ownData = true;
        }

        Nektar::Array<OneD, NekDouble>& SegGeom::v_UpdatePhys(const int i)
        {
            return m_xmap[i]->UpdatePhys();
        }

        StdRegions::Orientation SegGeom::v_GetPorient(const int i) const
        {
            //cout << "StdRegions::PointOrientation GetPorient"<<endl;
            ASSERTL2((i >=0) && (i <= 1),"Point id must be between 0 and 1");

            if (i%2==0)
            {
                return StdRegions::eBwd;
            }
            else 
            {
                return StdRegions::eFwd;
            }

            //return m_porient[i];
        }


        StdRegions::ExpansionType SegGeom::v_DetExpansionType() const
        {
            return StdRegions::eSegment;
        }

        int SegGeom::v_GetNumVerts() const
        {
            return kNverts;
        }

        int SegGeom::v_GetNumEdges() const
        {
            return kNedges;
        }


    }; //end of namespace
}; //end of namespace

