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
#include <SpatialDomains/GeomFactors.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

namespace Nektar
{
    namespace SpatialDomains
    {
        SegGeom::SegGeom()
        {
            m_shapeType = LibUtilities::eSegment;
        }

        SegGeom::SegGeom(int id, const int coordim):
            Geometry1D(coordim)
        {
            const LibUtilities::BasisKey B(
                LibUtilities::eModified_A, 2,
                LibUtilities::PointsKey(3, LibUtilities::eGaussLobattoLegendre));

            m_shapeType = LibUtilities::eSegment;
            m_eid = id;
            m_globalID = id;
            m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        SegGeom::SegGeom(
                int id,
                const int coordim,
                const PointGeomSharedPtr vertex[]):
            Geometry1D(coordim)
        {
            m_shapeType = LibUtilities::eSegment;
            m_eid   = id;
            m_globalID = id;
            m_state = eNotFilled;

            if (coordim > 0)
            {
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                               LibUtilities::PointsKey(3,
                                                  LibUtilities::eGaussLobattoLegendre
                                               )
                                              );

                m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                SetUpCoeffs(m_xmap->GetNcoeffs());
            }

            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }

        SegGeom::SegGeom(
                int id,
                const int coordim,
                const PointGeomSharedPtr vertex[],
                const CurveSharedPtr& curve):
            Geometry1D(coordim)
        {
            m_shapeType = LibUtilities::eSegment;
            m_eid = id;
            m_globalID = id; 
            m_state = eNotFilled;

            if (coordim > 0)
            {
                int npts = curve->m_points.size();
                LibUtilities::PointsKey pkey(npts+1,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, npts, pkey);

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

                m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                SetUpCoeffs(m_xmap->GetNcoeffs());

                for(int i = 0; i < m_coordim; ++i)
                {
                    // Load up coordinate values into tmp
                    for(int j = 0; j < npts; ++j)
                    {
                        tmp[j] = (curve->m_points[j]->GetPtr())[i];
                    }

                    // Interpolate to GLL points
                    DNekMatSharedPtr I0;

                    LibUtilities::PointsKey fkey(npts,curve->m_ptype);
                    I0 = LibUtilities::PointsManager()[fkey]->GetI(pkey);

                    NekVector<NekDouble> in (npts, tmp, eWrapper);
                    NekVector<NekDouble> out(npts+1);
                    out = (*I0)*in;

                    m_xmap->FwdTrans(out.GetPtr(), m_coeffs[i]);
                }
            }

            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }


        SegGeom::SegGeom(
                const int id,
                const PointGeomSharedPtr& vert1,
                const PointGeomSharedPtr& vert2):
            Geometry1D(vert1->GetCoordim())
        {
            m_shapeType = LibUtilities::eSegment;

            m_verts[0] = vert1;
            m_verts[1] = vert2;

            m_state = eNotFilled;

            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,
                                              LibUtilities::eGaussLobattoLegendre
                                           )
                                          );
            m_eid = id;
            m_globalID = id;
            m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        SegGeom::SegGeom(const SegGeom &in)
        {
            // From Geometry class
            m_shapeType = in.m_shapeType;

            // info from EdgeComponent class
            m_eid     = in.m_eid;
            m_globalID = in.m_globalID;

            std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtMap.begin(); def != in.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }
            m_xmap = in.m_xmap;
            SetUpCoeffs(m_xmap->GetNcoeffs());

            // info from SegGeom class
            m_coordim  = in.m_coordim;
            m_verts[0] = in.m_verts[0];
            m_verts[1] = in.m_verts[1];

            m_state = in.m_state;
        }


        /** \brief Generate a one dimensional space segment geometry where
            the vert[0] has the same x value and vert[1] is set to
            vert[0] plus the length of the original segment  
            
        **/
        SegGeomSharedPtr SegGeom::GenerateOneSpaceDimGeom(void)
        {
            SegGeomSharedPtr returnval = MemoryManager<SegGeom>::AllocateSharedPtr();
            
            // info about numbering 
            returnval->m_eid       = m_eid;
            returnval->m_globalID  = m_globalID;
            returnval->m_elmtMap   = m_elmtMap; 
            

            // geometric information. 
            returnval->m_coordim = 1;
            NekDouble x0 = (*m_verts[0])[0];
            PointGeomSharedPtr vert0 = MemoryManager<PointGeom>::AllocateSharedPtr(1,m_verts[0]->GetVid(),x0,0.0,0.0);
            vert0->SetGlobalID(vert0->GetVid());
            returnval->m_verts[0] = vert0;
            
            // Get information to calculate length. 
            const Array<OneD, const LibUtilities::BasisSharedPtr> base =  m_xmap->GetBase();
            LibUtilities::PointsKeyVector v;
            v.push_back(base[0]->GetPointsKey());
            v_GenGeomFactors();

            const Array<OneD, const NekDouble> jac = m_geomFactors->GetJac(v);
            
            NekDouble len;
            if(jac.num_elements() == 1)
            {
                len = jac[0]*2.0;
            }
            else
            {
                Array<OneD, const NekDouble> w0 = base[0]->GetW();
                len = 0.0;

                for(int i = 0; i < jac.num_elements(); ++i)
                {
                    len += jac[i]*w0[i];
                }
            }
            // Set up second vertex. 
            PointGeomSharedPtr vert1 = MemoryManager<PointGeom>::AllocateSharedPtr(1,m_verts[1]->GetVid(),x0+len,0.0,0.0);
            vert1->SetGlobalID(vert1->GetVid());

            returnval->m_verts[1] = vert1;
            
            // at present just use previous m_xmap[0]; 
            returnval->m_xmap    = m_xmap;
            returnval->SetUpCoeffs(m_xmap->GetNcoeffs());
            returnval->m_state   = eNotFilled;

            return returnval;
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
            ASSERTL1(m_state == ePtsFilled,
                "Geometry is not in physical space");

            Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
            m_xmap->BwdTrans(m_coeffs[i], tmp);

            return m_xmap->PhysEvaluate(Lcoord, tmp);
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
        void SegGeom::v_GenGeomFactors()
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                SpatialDomains::GeomType gType = eRegular;

                SegGeom::v_FillGeom();

                if(m_xmap->GetBasisNumModes(0)!=2)
                {
                    gType = eDeformed;
                }

                m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
                    gType, m_coordim, m_xmap, m_coeffs);
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

                for (i = 0; i < m_coordim; ++i)
                {
                    m_coeffs[i][0] = (*m_verts[0])[i];
                    m_coeffs[i][1] = (*m_verts[1])[i];
                }

                m_state = ePtsFilled;
            }
        }

        NekDouble SegGeom::v_GetLocCoords(
                const Array<OneD, const NekDouble>& coords,
                      Array<OneD,NekDouble>& Lcoords)
        {
            int i;
            NekDouble resid = 0.0;
            SegGeom::v_FillGeom();

            // calculate local coordinate for coord
            if(GetMetricInfo()->GetGtype() == eRegular)
            {
                NekDouble len = 0.0;
                NekDouble xi  = 0.0;

                const int npts = m_xmap->GetTotPoints();
                Array<OneD, NekDouble> pts(npts);

                for(i = 0; i < m_coordim; ++i)
                {
                    m_xmap->BwdTrans(m_coeffs[i], pts);
                    len  += (pts[npts-1]-pts[0])*(pts[npts-1]-pts[0]);
                    xi   += (coords[i]-pts[0])*(coords[i]-pts[0]);
                }

                len = sqrt(len);
                xi  = sqrt(xi);

                Lcoords[0] = 2*xi/len-1.0;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "inverse mapping must be set up to use this call");
            }
            return resid;
        }

        /**
         * @brief Determines if a point specified in global coordinates is
         * located within this tetrahedral geometry.
         */
        bool SegGeom::v_ContainsPoint(
            const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            Array<OneD,NekDouble> locCoord(GetCoordim(),0.0);
            return v_ContainsPoint(gloCoord,locCoord,tol);

        }

        bool SegGeom::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                Array<OneD, NekDouble> &stdCoord,
                NekDouble tol)
        {
            NekDouble resid; 
            return v_ContainsPoint(gloCoord,stdCoord,tol,resid);
        }
        
        bool SegGeom::v_ContainsPoint(
                const Array<OneD, const NekDouble>& gloCoord,
                Array<OneD, NekDouble> &stdCoord,
                NekDouble tol,
                NekDouble &resid)
        {
            resid = GetLocCoords(gloCoord, stdCoord);
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

        PointGeomSharedPtr SegGeom::v_GetVertex(const int i) const
        {
            PointGeomSharedPtr returnval;

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

        const LibUtilities::BasisSharedPtr SegGeom::v_GetBasis(const int i)
        {
            return m_xmap->GetBasis(i);
        }

        StdRegions::StdExpansionSharedPtr SegGeom::v_GetXmap() const
        {
            return m_xmap;
        }

        void SegGeom::v_SetOwnData()
        {
            m_ownData = true;
        }

        StdRegions::Orientation SegGeom::v_GetPorient(const int i) const
        {
            //cout << "StdRegions::PointOrientation GetPorient"<<endl;
            ASSERTL2((i >=0) && (i <= 1),"Point id must be between 0 and 1");

            if (i % 2 == 0)
            {
                return StdRegions::eBwd;
            }
            else 
            {
                return StdRegions::eFwd;
            }
        }


        LibUtilities::ShapeType SegGeom::v_DetShapeType() const
        {
            return  LibUtilities::eSegment;
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

