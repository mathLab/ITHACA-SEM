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
#include "pchSpatialDomains.h"

#include <SpatialDomains/SegGeom.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>

#include <fstream>

namespace Nektar
{
    namespace SpatialDomains
    {
        SegGeom::SegGeom()
        {
            m_GeomShapeType = eSegment;
        }
        
        SegGeom::SegGeom(int id, const int coordim):
            Geometry1D(coordim),
            m_xmap(coordim)
        {
            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
            m_GeomShapeType = eSegment;
            m_eid = id;
            
            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }
        }
        
        SegGeom::SegGeom(int id, const int coordim,
                         const VertexComponentSharedPtr vertex[]): 
            Geometry1D(coordim)
        {
            m_GeomShapeType = eSegment;
            m_eid   = id;
            m_state = eNotFilled;
            
            if (coordim > 0)
            {
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                               LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
                
                m_xmap = Array<OneD, StdRegions::StdExpansion1DSharedPtr>(m_coordim);
                
                for(int i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                }
            }
            
            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }

        SegGeom::SegGeom(int id, const int coordim,
                         const VertexComponentSharedPtr vertex[], 
                         const CurveSharedPtr &curve): 
            Geometry1D(coordim)
        {
            m_GeomShapeType = eSegment;
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
                    std::string err = "Vertex 0 is separated from first point by more than ";
                    std::stringstream strstrm;
                    strstrm << NekConstants::kVertexTheSameDouble; 
                    err += strstrm.str();
                    NEKERROR(ErrorUtil::ewarning, err.c_str()); 
                }

                if(vertex[1]->dist(*(curve->m_points[npts-1])) > NekConstants::kVertexTheSameDouble)
                {
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
                    
                    NekVector<const NekDouble> in(npts,tmp,eWrapper);
                    NekVector<NekDouble>       out(npts+1,m_xmap[i]->UpdatePhys(),eWrapper);
                    out  = (*I0)*in;

                    m_xmap[i]->FwdTrans(m_xmap[i]->GetPhys(),m_xmap[i]->UpdateCoeffs());
                }
            }
            
            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];
        }
        

        SegGeom::SegGeom(const int id, const VertexComponentSharedPtr vert1, 
                         const VertexComponentSharedPtr  vert2):
            Geometry1D(vert1->GetCoordim()), m_xmap(vert1->GetCoordim())
        {
            m_GeomShapeType = eSegment;
            
            m_verts[0] = vert1; 
            m_verts[1] = vert2;
            
            m_state = eNotFilled;
            
            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
            m_eid = id;
            
            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }
        }
        
        SegGeom::SegGeom(const SegGeom &in)
        {
            // From Geometry class
            m_GeomShapeType = in.m_GeomShapeType;
            
            // info from EdgeComponent class
            m_eid     = in.m_eid;
            std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtmap.begin(); def != in.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);    
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

        NekDouble SegGeom::GetCoord(const int i, 
									const Array<OneD, const NekDouble> &Lcoord) 
        {

            ASSERTL1(m_state == ePtsFilled, "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        void SegGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int SegGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool SegGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtmap.end())
            {
                return(true);
            }

            return(false);
        }

        /// \brief Get the orientation of edge1.
        ///
        /// If edge1 is connected to edge2 in the same direction as
        /// the points comprising edge1 then it is forward, otherwise
        /// it is backward.  
        ///
        /// For example, assume edge1 is comprised of points 1 and 2,
        /// and edge2 is comprised of points 2 and 3, then edge1 is
        /// forward.
        ///
        /// If edge1 is comprised of points 2 and 1 and edge2 is
        /// comprised of points 3 and 2, then edge1 is backward.

        StdRegions::EdgeOrientation SegGeom::GetEdgeOrientation(const SegGeom &edge1,  const SegGeom &edge2)
        {
            StdRegions::EdgeOrientation returnval = StdRegions::eForwards;

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

        // Set up GeoFac for this geometry using Coord quadrature distribution
        void SegGeom::GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            SpatialDomains::GeomType gType = eRegular;
            const SpatialDomains::GeomType kDeformedType = eDeformed;
            
            FillGeom();
            
            if(m_xmap[0]->GetBasisNumModes(0)!=2)
            {
                gType = eDeformed;
            }

            m_geomfactors = MemoryManager<GeomFactors1D>::AllocateSharedPtr(gType, m_coordim, m_xmap, tbasis);
        }
        
        
        /** \brief put all quadrature information into edge structure and 
            backward transform */
        
        void SegGeom::FillGeom()
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

        void SegGeom::GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
        {
            int i;
            
            FillGeom();
            
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
        
        void SegGeom::WriteToFile(std::ofstream &outfile, const int dumpVar)
        {
            
            int i,j;
            int  nquad = m_xmap[0]->GetNumPoints(0);
            double *coords[3];
            
            FillGeom();
            
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
                const Array<OneD, const NekDouble> &gloCoord)
        {
            Array<OneD,NekDouble> stdCoord(GetCoordim(),0.0);
            GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -1 && stdCoord[0] <= 1)
            {
                return true;
            }
            return false;
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: SegGeom.cpp,v $
// Revision 1.27  2009/12/15 18:09:02  cantwell
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
// Revision 1.26  2009/07/02 13:32:24  sehunchun
// *** empty log message ***
//
// Revision 1.25  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.24  2008/12/18 14:08:59  pvos
// NekConstants update
//
// Revision 1.23  2008/09/09 14:22:39  sherwin
// Added curved segment constructor methods
//
// Revision 1.22  2008/06/11 21:34:42  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.21  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.20  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.19  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.18  2007/07/27 21:08:30  jfrazier
// Removed use of temporary in call to memory manager allocate.
//
// Revision 1.17  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.16  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.15  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.14  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.13  2007/05/17 18:45:25  jfrazier
// Minor changes to accommodate Array class.
//
// Revision 1.12  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.11  2007/03/29 19:26:36  bnelson
// *** empty log message ***
//
// Revision 1.10  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.9  2007/03/20 09:17:40  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.8  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.7  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.6  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.5  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.4  2006/05/16 20:12:59  jfrazier
// Minor fixes to correct bugs.
//
// Revision 1.3  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.2  2006/05/06 20:36:16  sherwin
// Modifications to get LocalRegions/Project1D working
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.26  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.25  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.24  2006/04/04 23:12:38  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.23  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.22  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.21  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.20  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.19  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
