////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/QuadGeom.cpp,v $
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

#include <StdRegions/StdRegions.hpp>

#include <SpatialDomains/EdgeComponent.h>
#include <SpatialDomains/QuadGeom.h>

#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdExpMap.h>
#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /// Default constructor

        QuadGeom::QuadGeom()
        {
        }

        QuadGeom::QuadGeom(const VertexComponentSharedPtr verts[], 
			   const EdgeComponentSharedPtr edges[], 
			   StdRegions::EdgeOrientation * eorient):
	  QuadFaceComponent(verts[0]->GetCoordim())
        {
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+QuadGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);


            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");
        }
	
	QuadGeom::QuadGeom(const EdgeComponentSharedPtr edges[], 
			   StdRegions::EdgeOrientation * eorient)
        {
	    int j;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);
	    
	    for(j=0; j <kNedges; ++j)
	    {
		if(eorient[j] == StdRegions::eForwards)
		{
		    m_verts.push_back(edges[j]->GetVertex(0));
		}
		else
		{
		    m_verts.push_back(edges[j]->GetVertex(1));
		}
	    }
            
            for (j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");
        }

        QuadGeom::~QuadGeom()
        {
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        GeoFacSharedPtr QuadGeom::GenXGeoFac()
        {
            GeoFacSharedPtr gfac;
            int i;
            StdRegions::GeomType Gtype = StdRegions::eRegular;

            StdRegions::StdExpansion2D ** xmaptemp = new 
                StdRegions::StdExpansion2D*[m_coordim];

            for(i=0;i<m_coordim;i++)
            {
	      xmaptemp[i] = m_xmap[i];
            }

            FillGeom();

            // check to see if expansions are linear
            for(i = 0; i < m_coordim; ++i)
            {
                if((m_xmap[i]->GetBasisOrder(0) != 2)||
                    (m_xmap[i]->GetBasisOrder(1) != 2))
                {
                    Gtype = StdRegions::eDeformed;
                }
            }

            // check to see if all angles are 90 degrees
            if(Gtype == StdRegions::eRegular)
            {
                double dx1,dx2,dy1,dy2;

                for(i = 0; i < 3; ++i)
                {
                    dx1 = m_verts[i+1]->x() - m_verts[i]->x();
                    dy1 = m_verts[i+1]->y() - m_verts[i]->y();

                    dx2 = m_verts[((i+3)%4)]->x() - m_verts[i]->x();
                    dy2 = m_verts[((i+3)%4)]->y() - m_verts[i]->y();

                    if(fabs(dx1*dx2 + dy1*dy2) > 
                        sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2))*
                        kGeomRightAngleTol)
                    {
                        Gtype = StdRegions::eDeformed;
                        break;
                    }
                }
            }

            gfac.reset(new GeoFac(Gtype,  m_coordim, 
                (const StdRegions::StdExpansion2D **) xmaptemp));

            delete[] xmaptemp;

            return gfac; 
        }


        /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts and edges are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */

        void QuadGeom::FillGeom(){

            // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j; 
                int order0 = m_xmap[0]->GetBasisOrder(0);
                int order1 = m_xmap[0]->GetBasisOrder(1);
                double *coef;
                StdRegions::StdExpMap Map,MapV;

                // set side 0 
		m_xmap[0]->MapTo((*m_edges[0])[0]->GetNcoeffs(),
				 (*m_edges[0])[0]->GetBasisType(0),
				 0,m_eorient[0],Map);

                for(i = 0; i < m_coordim; ++i)
                {
                    coef  = (*m_edges[0])[i]->GetCoeffs();
                    for(j = 0; j < order0; ++j)
                    {
                        m_xmap[i]->SetCoeff(Map[j],coef[j]);
                    }
                }

                // set side 2
		m_xmap[0]->MapTo((*m_edges[2])[0]->GetNcoeffs(),
				 (*m_edges[2])[0]->GetBasisType(0),
				 2,m_eorient[2],Map);


                for(i = 0; i < m_coordim; ++i)
                {
                    coef  = (*m_edges[2])[i]->GetCoeffs();
                    for(j = 0; j < order0; ++j)
                    {
                        m_xmap[i]->SetCoeff(Map[j],coef[j]);
                    }
                }

                // set Vertices into side 1
                (*m_edges[1])[0]->MapTo(m_eorient[1],MapV);

                for(i = 0; i < m_coordim; ++i)
                {
                    (*m_edges[1])[i]->SetCoeff(MapV[0],(*m_verts[1])[i]);
                    (*m_edges[1])[i]->SetCoeff(MapV[1],(*m_verts[2])[i]);
                }

                // set Vertices into side 3
                (*m_edges[3])[0]->MapTo(m_eorient[3],MapV);

                for(i = 0; i < m_coordim; ++i)
                {
                    (*m_edges[3])[i]->SetCoeff(MapV[0],(*m_verts[0])[i]);
                    (*m_edges[3])[i]->SetCoeff(MapV[1],(*m_verts[3])[i]);
                }

                // set side 1
		m_xmap[0]->MapTo((*m_edges[1])[0]->GetNcoeffs(),
				 (*m_edges[1])[0]->GetBasisType(0),
				 1,m_eorient[1],Map);

                for(i = 0; i < m_coordim; ++i)
                {
                    coef  = (*m_edges[1])[i]->GetCoeffs();
                    for(j = 0; j < order1; ++j)
                    {
                        m_xmap[i]->SetCoeff(Map[j],coef[j]);
                    }
                }

                // set side 3
		m_xmap[0]->MapTo((*m_edges[3])[0]->GetNcoeffs(),
				 (*m_edges[3])[0]->GetBasisType(0),
				 3,m_eorient[3],Map);

                for(i = 0; i < m_coordim; ++i)
                {
                    coef  = (*m_edges[3])[i]->GetCoeffs();
                    for(j = 0; j < order1; ++j)
                    {
                        m_xmap[i]->SetCoeff(Map[j],coef[j]);
                    }
                }

                for(i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i]->BwdTrans(m_xmap[i]->GetPhys());
                }

                m_state = ePtsFilled;
            }
        }

        void QuadGeom::GetLocCoords(double *Lcoords, const double *coords)
        {
            int i;

            FillGeom();

            // calculate local coordinate for coord 
            if(GetGtype() == StdRegions::eRegular)
            { // can assume it is right angled rectangle
                double len0 = 0.0 ;
                double len1 = 0.0;
                double xi0 = 0.0;
                double xi1 = 0.0;
                const double *pts;
                int nq0, nq1;

                // get points;
                //find end points
                for(i = 0; i < m_coordim; ++i)
                {
                    nq0 = m_xmap[i]->GetPointsOrder(0);
                    nq1 = m_xmap[i]->GetPointsOrder(1);

                    pts = m_xmap[i]->GetPhys();

                    // use projection to side 1 to determine xi_1 coordinate based on length
                    len0 += (pts[nq0-1]-pts[0])*(pts[nq0-1]-pts[0]);	
                    xi0  += (coords[i] -pts[0])*(pts[nq0-1]-pts[0]);

                    // use projection to side 4 to determine xi_2 coordinate based on length
                    len1 += (pts[nq0*(nq1-1)]-pts[0])*(pts[nq0*(nq1-1)]-pts[0]); 
                    xi1  += (coords[i] -pts[0])*(pts[nq0*(nq1-1)]-pts[0]);	
                }

                Lcoords[0] =  2*xi0/len0-1.0;
                Lcoords[1] =  2*xi1/len1-1.0;

            }
            else
            {
	      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
                    "inverse mapping must be set up to use this call");
            }
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: QuadGeom.cpp,v $
// Revision 1.5  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.4  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.3  2006/05/16 22:28:31  sherwin
// Updates to add in FaceComponent call to constructors
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:03  kirby
// *** empty log message ***
//
// Revision 1.23  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.22  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.21  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.20  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.19  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.18  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.17  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.16  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.15  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
