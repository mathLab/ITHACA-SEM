////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/SegGeom.cpp,v $
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
        SegGeom::SegGeom():
            m_owndata(false)
        {
        }

        SegGeom::SegGeom(int id, VertexComponentSharedPtr vert1, VertexComponentSharedPtr  vert2):
            EdgeComponent(id,vert1->GetCoordim())
        {
            m_owndata = false;
            m_verts[0] = vert1; 
            m_verts[1] = vert2;

            m_state = eNotFilled;
        }

        SegGeom::SegGeom(SegGeom &in)
        {

            // info from EdgeComponent class
            m_eid     = in.m_eid;
            m_elmtmap = in.m_elmtmap;
            m_xmap    = in.m_xmap;

            // info from SegGeom class
            m_owndata  = false; 
            m_coordim  = in.m_coordim;
            m_verts[0] = in.m_verts[0]; 
            m_verts[1] = in.m_verts[1];

            m_state = in.m_state;
        }

        SegGeom::~SegGeom()
        {
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution
        void SegGeom::GenGeomFactors(void)
        {
            GeomFactorsSharedPtr gfac;

            StdRegions::StdExpansion1D ** xmaptemp =
                new StdRegions::StdExpansion1D*[m_coordim];

            for(int i=0;i<m_coordim;++i)
            {
                xmaptemp[i] = m_xmap[i];
            }

            FillGeom();

            if(m_xmap[0]->GetBasisNumModes(0)==2)
            {// assumes all direction have same order
                m_geomfactors.reset( new GeomFactors(eRegular, m_coordim, 
                    (const StdRegions::StdExpansion1D **) xmaptemp));
            }
            else
            {
                m_geomfactors.reset(new GeomFactors(eDeformed, m_coordim, 
                    (const StdRegions::StdExpansion1D **) xmaptemp));
            }

            delete[] xmaptemp;
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
                    m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(),m_xmap[i]->GetPhys());
                }
                m_state = ePtsFilled;
            }
        }


        void SegGeom::GetLocCoords(double *Lcoords,const double *coords)
        {
            int i;
	    
            FillGeom();
	    
            // calculate local coordinate for coord
            if(GetGtype() == eRegular)
            {
                const double *pts;
                double len = 0.0;
                double xi  = 0.0;
                int nq;

                // get points;
                //find end points 
                for(i = 0; i < m_coordim; ++i)
                {
                    nq   = m_xmap[i]->GetNumPoints(0);
                    pts  = &(m_xmap[i]->GetPhys())[0];
                    len  += (pts[nq-1]-pts[0])*(pts[nq-1]-pts[0]);	
                    xi   += (coords[i]-pts[0])*(coords[i]-pts[0]);
                }

                len = sqrt(len);
                xi  = sqrt(xi);

                Lcoords[0] =  2*xi/len-1.0;
            }
            else
            {
                ErrorUtil::Error(ErrorUtil::efatal,__FILE__, __LINE__,
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
                coords[i] = &(m_xmap[i]->GetPhys())[0];
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

    }; //end of namespace
}; //end of namespace

//
// $Log: SegGeom.cpp,v $
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
