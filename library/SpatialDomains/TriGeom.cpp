////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/TriGeom.cpp,v $
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

#include <SpatialDomains/TriGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        TriGeom::TriGeom()
        {
        }

        TriGeom::TriGeom(const VertexComponentSharedPtr verts[], 
             const EdgeComponentSharedPtr edges[], 
             const StdRegions::EdgeOrientation eorient[]):
        TriFaceComponent(verts[0]->GetCoordim())
        {
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+TriGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TriGeom::kNedges);

            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");
        }

        TriGeom::TriGeom(const EdgeComponentSharedPtr edges[], 
                         const StdRegions::EdgeOrientation eorient[]) :
            TriFaceComponent(edges[0]->GetCoordim())
        {
            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TriGeom::kNedges);

            for(int j=0; j <kNedges; ++j)
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

            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,"Cannot call function with dim == 1");
        }

        TriGeom::~TriGeom()
        {
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        void TriGeom::GenGeomFactors(void)
        {
            GeomType Gtype = eRegular;
            
            FillGeom();

            // check to see if expansions are linear
            for(int i = 0; i < m_coordim; ++i)
            {
                if((m_xmap[i]->GetBasisNumModes(0) != 2)||
                   (m_xmap[i]->GetBasisNumModes(1) != 2))
                {
                    Gtype = eDeformed;
                }
            }

            m_geomfactors = MemoryManager<GeomFactors>::AllocateSharedPtr(Gtype, m_coordim, m_xmap);

        }

        /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts and edges are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */

        void TriGeom::FillGeom()
        {
            // check to see if geometry structure is already filled
            if(m_state == ePtsFilled)
            {
                return;
            }

            int i,j; 
            int order0 = m_xmap[0]->GetBasisNumModes(0);
            int order1 = m_xmap[0]->GetBasisNumModes(1);
            ConstArray<OneD,NekDouble> coef;
            StdRegions::StdExpMap Map,MapV;

            // set side 1 
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

            // set Vertices into side 2 
            (*m_edges[1])[0]->MapTo(m_eorient[1],MapV);

            for(i = 0; i < m_coordim; ++i)
            {
                (*m_edges[1])[i]->SetCoeff(MapV[0],(*m_verts[1])[i]);
                (*m_edges[1])[i]->SetCoeff(MapV[1],(*m_verts[2])[i]);
            }

            // set Vertices into side 3
            (*m_edges[2])[0]->MapTo(m_eorient[2],MapV);

            for(i = 0; i < m_coordim; ++i)
            {
                (*m_edges[2])[i]->SetCoeff(MapV[0],(*m_verts[0])[i]);
                (*m_edges[2])[i]->SetCoeff(MapV[1],(*m_verts[2])[i]);
            }

            // set side 2
            m_xmap[0]->MapTo((*m_edges[1])[0]->GetNcoeffs(),
                 (*m_edges[1])[0]->GetBasisType(0),
                 1,(m_eorient[1]),Map);
        
            for(i = 0; i < m_coordim; ++i)
            {
                coef  = (*m_edges[1])[i]->GetCoeffs();
                for(j = 0; j < order1; ++j)
                {
                    m_xmap[i]->SetCoeff(Map[j],coef[j]);
                }
            }

            // set side 3
            m_xmap[0]->MapTo((*m_edges[2])[0]->GetNcoeffs(),
                 (*m_edges[2])[0]->GetBasisType(0),
                 2,(m_eorient[2]),Map);

            for(i = 0; i < m_coordim; ++i)
            {
                coef  = (*m_edges[2])[i]->GetCoeffs();
                for(j = 0; j < order0; ++j)
                {
                    m_xmap[i]->SetCoeff(Map[j],coef[j]);
                }
            }

            for(i = 0; i < m_coordim; ++i)
            {
                m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(),
                                    m_xmap[i]->UpdatePhys());
            }

            m_state = ePtsFilled;

        }

        void TriGeom::GetLocCoords(const ConstArray<OneD,NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
        {
            FillGeom();

            // calculate local coordinate for coord 
            if(GetGtype() == eRegular)
            { // can assume it is right angled rectangle
                VertexComponent dv1, dv2, norm, orth1, orth2;
                VertexComponent *xin;
                
                switch(m_coordim)
                {
                case 2:
                    xin = new VertexComponent (m_coordim,0,coords[0], coords[1], 0.0);
                    break;
                case 3:
                    xin = new VertexComponent (m_coordim,0,coords[0], coords[1], coords[2]);
                    break;
                }

                dv1.Sub(*m_verts[1],*m_verts[0]);
                dv2.Sub(*m_verts[2],*m_verts[0]);

                norm.Mult(dv1,dv2);

                orth1.Mult(norm,dv1);
                orth2.Mult(norm,dv2);

                xin[0] *= 2.0;
                xin[0] -= *m_verts[1];
                xin[0] -= *m_verts[2];

                Lcoords[0] = xin->dot(orth2)/dv1.dot(orth2);
                Lcoords[1] = xin->dot(orth1)/dv2.dot(orth1);

                delete xin;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    "inverse mapping must be set up to use this call");
            }
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: TriGeom.cpp,v $
// Revision 1.10  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.9  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.8  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.7  2006/07/02 17:16:18  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.6  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.5  2006/06/01 14:15:31  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/23 19:56:33  jfrazier
// These build and run, but the expansion pieces are commented out
// because they would not run.
//
// Revision 1.3  2006/05/16 22:28:31  sherwin
// Updates to add in FaceComponent call to constructors
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:05  kirby
// *** empty log message ***
//
// Revision 1.19  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.18  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.17  2006/04/09 02:08:36  jfrazier
// Added precompiled header.
//
// Revision 1.16  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.15  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.14  2006/03/12 14:20:44  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.12  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.11  2006/02/19 01:37:35  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
