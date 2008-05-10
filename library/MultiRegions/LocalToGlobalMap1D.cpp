///////////////////////////////////////////////////////////////////////////////
//
// File LocaltoGlobalMap1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Local to Global mapping routines in 1D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalMap1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalMap1D::LocalToGlobalMap1D()
        {
        }

        LocalToGlobalMap1D::LocalToGlobalMap1D(const int loclen, 
                                               const StdRegions::StdExpansionVector &locexp, 
                                               const SpatialDomains::MeshGraph1D &graph1D)
        {
            int i,j;
            int vid = 0;
            int gid = 0;
            int cnt = 0;
            int cnt1 = 0;
            
            // set up Local to Continuous mapping 
            StdRegions::StdExpMap vmap;
            LocalRegions::SegExpSharedPtr locSegExp;
            
            m_totLocDofs      = loclen;        
            m_totLocBndDofs   = 2*locexp.size();
            m_locToContMap    = Array<OneD, int>(m_totLocDofs,-1);
            m_locToContBndMap = Array<OneD, int>(m_totLocBndDofs,-1);

            // Reserve storage for the re-ordering (give a new vertex) of the vertices.
            // This is needed because the set-up of the mapping is based on the vertex id's, but
            // because the domain does not necessarily encompassed the entire meshgraph, the original id's
            // might be unsuitable.
            Array<OneD, int> renumbVerts(graph1D.GetNvertices(),-1);
            
            // set up simple map based on vertex and edge id's
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locexp[i]))
                {
                    locSegExp->MapTo(StdRegions::eForwards,vmap);
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {   
                        vid = (locSegExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = gid++;
                        }   
                        m_locToContMap[cnt     + vmap[j]] =  renumbVerts[vid];
                        m_locToContBndMap[cnt1 + vmap[j]] =  renumbVerts[vid];
                    }    
                    cnt  += locSegExp->GetNcoeffs();
                    cnt1 += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }  

            m_totGloBndDofs = gid;
            
            for(i = 0; i < m_totLocDofs; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = gid++;
                }
            }
            m_totGloDofs = gid;
        }
        

        LocalToGlobalMap1D::LocalToGlobalMap1D(const int loclen, 
                                               const StdRegions::StdExpansionVector &locexp, 
                                               const SpatialDomains::MeshGraph1D &graph1D,
                                               const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                               const Array<OneD, const SpatialDomains::BoundaryConditionType> &bndCondTypes)
        {
            int i,j;
            int vid = 0;
            int gid = 0;
            int cnt = 0;
            int cnt1 = 0;
            int nbnd = bndCondExp.num_elements();
            
            // set up Local to Continuous mapping 
            StdRegions::StdExpMap vmap;
            LocalRegions::SegExpSharedPtr locSegExp;
            
            m_totLocDofs      = loclen;        
            m_totLocBndDofs   = 2*locexp.size();
            m_locToContMap    = Array<OneD, int>(m_totLocDofs,-1);
            m_locToContBndMap = Array<OneD, int>(m_totLocBndDofs,-1);

            m_locToContBndCondMap = Array<OneD, int>(nbnd);   

            // re-order the vertices (as domain does not necessarily contains the entire meshgraph)
            Array<OneD, int> renumbVerts(graph1D.GetNvertices(),-1);

            // This array is used to indicate whether an vertex is part of the boundary of the domain.
            Array<OneD, int> bndCondVertID(graph1D.GetNvertices(),-1);

            // Order the Dirichlet vertices first.
            m_numDirichletDofs = 0;
            for(i = 0; i < nbnd; i++)
            {
                vid = ((bndCondExp[i])->GetVertex())->GetVid();
                bndCondVertID[vid] = i;

                if(bndCondTypes[i]==SpatialDomains::eDirichlet)
                {
                    m_numDirichletDofs++;
                    if(renumbVerts[vid]==-1)
                    {
                        renumbVerts[vid] = gid++;
                    }                    
                }
            }
            
            // set up simple map based on vertex and edge id's
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locexp[i]))
                {
                    locSegExp->MapTo(StdRegions::eForwards,vmap);
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {   
                        vid = (locSegExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = gid++;
                        }   
                        if(bndCondVertID[vid] != -1)
                        {
                            m_locToContBndCondMap[bndCondVertID[vid]] = renumbVerts[vid];
                        }
                        m_locToContMap[cnt     + vmap[j]] =  renumbVerts[vid];
                        m_locToContBndMap[cnt1 + j] =  renumbVerts[vid];
                    }    
                    cnt  += locSegExp->GetNcoeffs();
                    cnt1 += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }  

            m_totGloBndDofs = gid;
            
            for(i = 0; i < m_totLocDofs; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = gid++;
                }
            }
            m_totGloDofs = gid;
        }
        
        LocalToGlobalMap1D::~LocalToGlobalMap1D()
        {
        }
    }
}

/**
* $Log: LocalToGlobalMap1D.cpp,v $
* Revision 1.23  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.22  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.21  2007/11/20 16:27:16  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
* Revision 1.20  2007/10/04 13:57:01  pvos
* fixed some more errors
*
* Revision 1.19  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.18  2007/10/04 11:01:31  pvos
* fixed some errors
*
* Revision 1.17  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.16  2007/09/25 14:25:30  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.15  2007/08/13 14:36:36  pvos
* Neumann BC update
*
* Revision 1.14  2007/08/13 11:09:42  pvos
* Implementation of Neumann BC
*
* Revision 1.13  2007/07/26 08:40:50  sherwin
* Update to use generalised i/o hooks in Helmholtz1D
*
* Revision 1.12  2007/07/23 09:13:57  sherwin
* Update for name change where we removed 'type' from the end
*
* Revision 1.11  2007/07/22 23:04:21  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.10  2007/07/20 02:04:13  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/07/10 08:54:30  pvos
* Updated ContField1D constructor
*
* Revision 1.8  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.7  2007/06/17 19:01:29  bnelson
* Removed unused variables.
*
* Revision 1.6  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
