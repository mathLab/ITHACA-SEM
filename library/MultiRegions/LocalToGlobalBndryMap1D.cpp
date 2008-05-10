///////////////////////////////////////////////////////////////////////////////
//
// File LocaltoGlobalBndryMap1D.cpp
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

#include <MultiRegions/LocalToGlobalBndryMap1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        LocalToGlobalBndryMap1D::LocalToGlobalBndryMap1D(const StdRegions::StdExpansionVector &locexp, 
                                                         const SpatialDomains::MeshGraph1D &graph1D,
                                                         const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                                         const Array<OneD, const SpatialDomains::BoundaryConditionType> &bndCondTypes)
        {
            int i,j;
            int vid  = 0;
            int gid  = 0;
            int cnt  = 0;
            int nbnd = bndCondExp.num_elements();
            
            // set up Local to Continuous mapping 
            StdRegions::StdExpMap vmap;
            LocalRegions::SegExpSharedPtr locSegExp;
              
            m_totLocBndDofs   = 2*locexp.size();
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
                        m_locToContBndMap[cnt + j] =  renumbVerts[vid];
                    }    
                    cnt += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }  

            m_totGloBndDofs = gid;
        }
 
        LocalToGlobalBndryMap1D::~LocalToGlobalBndryMap1D()
        {
        }
    }
}

/**
* $Log: LocalToGlobalBndryMap1D.cpp,v $
* Revision 1.3  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.2  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.1  2007/11/20 16:27:16  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
**/
