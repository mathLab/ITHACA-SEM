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

        LocalToGlobalBndryMap1D::LocalToGlobalBndryMap1D(const int NumDirichlet,
                                SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                StdRegions::StdExpansionVector &locexp,
                                const SpatialDomains::CompositeVector &domain)
        {
            int i,j,k;
            int gid = 0;
            int cnt = 0;
            int cnt1 = 0;
            int nbnd; 

            StdRegions::StdExpMap vmap;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            
            SpatialDomains::BoundaryRegionCollection    &bregions    = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            Array<OneD,int> DirMap(2,-1);

            m_totLocBndDofs    = 2*locexp.size();
            m_locToContBndMap  = Array<OneD, int>(m_totLocBndDofs,-1);
            
            m_numDirichletBCs = NumDirichlet; 
            
            nbnd = bregions.size();

            // set up simple map based on vertex and edge id's
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];

                for(j = 0; j < comp->size(); ++j)
                {
                    locexp[cnt1]->MapTo(StdRegions::eForwards,vmap);
                
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        for(k = 0; k < 2; ++k)
                        {
                            m_locToContBndMap[cnt + vmap[k]] =  SegmentGeom->GetVid(k);
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }
                    cnt += locexp[cnt1]->NumBndryCoeffs();
                    ++cnt1;
                }
            }
            
            m_totGloBndDofs = locexp.size()+1;

            // reshuffle Dirichlet boundary conditions
            for(i = cnt = 0; i < nbnd; ++i)
            {
                if(((*(bconditions[i]))[variable])->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    SpatialDomains::VertexComponentSharedPtr vert;
                    
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        DirMap[cnt++] = vert->GetVid(); 
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex failed");
                    }
                }
            }

            // Reset mapping 
            for(i = 0; i < NumDirichlet; ++i)
            {
                for(j = 0; j < m_totLocBndDofs; ++j)
                {
                    if(m_locToContBndMap[j] == DirMap[i])
                    {
                        m_locToContBndMap[j] = i + m_totLocBndDofs; 
                    }
                    else if(m_locToContBndMap[j] < DirMap[i])
                    {
                        m_locToContBndMap[j] += 1;
                    }
                }
            }

            for(i = 0; i < m_totLocBndDofs; ++i)
            {
                if(m_locToContBndMap[i] >= m_totLocBndDofs)
                {
                    m_locToContBndMap[i] -= m_totLocBndDofs;
                }
            }
        }

        LocalToGlobalBndryMap1D::~LocalToGlobalBndryMap1D()
        {
        }
    }
}

/**
* $Log: LocalToGlobalMap1D.cpp,v $
**/
