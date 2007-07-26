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
        LocalToGlobalMap1D::LocalToGlobalMap1D(const int loclen, 
                             const StdRegions::StdExpansionVector &locexp, 
                             const SpatialDomains::Composite &cmps)
    {
        int i,j,gid,cnt;
        
        // set up Local to Continuous mapping 
        StdRegions::StdExpMap vmap;

            m_totLocLen    = loclen;        
            m_locToContMap = Array<OneD, int>(m_totLocLen,-1);

        // set up simple map based on vertex and edge id's
        for(i = cnt = gid = 0; i < locexp.size(); ++i)
        {
        locexp[i]->MapTo(StdRegions::eForwards,vmap);
                
                SpatialDomains::SegGeomSharedPtr SegmentGeom;
                
                if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*cmps)[i]))
                {
                    for(j = 0; j < 2; ++j)
                    {
                        m_locToContMap[cnt + vmap[j]] =  SegmentGeom->GetVid(j);
                        gid = max(gid,m_locToContMap[cnt + vmap[j]]);
                    }
        }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }
                cnt += locexp[i]->GetNcoeffs();
        }
        
            for(i = 0; i < m_totLocLen; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = ++gid;
                }
            }
        m_totGloLen = ++gid;
    }
    
    
    LocalToGlobalMap1D::~LocalToGlobalMap1D()
    {
    }

        void LocalToGlobalMap1D::ResetMapping(const int NumDirichlet, 
                                         SpatialDomains::BoundaryConditions &bcs,
                                              const std::string variable)
        {
            m_numDirichletBCs = NumDirichlet;

            if(NumDirichlet)
            {
                int i,j, nbnd, cnt;
                Array<OneD, int> dbc_id(NumDirichlet);                

                // Find the global ID of the vertices where dirichlet BC are imposed
                SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
                SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

                nbnd = bregions.size();

                for(i = cnt = 0; i < nbnd; ++i)
                {
                    if(  ((*(bconditions[i]))[variable])->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        SpatialDomains::VertexComponentSharedPtr vert;
                        
                        if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                        {
                            dbc_id[cnt++] = vert->GetVid(); 
                        }
                        else
                        {
                            ASSERTL0(false,"dynamic cast to a vertex failed");
                        }                        
                        if(cnt==NumDirichlet)
                        {
                            break;
                        }
                    }
                }

                bool check;
                int incr;
                // Modify the numbering 
                for(i = 0; i < m_totLocLen; ++i)
                {
                    check=true;
                    incr=0;
                    for(j = 0; j<NumDirichlet; ++j)
                    {
                        if(m_locToContMap[i] == dbc_id[j])
                        {
                            m_locToContMap[i] = j;
                            check=false;
                            break;
                        }
                    }
                    if(check)
                    {
                        for(j = 0; j<NumDirichlet; ++j)
                        {
                            if(m_locToContMap[i] < dbc_id[j])
                            {
                                ++incr;
                            }
                        }
                        m_locToContMap[i]=m_locToContMap[i]+incr;
                    }
                }    
            }            
        }
    }
}

/**
* $Log: LocalToGlobalMap1D.cpp,v $
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
