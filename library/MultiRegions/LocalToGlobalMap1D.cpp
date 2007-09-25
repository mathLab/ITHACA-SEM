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
                                               const SpatialDomains::CompositeVector &domain)
        {
            int i,j,k;
            int gid = 0;
            int cnt = 0;
            int cnt2 = 0;
            
            // set up Local to Continuous mapping 
            StdRegions::StdExpMap vmap;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            
            m_totLocLen    = loclen;        
            m_locToContMap = Array<OneD, int>(m_totLocLen,-1);
            
            // set up simple map based on vertex and edge id's
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];

                for(j = 0; j < comp->size(); ++j)
                {
                    locexp[cnt2]->MapTo(StdRegions::eForwards,vmap);
                
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        for(k = 0; k < 2; ++k)
                        {
                            m_locToContMap[cnt + vmap[k]] =  SegmentGeom->GetVid(k);
                            gid = max(gid,m_locToContMap[cnt + vmap[k]]);
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }
                    cnt += locexp[cnt2]->GetNcoeffs();
                    ++cnt2;
                }
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
        
        void LocalToGlobalMap1D::ResetMapping(const int NumDirichlet, const int NumNeumann, const int NumRobin,
                                              const LocalRegions::PointExpVector &bndCond)
        {
            int i,j,cnt;
            int nbnd;
            m_numDirichletBCs = NumDirichlet; 
            m_numNeumannBCs = NumNeumann; 
            m_numRobinBCs = NumRobin; 

            nbnd = bndCond.size();

            Array<OneD, int> oldGlobalID(nbnd);   

            for(i = cnt = 0; i < nbnd; ++i)
            {
                SpatialDomains::VertexComponentSharedPtr vert = (bndCond[i])->GetVertex();
                oldGlobalID[cnt++] = vert->GetVid(); 
            } 

            // Find the index of the BCs entry in the mapping array
            Array<OneD, int> LocalID(nbnd);
            
            for(i = 0; i < m_totLocLen; ++i)
            {        
                for(j = 0; j<nbnd; ++j)
                {
                    if(m_locToContMap[i] == oldGlobalID[j])
                    {
                            LocalID[j] = i;
                            break;
                    }
                }
            }      

            // Reset the mapping   
            if(NumDirichlet)
            {
                bool check;
                int incr;
                for(i = 0; i < m_totLocLen; ++i)
                {
                    check=true;
                    incr=0;
                    for(j = 0; j<NumDirichlet; ++j)
                    {
                        if(m_locToContMap[i] == oldGlobalID[j])
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
                            if(m_locToContMap[i] < oldGlobalID[j])
                            {
                                ++incr;
                            }
                        }
                        m_locToContMap[i]+=incr;
                    }
                }    
            } 
            
            // Store the new global id of the vertices where the BC are imposed
            
            m_bndCondGlobalID = Array<OneD, int>(nbnd);   
            for(i = 0; i < nbnd; ++i)
            {
                m_bndCondGlobalID[i] = m_locToContMap[ LocalID[i] ];
            }
        }
    }
}

/**
* $Log: LocalToGlobalMap1D.cpp,v $
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
