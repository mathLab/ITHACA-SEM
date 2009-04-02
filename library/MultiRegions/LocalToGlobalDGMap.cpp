///////////////////////////////////////////////////////////////////////////////
//
// File LocToGlobalDGMap.cpp
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
// Description: Local to Global Base Class mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalDGMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalDGMap::LocalToGlobalDGMap():
            m_numDirichletBndPhys(0)
        {
        }
        
        LocalToGlobalDGMap::~LocalToGlobalDGMap()
        {
        }

        LocalToGlobalDGMap::LocalToGlobalDGMap( const SpatialDomains::MeshGraph1D &graph1D,
                                                const boost::shared_ptr<StdRegions::StdExpansionVector> &exp1D,
                                                const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond)
        {
            int i,j; 
            int cnt, vid, gid;
            int nbnd = bndCondExp.num_elements();
            
            // set up Local to Continuous mapping 
            Array<OneD,unsigned int> vmap;
            LocalRegions::SegExpSharedPtr locSegExp;
            
            m_numGlobalBndCoeffs  = exp1D->size()+1;

            m_numLocalBndCoeffs = 2*exp1D->size();
            m_localToGlobalBndMap   = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_localToGlobalBndSign  = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);

            m_signChange = true;

            map<int, int> MeshVertToLocalVert;
            
            // Order the Dirichlet vertices first.
            gid = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndCond[i]->GetBoundaryConditionType() ==SpatialDomains::eDirichlet)
                {
                    m_numDirichletBndPhys++;
                    vid = ((bndCondExp[i])->GetVertex())->GetVid();

                    MeshVertToLocalVert[vid] = gid++;
                }
            }
            
            // set up simple map based on vertex and edge id's
            cnt = 0;
            for(i = 0; i < exp1D->size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>((*exp1D)[i]))
                {
                    locSegExp->GetBoundaryMap(vmap);
                
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {   
                        vid = (locSegExp->GetGeom1D())->GetVid(j);

                        if(MeshVertToLocalVert.count(vid) == 0)
                        {
                            MeshVertToLocalVert[vid] = gid++;
                        }   
                        
                        m_localToGlobalBndMap[cnt + j] =  MeshVertToLocalVert.find(vid)->second;
                    }    
                    cnt += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }  

            // Set up boundary mapping
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int >(nbnd);
            m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD, NekDouble >(nbnd,1.0);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys = 0;

            for(i = 0; i < nbnd; ++i)
            {
                vid = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = MeshVertToLocalVert.find(vid)->second;
                
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    m_numLocalDirBndCoeffs += 1;
                    m_numDirichletBndPhys  += 1;
                }
            }

            CalculateBndSystemBandWidth(*exp1D);
        }

        LocalToGlobalDGMap::LocalToGlobalDGMap(SpatialDomains::MeshGraph2D &graph2D, 
                                               const GenExpList1DSharedPtr &trace, 
                                               const boost::shared_ptr<StdRegions::StdExpansionVector> &exp2D, 
                                               const Array<OneD, MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond, 
                                               const map<int,int> &periodicEdges)
        {
            int i,j,k,cnt,id, id1, order_e,gid;
            int ntrace_exp = trace->GetExpSize();
            int nel        = exp2D->size();
            int nbnd = bndCondExp.num_elements();
            LocalRegions::SegExpSharedPtr  locSegExp;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr  locTriExp;
            SpatialDomains::Geometry1DSharedPtr SegGeom;
            
            map<int, int> MeshEdgeId;

            m_signChange = true;
            
            // determine mapping from geometry edges to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(trace->GetExp(i)))
                {
                    id = (locSegExp->GetGeom1D())->GetEid();
                    
                    if(periodicEdges.count(id) > 0)
                    {
                        if(MeshEdgeId.count(id) == 0)
                        {
                            id1 = periodicEdges.find(id)->second;
                            MeshEdgeId[id] = i;
                            MeshEdgeId[id1] = i;
                        }
                    }
                    else
                    {
                        MeshEdgeId[id] = i;
                    }
                }
                else
                {
                    ASSERTL0(false,"Dynamics cast to segment expansion failed");
                }
            }

            // Count total number of edges edges
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*exp2D)[i]->GetNedges();
            }
            
            Array<OneD, StdRegions::StdExpansion1DSharedPtr> edgemap(cnt);
            m_elmtToTrace = Array<OneD, Array<OneD,StdRegions::StdExpansion1DSharedPtr> >(nel);

            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = edgemap + cnt; 

                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {   
                        SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);
                        
                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::GenSegExp> ((*trace).GetExp(MeshEdgeId.find(id)->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {    
                        SegGeom = (locTriExp->GetGeom2D())->GetEdge(j);

                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::GenSegExp> ((*trace).GetExp((MeshEdgeId.find(id))->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }
                
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }                    
                cnt += (*exp2D)[i]->GetNedges();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetExpSize();
            }
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int >(cnt);
            m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD,NekDouble >(cnt,1.0);
            m_bndExpAdjacentOrient = Array<OneD, AdjacentTraceOrientation > (cnt);
            
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys   = 0;
            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        
                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+j] = MeshEdgeId.find(id)->second;
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }

                        // Check to see which way boundar edge is
                        // orientated with respect to connecting
                        // element counter-clockwise convention.

                        SpatialDomains::ElementEdgeVectorSharedPtr con_elmt
                            = graph2D.GetElementsFromEdge(SegGeom);
                        
                        if((boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*con_elmt)[0]->m_Element))->GetEorient((*con_elmt)[0]->m_EdgeIndx) == StdRegions::eForwards)
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsForwards;
                        }
                        else
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsBackwards;
                        }
                        
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local Segment expansion failed");
                    }
                    
                    if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs  += locSegExp->GetNcoeffs();
                        m_numDirichletBndPhys   += locSegExp->GetTotPoints();
                    }
                }
                cnt += j;
            }
            
            // Set up integer mapping array and sign change for each
            // degree of freedom

            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                nbndry += (*exp2D)[i]->NumDGBndryCoeffs();
            }

            m_numLocalBndCoeffs = nbndry;
            m_localToGlobalBndMap  = Array<OneD, int > (nbndry);
            m_localToGlobalBndSign = Array<OneD, NekDouble > (nbndry,1);

            nbndry = cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                nbndry += (*exp2D)[i]->NumDGBndryCoeffs();

                for(j = 0; j < (*exp2D)[i]->GetNedges(); ++j)
                {   
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[i][j]);
                    SegGeom = locSegExp->GetGeom1D();
                    
                    id  = SegGeom->GetEid();
                    gid = trace->GetCoeff_Offset(MeshEdgeId.find(id)->second);
                    
                    order_e = locSegExp->GetNcoeffs();
                    
                    if((*exp2D)[i]->GetEorient(j) == StdRegions::eForwards)
                    {
                        for(k = 0; k < order_e; ++k)
                        {
                            m_localToGlobalBndMap[k+cnt] = gid + k;
                        }
                    }
                    else // backwards orientated
                    {
                        switch(locSegExp->GetBasisType(0))
                        {
                        case LibUtilities::eModified_A:
                            // reverse vertex order
                            m_localToGlobalBndMap[cnt] = gid + 1;
                            m_localToGlobalBndMap[cnt+1] = gid;
                            for(k = 2; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[k+cnt] = gid + k;
                            }

                            // negate odd modes
                            for(k = 3; k < order_e; k+=2)
                            {
                                m_localToGlobalBndSign[cnt+k] = -1.0;
                            }
                            
                            
                            break;
                        case LibUtilities::eGLL_Lagrange:
                            // reverse  order
                            for(k = 0; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                            }
                            break;
                        default:
                            ASSERTL0(false,"Boundary type not permitted");
                            
                        }
                    }
                    cnt += order_e;
                }
            }
            
            m_numGlobalBndCoeffs = trace->GetNcoeffs();

            CalculateBndSystemBandWidth(*exp2D);
        }



        // ----------------------------------------------------------------
        // Calculation of the bandwith ---- The bandwidth here
        // calculated corresponds to what is referred to as
        // half-bandwidth.  If the elements of the matrix are
        // designated as a_ij, it corresponds to the maximum value of
        // |i-j| for non-zero a_ij.  As a result, the value also
        // corresponds to the number of sub or superdiagonals.
        //
        // The bandwith can be calculated elementally as it
        // corresponds to the maximal elemental bandwith (i.e. the
        // maximal difference in global DOF index for every element)
        //
        // 2 different bandwiths can be calculated: - the bandwith of
        // the full global system - the bandwith of the global
        // boundary system (as used for static condensation)
        void LocalToGlobalDGMap::CalculateBndSystemBandWidth(const StdRegions::StdExpansionVector &locExpVector)
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;

            for(i = 0; i < locExpVector.size(); ++i)
            {
                locSize = locExpVector[i]->NumDGBndryCoeffs();
                maxId = -1;
                minId = m_numLocalBndCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalBndMap[cnt+j] >= m_numLocalDirBndCoeffs)
                    {
                        if(m_localToGlobalBndMap[cnt+j] > maxId)
                        {
                            maxId = m_localToGlobalBndMap[cnt+j];
                        }
                        
                        if(m_localToGlobalBndMap[cnt+j] < minId)
                        {
                            minId = m_localToGlobalBndMap[cnt+j];
                        }
                    }
                }
                bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

                cnt+=locSize;
            }

            m_bndSystemBandWidth = bwidth;
        }
        
    }


}
