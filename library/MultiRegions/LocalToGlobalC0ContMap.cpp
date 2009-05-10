///////////////////////////////////////////////////////////////////////////////
//
// File LocToGlobalC0ContMap.cpp
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
// Description: C0-continuous Local to Global mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(void):
            m_numLocalCoeffs(0),
            m_numGlobalCoeffs(0)
        {
        }

        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                                       const StdRegions::StdExpansionVector &locExpVector)
        {
            switch(locExpVector[0]->GetShapeDimension())
            {
            case 1:
                {
                    SetUp1DExpansionC0ContMap(numLocalCoeffs, locExpVector);
                }
                break;
            case 2:
                {
                    SetUp2DExpansionC0ContMap(numLocalCoeffs, locExpVector);
                }
                break;
            case 3:
                {
                    SetUp3DExpansionC0ContMap(numLocalCoeffs, locExpVector);
                }
                break;
            default:
                {
                    ASSERTL0(false,"Local To Global map not defined for this dimension");
                }
            }

            CalculateBndSystemBandWidth(locExpVector);
            CalculateFullSystemBandWidth(locExpVector);

            //cout << "bwidth: " << m_bndSystemBandWidth << endl;
        }

        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                                       const StdRegions::StdExpansionVector &locExpVector, 
                                                       const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                                       const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                       const map<int,int>& periodicVerticesId)
        {
            SetUp1DExpansionC0ContMap(numLocalCoeffs, locExpVector, bndCondExp, bndConditions, periodicVerticesId);
            CalculateBndSystemBandWidth(locExpVector);
            CalculateFullSystemBandWidth(locExpVector);
            //cout << "bwidth: " << m_bndSystemBandWidth << endl;
        }

        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                                       const StdRegions::StdExpansionVector &locExpVector, 
                                                       const Array<OneD, const ExpList1DSharedPtr> &bndCondExp,
                                                       const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                       const vector<map<int,int> >& periodicVerticesId,
                                                       const map<int,int>& periodicEdgesId)
        {
            SetUp2DExpansionC0ContMap(numLocalCoeffs, locExpVector, bndCondExp, bndConditions, periodicVerticesId, periodicEdgesId);
            CalculateBndSystemBandWidth(locExpVector);
            CalculateFullSystemBandWidth(locExpVector);
            //cout << "bwidth: " << m_bndSystemBandWidth << endl;
        }

        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                                       const StdRegions::StdExpansionVector &locExpVector, 
                                                       const Array<OneD, const ExpList2DSharedPtr> &bndCondExp,
                                                       const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                       const map<int,int>& periodicVerticesId,
                                                       const map<int,int>& periodicEdgesId,
                                                       const map<int,int>& periodicFacesId)
        {
            SetUp3DExpansionC0ContMap(numLocalCoeffs, locExpVector, bndCondExp, bndConditions, periodicVerticesId, periodicEdgesId, periodicFacesId);
            CalculateBndSystemBandWidth(locExpVector);
            CalculateFullSystemBandWidth(locExpVector);
        }
        
        LocalToGlobalC0ContMap::~LocalToGlobalC0ContMap(void)
        {
        }
            
        void LocalToGlobalC0ContMap::SetUp1DExpansionC0ContMap(const int numLocalCoeffs, 
                                                               const StdRegions::StdExpansionVector &locExpVector, 
                                                               const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                                               const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                               const map<int,int>& periodicVerticesId)
        {
            int i,j;
            int cnt = 0;
            int meshVertId;
            int meshVertId2;
            int graphVertId = 0;
            int globalId;

            LocalRegions::SegExpSharedPtr locSegExp;

            int nbnd = bndCondExp.num_elements();
            
            m_signChange                     = false;
            m_numLocalCoeffs                 = numLocalCoeffs;        
            m_numLocalBndCoeffs              = 2*locExpVector.size();
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(nbnd);   

            // The only unique identifiers of the vertices of the mesh are the
            // vertex id (stored in their corresponding Geometry object).
            // However, setting up a global numbering based on these id's will not lead to 
            // a suitable or optimal numbering. Mainly because:
            // - we want the Dirichlet DOF's to be listed first
            // - we want an optimal global numbering of the remaining DOF's (strategy still
            //   need to be defined but can for example be: minimum bandwith or minimum fill-in
            //   of the resulting global system matrix)
            //
            // That's why the vertices be rearranged. Currently, this is done 
            // in the following way: The vertices of the mesh are considered as vertices 
            // of a graph. (We then will use algorithms to reorder these 
            // graph-vertices - although not implemented yet in 1D).
            //
            // A container is used to store the graph vertex id's of the different
            // mesh vertices. This is implemented as a STL map such that the graph vertex id
            // can later be retrieved by the unique mesh vertex id's which serve as the key of the map.  
            map<int, int> vertReorderedGraphVertId;
            map<int,int>::const_iterator mapConstIt;

            // STEP 1: Order the Dirichlet vertices first
            m_numGlobalDirBndCoeffs = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                {
                    meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                    vertReorderedGraphVertId[meshVertId] = graphVertId++; 
                    m_numGlobalDirBndCoeffs++;
                    m_numLocalDirBndCoeffs++;
                }
            }

            // STEP 2: Order the periodic vertices next
            // This allows to give corresponding DOF's the same
            // global ID
            for(mapConstIt = periodicVerticesId.begin(); mapConstIt != periodicVerticesId.end(); mapConstIt++)
            {
                meshVertId  = mapConstIt->first;
                meshVertId2 = mapConstIt->second;
                
                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");
                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");

                vertReorderedGraphVertId[meshVertId]  = graphVertId;
                vertReorderedGraphVertId[meshVertId2] = graphVertId++;
            }
            

            // STEP 3: List the other vertices
            // STEP 4: Set up simple map based on vertex and edge id's
            for(i = 0; i < locExpVector.size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locExpVector[i]))
                {
                    for(j = 0; j < locSegExp->GetNverts(); ++j) 
                    {
                        meshVertId = (locSegExp->GetGeom1D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                        {
                            vertReorderedGraphVertId[meshVertId] = graphVertId++;
                        }

                        m_localToGlobalMap[cnt+locSegExp->GetVertexMap(j)] = 
                            vertReorderedGraphVertId[meshVertId];
                    }
                    cnt  += locSegExp->GetNcoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            } 

            // Set up the mapping for the boundary conditions
            for(i = 0; i < nbnd; i++)
            {
                meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = vertReorderedGraphVertId[meshVertId];
            }

            // Setup interior mapping and the boundary map 
            globalId = graphVertId;
            m_numGlobalBndCoeffs = globalId;
            
            cnt=0;
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }
            m_numGlobalCoeffs = globalId;  
        }

        void LocalToGlobalC0ContMap::SetUp2DExpansionC0ContMap(const int numLocalCoeffs, 
                                                               const StdRegions::StdExpansionVector &locExpVector, 
                                                               const Array<OneD, const MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                                               const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                               const vector<map<int,int> >& periodicVerticesId,
                                                               const map<int,int>& periodicEdgesId)
        {
            int i,j,k;
            int cnt = 0;
            int bndEdgeCnt;
            int meshVertId, meshVertId2;
            int meshEdgeId, meshEdgeId2;
            int globalId;
            int nEdgeCoeffs;
            int nEdgeInteriorCoeffs;
            int firstNonDirGraphVertId;
            int nLocBndCondDofs = 0;
            int nLocDirBndCondDofs = 0;
            int graphVertId = 0;
            StdRegions::StdExpansion2DSharedPtr locExpansion;  
            LocalRegions::SegExpSharedPtr       bndSegExp;
            LibUtilities::BasisType             bType;  
            StdRegions::EdgeOrientation         edgeOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;  
            Array<OneD, int>                    edgeInteriorSign;  

            m_signChange = false;

            // The only unique identifiers of the vertices and edges
            // of the mesh are the vertex id and the mesh id (stored
            // in their corresponding Geometry object).  However,
            // setting up a global numbering based on these id's will
            // not lead to a suitable or optimal numbering. Mainly
            // because: - we want the Dirichlet DOF's to be listed
            // first - we want an optimal global numbering of the
            // remaining DOF's (strategy still need to be defined but
            // can for example be: minimum bandwith or minimum fill-in
            // of the resulting global system matrix)
            //
            // That's why the vertices annd egdes will be
            // rearranged. Currently, this is done in the following
            // way: The vertices and edges of the mesh are considered
            // as vertices of a graph (in a computer science
            // way)(equivalently, they can also be considered as
            // boundary degrees of freedom, whereby all boundary modes
            // of a single edge are considered as a single DOF). We
            // then will use algorithms to reorder these
            // graph-vertices (or boundary DOF's).
            //
            // Two different containers are used to store the graph
            // vertex id's of the different mesh vertices and
            // edges. They are implemented as a STL map such that the
            // graph vertex id can later be retrieved by the unique
            // mesh vertex or edge id's which serve as the key of the
            // map.
            map<int,int> vertReorderedGraphVertId;
            map<int,int> edgeReorderedGraphVertId; 
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;
           
            // STEP 1: Order the Dirichlet vertices and edges first
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {                               
                    bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {          
                        meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();  
                        edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;  
                        for(k = 0; k < 2; k++)
                        {
                            meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                        }
                        nLocDirBndCondDofs += bndSegExp->GetNcoeffs();
                    }
                    nLocBndCondDofs += bndSegExp->GetNcoeffs();
                }
            }
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs;
            firstNonDirGraphVertId = graphVertId;

            // STEP 2: Now order all other vertices and edges in the graph
#ifdef NEKTAR_USING_METIS
            // Possibility 1: Use METIS to obtain a gloabl matrix ordering which leads 
            // to minimimal fill in of the factorised matrix (by Cholesky decomposition)
            {      
                int tempGraphVertId = 0;
                int adjncySize = 0;
                int nVerts;
                int vertCnt;
                int edgeCnt;
                map<int, int>          vertTempGraphVertId;
                map<int, int>          edgeTempGraphVertId;
                map<int, set<int> >    graphAdjncySet;  
                set<int>::iterator     setIt;
                Array<OneD, int>       localVerts;
                Array<OneD, int>       localEdges;               

                m_numLocalBndCoeffs = 0;
                // First we are going to set up a temporary ordering
                // of the mesh vertices and edges in the graph. We
                // will then later use METIS to optimise this ordering

                // List the periodic vertices and edges next.  This
                // allows us to give corresponding DOF's the same
                // global ID a) periodic vertices

                for(k = 0; k < PeriodicVerticesId.size(); ++k)
                {
                    for(mapConstIt = periodicVerticesId[k].begin(); 
                        mapConstIt != periodicVerticesId[k].end();
                        mapConstIt++)
                    {
                        meshVertId  = mapConstIt->first;
                        meshVertId2 = mapConstIt->second;

                        if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                        {
                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                            {
                                vertTempGraphVertId[meshVertId]  = tempGraphVertId;
                                vertTempGraphVertId[meshVertId2] = tempGraphVertId++;
                            }
                            else
                            {
                                vertReorderedGraphVertId[meshVertId] = vertReorderedGraphVertId[meshVertId2];
                            }
                        }
                        else
                        {
                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId2] = vertReorderedGraphVertId[meshVertId];
                            }
                            else // Doubly periodic vertex case. 
                            {
                                int id1 = vertTempGraphVertId[meshVertId];
                                int id2 = vertTempGraphVertId[meshVertId2];
                                int id;
                                
                                if(id1 != id2)
                                {
                                    // Reset any values set to id2 to id1. In
                                    // addition if local id is greater than
                                    // id2 decrement list
                                    for(i = 0; i < vertTempGraphVertId.size(); ++i)
                                    {
                                        id = vertTempGraphVertId[i];
                                        if(id == id2)
                                        {
                                            vertTempGraphVertId[i] = id1;
                                        }
                                        else if (id > id2)
                                        {
                                            vertTempGraphVertId[i] = id-1;
                                        }
                                    }
                                    tempGraphVertId--;
                                }
                            }
                        }
                    }
                }
                    
                // b) periodic edges
                for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
                {
                    meshEdgeId  = mapConstIt->first;
                    meshEdgeId2 = mapConstIt->second;
                    
                    if(meshEdgeId < meshEdgeId2)
                    {
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId) == 0,
                                 "This periodic boundary edge has been specified before");
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId2) == 0,
                                 "This periodic boundary edge has been specified before");
                        
                        edgeTempGraphVertId[meshEdgeId]  = tempGraphVertId;
                        edgeTempGraphVertId[meshEdgeId2] = tempGraphVertId++;
                    }
                }

                // List all other vertices and edges
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();
                        
                        nVerts = locExpansion->GetNverts();
                        // For element i, store the temporary graph vertex id's of all element
                        // edges and verices in these 2 arrays below
                        localVerts = Array<OneD, int>(nVerts,-1);
                        localEdges = Array<OneD, int>(nVerts,-1);
                        vertCnt = 0;
                        edgeCnt = 0;
                        
                        // List the (non-dirichlet) vertices and edges
                        // of the mesh as the vertices of the
                        // temporary graph
                        //
                        // At the same time, we can already start
                        // setting ip the adjacency information of the
                        // graph.  This is required as input of the
                        // METIS routines. This adjacency information
                        // should list all adjacent graph-vertices for
                        // all the separate graph-vertices
                        for(j = 0; j < nVerts; ++j) 
                        {   
                            meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {                                
                                if(vertTempGraphVertId.count(meshVertId) == 0)
                                {
                                    vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                                }
                                localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                            }                             
                            
                            meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {                            
                                if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                                {
                                    edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                                }
                                localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            }  
                        }

                        // Now loop over all local edges and vertices
                        // of this element and define that all other
                        // edges and verices of this element are
                        // adjacent to them.  To do so, we use an STL
                        // map where the key is the unique temporary
                        // graph-vertex id and the mapped values is an
                        // STL set which contains the temporary
                        // graph-vertex id's of all adajacent graph
                        // vertices. By looping over all elements in
                        // the mesh and everytime appending the
                        // adjacent grapgh vertices, we will make sure
                        // that all the adjacency information will be
                        // complete.
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localVerts[j]==-1)
                            {
                                break;
                            }
                            for(k = 0; k < nVerts; k++)
                            {
                                if(localVerts[k]==-1)
                                {
                                    break;
                                }
                                if(k!=j)
                                {
                                    graphAdjncySet[ localVerts[j] ].insert(localVerts[k]);
                                }
                            }
                            for(k = 0; k < nVerts; k++)
                            {     
                                if(localEdges[k]==-1)
                                {
                                    break;
                                }                       
                                graphAdjncySet[ localVerts[j] ].insert(localEdges[k]);                            
                            }
                        }
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localEdges[j]==-1)
                            {
                                break;
                            }
                            for(k = 0; k < nVerts; k++)
                            {
                                if(localEdges[k]==-1)
                                {
                                    break;
                                }
                                if(k!=j)
                                {
                                    graphAdjncySet[ localEdges[j] ].insert(localEdges[k]);
                                }
                            }
                            for(k = 0; k < nVerts; k++)
                            {                            
                                if(localVerts[k]==-1)
                                {
                                    break;
                                }
                                graphAdjncySet[ localEdges[j] ].insert(localVerts[k]);                            
                            }
                        }                        
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                    }
                }
                
                for( i = 0; i < tempGraphVertId; i++ )
                {
                    adjncySize += graphAdjncySet[i].size();
                }

                // Now, translate this adjacency information to a
                // format to the format that METIS understands.
                cnt = 0;
                Array<OneD, int> xadj(tempGraphVertId+1);
                Array<OneD, int> adjncy(adjncySize);
                xadj[0]=0;

                for( i = 0; i < tempGraphVertId; i++ )
                {
                    xadj[i+1] = xadj[i] + graphAdjncySet[i].size();
                    for(setIt=(graphAdjncySet[i]).begin(); setIt!=(graphAdjncySet[i]).end(); setIt++)
                    {
                        adjncy[cnt++] = *setIt;
                    }
                }

                // Call METIS to reorder the graph-vertices
                Array<OneD,int> perm(tempGraphVertId);
                Array<OneD,int> iperm(tempGraphVertId);

                Metis::onmetis(tempGraphVertId,xadj,adjncy,perm,iperm);

                // Fill the vertReorderedGraphVertId and
                // edgeReorderedGraphVertId with the optimal ordering
                // from METIS
                for(mapIt = vertTempGraphVertId.begin(); mapIt != vertTempGraphVertId.end(); mapIt++)
                {
                    vertReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId; 
                }
                for(mapIt = edgeTempGraphVertId.begin(); mapIt != edgeTempGraphVertId.end(); mapIt++)
                {
                    edgeReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId; 
                }            
            } 
#else //NEKTAR_USING_METIS
            // Possibility 2: Minimize the bandwith using a reversed
            // Cuthill-McKee algorithm (as provided by the Boost Graph
            // Library)
            if(true)
            {
                // the first template parameter (=OutEdgeList) is
                // chosen to be of type std::set as in the set up of
                // the adjacency, a similar edge might be created
                // multiple times.  And to prevent the definition of
                // parallell edges, we use std::set (=boost::setS)
                // rather than std::vector (=boost::vecS)
                typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
                typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
                BoostGraph boostGraphObj;

                int tempGraphVertId = 0;
                int nVerts;
                int vertCnt;
                int edgeCnt;
                map<int, int>          vertTempGraphVertId;
                map<int, int>          edgeTempGraphVertId;
                Array<OneD, int>       localVerts;
                Array<OneD, int>       localEdges; 
 
                m_numLocalBndCoeffs = 0;
                // First we are going to set up a temporary ordering
                // of the mesh vertices and edges in the graph. We
                // will then later use boost::cuthill_mckee_ordering
                // (reverse Cuthill-McKee algorithm) to minimize the
                // bandwidth.
                //
                // List the (non-dirichlet) vertices and edges of the
                // mesh as the vertices of the temporary graph.  Also
                // define the adjancency between the different
                // vertices of the graph.

                // List the periodic vertices and edges next.  This
                // allows to give corresponding DOF's the same global
                // ID a) periodic vertices
                for(k = 0; k < periodicVerticesId.size(); ++k)
                {
                    for(mapConstIt = periodicVerticesId[k].begin(); mapConstIt != periodicVerticesId[k].end(); mapConstIt++)
                    {
                        meshVertId  = mapConstIt->first;
                        meshVertId2 = mapConstIt->second;
                        

                        // Not entirely sure why we need to have
                        // vertReorderedGraphVertId and
                        // vertTempGraphVertId
                        if(vertReorderedGraphVertId.count(meshVertId) == 0) 
                        {

                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)
                            {
                                
                                if(vertTempGraphVertId.count(meshVertId) == 0)
                                {
                                    vertTempGraphVertId[meshVertId]  = tempGraphVertId;
                                    if(vertTempGraphVertId.count(meshVertId2) == 0)
                                    {
                                        vertTempGraphVertId[meshVertId2] = tempGraphVertId++;
                                    }
                                    else
                                    {
                                        ASSERTL0(false,"Unexplained Periodicity connectivity");
                                    }
                                }
                                else
                                {
                                    if(vertTempGraphVertId.count(meshVertId2) == 0)
                                    {
                                        ASSERTL0(false,"Unexplained Periodicity connectivity");
                                    }
                                    else // Doubly periodic region
                                    {
                                        int id1 = vertTempGraphVertId[meshVertId];
                                        int id2 = vertTempGraphVertId[meshVertId2];
                                        int id;
                                
                                        if(id1 != id2)
                                        {
                                            // Reset any values set to
                                            // id2 to id1. In addition
                                            // if local id is greater
                                            // than id2 decrement list
                                            for(mapIt = vertTempGraphVertId.begin(); 
                                                mapIt != vertTempGraphVertId.end(); mapIt++)
                                            {
                                                id = mapIt->second;
                                                if(id == id2)
                                                {
                                                    vertTempGraphVertId[mapIt->first] = id1;
                                                }
                                                else if (id > id2)
                                                {
                                                    vertTempGraphVertId[mapIt->first] = id-1;
                                                }
                                            }
                                            tempGraphVertId--;
                                        }
                                    }            
                                    

                                }
                                
                            }
                            else
                            {
                                vertReorderedGraphVertId[meshVertId] = vertReorderedGraphVertId[meshVertId2];
                            }
                        }
                        else
                        {
                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)
                            {
                                vertReorderedGraphVertId[meshVertId2] = vertReorderedGraphVertId[meshVertId];
                            }
                        }
                    }
                }

                // b) periodic edges
                for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
                {
                    meshEdgeId  = mapConstIt->first;
                    meshEdgeId2 = mapConstIt->second;

                    if(meshEdgeId < meshEdgeId2)
                    {
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId) == 0,
                                 "This periodic boundary edge has been specified before");
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId2) == 0,
                                 "This periodic boundary edge has been specified before");
                    
                        edgeTempGraphVertId[meshEdgeId]  = tempGraphVertId;
                        edgeTempGraphVertId[meshEdgeId2] = tempGraphVertId++;
                    }
                }

                // List all other vertices and edges
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();

                        nVerts = locExpansion->GetNverts();
                        // For element i, store the temporary graph
                        // vertex id's of all element edges and
                        // verices in these 2 arrays below
                        localVerts = Array<OneD, int>(nVerts,-1);
                        localEdges = Array<OneD, int>(nVerts,-1);
                        vertCnt = 0;
                        edgeCnt = 0;
                        
                        for(j = 0; j < nVerts; ++j) 
                        {   
                            meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {                                
                                if(vertTempGraphVertId.count(meshVertId) == 0)
                                {
                                    vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                                }
                                localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                            }                             
                            
                            meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {                            
                                if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                                {
                                    edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                                }
                                localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            }  
                        }

                        // Now loop over all local edges and vertices
                        // of this element and define that all other
                        // edges and verices of this element are
                        // adjacent to them.
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localVerts[j]==-1)
                            {
                                break;
                            }
                            for(k = 0; k < nVerts; k++)
                            {
                                if(localVerts[k]==-1)
                                {
                                    break;
                                }
                                if(k!=j)
                                {
                                    boost::add_edge( (size_t) localVerts[j], (size_t) localVerts[k],boostGraphObj);
                                }
                            }
                            for(k = 0; k < nVerts; k++)
                            {     
                                if(localEdges[k]==-1)
                                {
                                    break;
                                }    
                                boost::add_edge( (size_t) localVerts[j], (size_t) localEdges[k],boostGraphObj);                       
                            }
                        }
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localEdges[j]==-1)
                            {
                                break;
                            }
                            for(k = 0; k < nVerts; k++)
                            {
                                if(localEdges[k]==-1)
                                {
                                    break;
                                }
                                if(k!=j)
                                {
                                    boost::add_edge( (size_t) localEdges[j], (size_t) localEdges[k],boostGraphObj);   
                                }
                            }
                            for(k = 0; k < nVerts; k++)
                            {                            
                                if(localVerts[k]==-1)
                                {
                                    break;
                                }
                                boost::add_edge( (size_t) localEdges[j], (size_t) localVerts[k],boostGraphObj);                       
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                    }
                }

                // Call boost::cuthill_mckee_ordering to reorder the graph-vertices
                // using the reverse Cuthill-Mckee algorithm
                vector<BoostVertex> inv_perm(tempGraphVertId);
                boost::cuthill_mckee_ordering(boostGraphObj, inv_perm.rbegin());

                // Recast the inverse permutation such that it can be used as a map
                // from old vertex ID to new vertex ID
                Array<OneD, int> iperm(tempGraphVertId);
                for(i = 0; i < tempGraphVertId; i++)
                {
                    iperm[inv_perm[i]]=i;
                }              

                // Fill the vertReorderedGraphVertId and edgeReorderedGraphVertId with the optimal ordering from boost
                for(mapIt = vertTempGraphVertId.begin(); mapIt != vertTempGraphVertId.end(); mapIt++)
                {
                    vertReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId; 
                }
                for(mapIt = edgeTempGraphVertId.begin(); mapIt != edgeTempGraphVertId.end(); mapIt++)
                {
                    edgeReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId; 
                }  
            }
            else
                // Possibility 3: Do not use any optomisation at
                // all. Just list the edges and verices in the order
                // they appear when looping over all the elements in
                // the domain
            {
                m_numLocalBndCoeffs = 0;

                // List the periodic vertices and edges next.
                // This allows to give corresponding DOF's the same
                // global ID
                // a) periodic vertices
                
                for(k = 0; k < periodicVerticesId.size(); ++k)
                {
                    for(mapConstIt = periodicVerticesId[k].begin(); mapConstIt != periodicVerticesId[k].end(); mapConstIt++)
                    {
                        meshVertId  = mapConstIt->first;
                        meshVertId2 = mapConstIt->second;
                        
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                        {
                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId]  = graphVertId;
                                vertReorderedGraphVertId[meshVertId2] = graphVertId++;
                            }
                            else
                            {
                                vertReorderedGraphVertId[meshVertId] = vertReorderedGraphVertId[meshVertId2];
                            }
                        }
                        else
                        {
                            if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId2] = vertReorderedGraphVertId[meshVertId];
                            }
                            else // Doubly periodic vertex case. 
                            {
                                int id1 = vertReorderedGraphVertId[meshVertId];
                                int id2 = vertReorderedGraphVertId[meshVertId2];
                                int id;
                                
                                if(id1 != id2)
                                {
                                    // Reset any values set to id2 to id1. In
                                    // addition if local id is greater than
                                    // id2 decrement list
                                    for(mapIt = vertReorderedGraphVertId.begin(); 
                                        mapIt != vertReorderedGraphVertId.end(); mapIt++)
                                    {
                                        id = mapIt->second;
                                        if(id == id2)
                                        {
                                            vertReorderedGraphVertId[mapIt->first] = id1;
                                        }
                                        else if (id > id2)
                                        {
                                            vertReorderedGraphVertId[mapIt->first] = id-1;
                                        }
                                    }
                                    graphVertId--;
                                }
                            }
                        }
                    }
                }

                // b) periodic edges
                for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
                {
                    meshEdgeId  = mapConstIt->first;
                    meshEdgeId2 = mapConstIt->second;

                    if(meshEdgeId < meshEdgeId2)
                    {
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId) == 0,
                                 "This periodic boundary edge has been specified before");
                        ASSERTL0(edgeReorderedGraphVertId.count(meshEdgeId2) == 0,
                                 "This periodic boundary edge has been specified before");
                    
                        edgeReorderedGraphVertId[meshEdgeId]  = graphVertId;
                        edgeReorderedGraphVertId[meshEdgeId2] = graphVertId++;
                    }
                }

                // List all other vertices and edges
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        for(j = 0; j < locExpansion->GetNverts(); ++j) // number verts = number edges for 2D geom
                        {   
                            meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                            
                            meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {
                                edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                            }                      
                        }
                        m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                    }
                }
            }
#endif //NEKTAR_USING_METIS            
            // Set up an array wich contains the offset information of the different
            // graph vertices. This basically means to identify to how many global degrees
            // of freedom the individual graph vertices correspond. Obviously, the graph vertices
            // corresponding to the mesh-vertices account for a single global DOF. However, the 
            // graph vertices corresponding to the element edges correspond to N-2 global DOF
            // where N is equal to the number of boundary modes on this edge
            Array<OneD, int> graphVertOffset(vertReorderedGraphVertId.size()+
                                             edgeReorderedGraphVertId.size()+1);
            graphVertOffset[0] = 0;
            
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]);
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {  
                    nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                    meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                    meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                    graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]+1] = nEdgeCoeffs-2;
                    graphVertOffset[vertReorderedGraphVertId[meshVertId]+1] = 1;
                    
                    bType = locExpansion->GetEdgeBasisType(j);
                    // need a sign vector for modal expansions if nEdgeCoeffs >=4 
                    if( (nEdgeCoeffs >= 4)&&
                        ( (bType == LibUtilities::eModified_A)||
                          (bType == LibUtilities::eModified_B) ) )
                    {
                        m_signChange = true;
                    }
                }                
            }            
            for(i = 1; i < graphVertOffset.num_elements(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }

            // Allocate the proper amount of space for the class-data 
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numGlobalDirBndCoeffs          = graphVertOffset[firstNonDirGraphVertId];
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int>(nLocBndCondDofs,-1);
            // If required, set up the sign-vector
            if(m_signChange)
            {
                m_localToGlobalSign = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
                m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
                m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD,NekDouble>(nLocBndCondDofs,1.0);
            }

            // Now, all ingredients are ready to set up the actual local to global mapping
            cnt = 0;
            // Loop over all the elements in the domain
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]);
                
                // Loop over all edges (and vertices) of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = (locExpansion->GetGeom2D())->GetEorient(j);                        
                    meshEdgeId          = (locExpansion->GetGeom2D())->GetEid(j);
                    meshVertId          = (locExpansion->GetGeom2D())->GetVid(j)
;
                    
                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                    // Set the global DOF for vertex j of element i
                    m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)] = 
                        graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] =  
                            graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+k;
                    }  
                    
                    // Fill the sign vector if required
                    if(m_signChange)
                    {   
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {                 
                            m_localToGlobalSign[cnt+edgeInteriorMap[k]] = (NekDouble) edgeInteriorSign[k];
                        }
                    }                       
                }                
                cnt += locExpVector[i]->GetNcoeffs();            
            }

            // Set up the mapping for the boundary conditions
            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndSegExp  = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
                
                    for(k = 0; k < 2; k++)
                    {
                        meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndSegExp->GetVertexMap(k)] = 
                            graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    }

                    meshEdgeId = (bndSegExp->GetGeom1D())->GetEid(); 
                    bndEdgeCnt = 0;
                    for(k = 0; k < bndSegExp->GetNcoeffs(); k++)
                    {
                        if(m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] = 
                                graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+bndEdgeCnt;
                            bndEdgeCnt++;
                        }
                    }
                    cnt += bndSegExp->GetNcoeffs();
                }
            }
    
            globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
            m_numGlobalBndCoeffs = globalId;
            
            cnt=0;
            // Setup interior mapping and the boundary map 
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    if(m_signChange)
                    {
                        m_localToGlobalBndSign[cnt]=m_localToGlobalSign[i];
                    }
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }
            m_numGlobalCoeffs = globalId;  
        }

        void LocalToGlobalC0ContMap::SetUp3DExpansionC0ContMap(const int numLocalCoeffs, 
                                                               const StdRegions::StdExpansionVector &locExpVector, 
                                                               const Array<OneD, const ExpList2DSharedPtr> &bndCondExp,
                                                               const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                               const map<int,int>& periodicVerticesId,
                                                               const map<int,int>& periodicEdgesId,
                                                               const map<int,int>& periodicFacesId)
        {
            int i,j,k,l;
            int cnt = 0;
            int intDofCnt;
            int meshVertId;
            int meshVertId2;
            int meshEdgeId;
            int meshEdgeId2;
            int meshFaceId;
            int meshFaceId2;
            int globalId;
            int nEdgeInteriorCoeffs;
            int nFaceInteriorCoeffs;
            int firstNonDirGraphVertId;
            int nLocBndCondDofs = 0;
            int nLocDirBndCondDofs = 0;
            int graphVertId = 0;
            StdRegions::StdExpansion3DSharedPtr locExpansion;  
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            LibUtilities::BasisType             bType;  
            StdRegions::EdgeOrientation         edgeOrient;
            StdRegions::FaceOrientation         faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap; 
            Array<OneD, int>                    edgeInteriorSign;   
            Array<OneD, unsigned int>           faceInteriorMap;  
            Array<OneD, int>                    faceInteriorSign;  

            m_signChange = false;

            // The only unique identifiers of the vertices and edges
            // of the mesh are the vertex id and the mesh id (stored
            // in their corresponding Geometry object).  However,
            // setting up a global numbering based on these id's will
            // not lead to a suitable or optimal numbering. Mainly
            // because: - we want the Dirichlet DOF's to be listed
            // first - we want an optimal global numbering of the
            // remaining DOF's (strategy still need to be defined but
            // can for example be: minimum bandwith or minimum fill-in
            // of the resulting global system matrix)
            //
            // That's why the vertices annd egdes will be
            // rearranged. Currently, this is done in the following
            // way: The vertices and edges of the mesh are considered
            // as vertices of a graph (in a computer science
            // way)(equivalently, they can also be considered as
            // boundary degrees of freedom, whereby all boundary modes
            // of a single edge are considered as a single DOF). We
            // then will use algorithms to reorder these
            // graph-vertices (or boundary DOF's).
            //
            // Two different containers are used to store the graph
            // vertex id's of the different mesh vertices and
            // edges. They are implemented as a STL map such that the
            // graph vertex id can later be retrieved by the unique
            // mesh vertex or edge id's which serve as the key of the
            // map.
            map<int,int> vertReorderedGraphVertId;
            map<int,int> edgeReorderedGraphVertId; 
            map<int,int> faceReorderedGraphVertId; 
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;
           
            // STEP 1: Order the Dirichlet vertices and edges first
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {                               
                    bndCondFaceExp = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
                    if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {          
                        meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid();  
                        faceReorderedGraphVertId[meshFaceId] = graphVertId++;  
                        for(k = 0; k < bndCondFaceExp->GetNverts(); k++)
                        {
                            meshVertId = (bndCondFaceExp->GetGeom2D())->GetVid(k);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                        }

                        for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                        {
                            meshEdgeId = (bndCondFaceExp->GetGeom2D())->GetEid(k);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {
                                edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                            }
                        }
                        nLocDirBndCondDofs += bndCondFaceExp->GetNcoeffs();
                    }
                    nLocBndCondDofs += bndCondFaceExp->GetNcoeffs();
                }
            }
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs;
            firstNonDirGraphVertId = graphVertId;

            // STEP 2: Now order all other vertices and edges in the graph
            // Possibility 3: Do not use any optomisation at all. Just list the edges and verices in the
            // order they appear when looping over all the elements in the domain
         
            m_numLocalBndCoeffs = 0;
            
            // List the periodic vertices and edges next.
            // This allows to give corresponding DOF's the same
            // global ID
            // a) periodic vertices
            for(mapConstIt = periodicVerticesId.begin(); mapConstIt != periodicVerticesId.end(); mapConstIt++)
            {
                meshVertId  = mapConstIt->first;
                meshVertId2 = mapConstIt->second;
                
                if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                {
                    if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                    {
                        vertReorderedGraphVertId[meshVertId]  = graphVertId;
                        vertReorderedGraphVertId[meshVertId2] = graphVertId++;
                    }
                    else
                    {
                        vertReorderedGraphVertId[meshVertId] = vertReorderedGraphVertId[meshVertId2];
                    }
                }
                else
                {
                    if(vertReorderedGraphVertId.count(meshVertId2) == 0)        
                    {
                        vertReorderedGraphVertId[meshVertId2] = vertReorderedGraphVertId[meshVertId];
                    }
                    else
                    {
                        ASSERTL0(vertReorderedGraphVertId[meshVertId2] == vertReorderedGraphVertId[meshVertId],
                                 "These values should be equal");
                    }
                }
            }
            
            // b) periodic edges
            for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
            {
                meshEdgeId  = mapConstIt->first;
                meshEdgeId2 = mapConstIt->second;

                if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                {
                    if(edgeReorderedGraphVertId.count(meshEdgeId2) == 0)        
                    {
                        edgeReorderedGraphVertId[meshEdgeId]  = graphVertId;
                        edgeReorderedGraphVertId[meshEdgeId2] = graphVertId++;
                    }
                    else
                    {
                        edgeReorderedGraphVertId[meshEdgeId] = edgeReorderedGraphVertId[meshEdgeId2];
                    }
                }
                else
                {
                    if(edgeReorderedGraphVertId.count(meshEdgeId2) == 0)        
                    {
                        edgeReorderedGraphVertId[meshEdgeId2] = edgeReorderedGraphVertId[meshEdgeId];
                    }
                    else
                    {
                        ASSERTL0(edgeReorderedGraphVertId[meshEdgeId2] == edgeReorderedGraphVertId[meshEdgeId],
                                 "These values should be equal");
                    }
                }
            }
            
            // b) periodic faces
            for(mapConstIt = periodicFacesId.begin(); mapConstIt != periodicFacesId.end(); mapConstIt++)
            {
                meshFaceId  = mapConstIt->first;
                meshFaceId2 = mapConstIt->second;

                if(meshFaceId < meshFaceId2)
                {
                    ASSERTL0(faceReorderedGraphVertId.count(meshFaceId) == 0,
                             "This periodic boundary face has been specified before");
                    ASSERTL0(faceReorderedGraphVertId.count(meshFaceId2) == 0,
                             "This periodic boundary face has been specified before");
                    
                    faceReorderedGraphVertId[meshFaceId]  = graphVertId;
                    faceReorderedGraphVertId[meshFaceId2] = graphVertId++;
                }
            }

            // List all other vertices and edges and faces
            for(i = 0; i < locExpVector.size(); ++i)
            { 
                if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[i]))
                {
                    for(j = 0; j < locExpansion->GetNverts(); ++j)
                    {   
                        meshVertId = (locExpansion->GetGeom3D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                        {
                            vertReorderedGraphVertId[meshVertId] = graphVertId++;
                        }
                    }
                    
                    for(j = 0; j < locExpansion->GetNedges(); ++j)
                    {     
                        meshEdgeId = (locExpansion->GetGeom3D())->GetEid(j);
                        if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                        {
                            edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                        }                      
                    }
                    
                    for(j = 0; j < locExpansion->GetNfaces(); ++j)
                    {     
                        meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);
                        if(faceReorderedGraphVertId.count(meshFaceId) == 0)        
                        {
                            faceReorderedGraphVertId[meshFaceId] = graphVertId++;
                        }                      
                    }
                    m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
            }     
            // Set up an array wich contains the offset information of the different
            // graph vertices. This basically means to identify to how many global degrees
            // of freedom the individual graph vertices correspond. Obviously, the graph vertices
            // corresponding to the mesh-vertices account for a single global DOF. However, the 
            // graph vertices corresponding to the element edges correspond to N-2 global DOF
            // where N is equal to the number of boundary modes on this edge
            Array<OneD, int> graphVertOffset(vertReorderedGraphVertId.size()+
                                             edgeReorderedGraphVertId.size()+
                                             faceReorderedGraphVertId.size()+1);
            graphVertOffset[0] = 0;
            
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[i]);

                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {  
                    meshVertId = (locExpansion->GetGeom3D())->GetVid(j);
                    graphVertOffset[vertReorderedGraphVertId[meshVertId]+1] = 1;
                }    

                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {  
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j) - 2;
                    meshEdgeId = (locExpansion->GetGeom3D())->GetEid(j);
                    graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]+1] = nEdgeInteriorCoeffs;
                    
                    bType = locExpansion->GetEdgeBasisType(j);
                    // need a sign vector for modal expansions if nEdgeCoeffs >=4 
                    if( (nEdgeInteriorCoeffs+2 >= 4)&&
                        ( (bType == LibUtilities::eModified_A)||
                          (bType == LibUtilities::eModified_B) ) )
                    {
                        m_signChange = true;
                    }
                }      

                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {  
                    nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);
                    graphVertOffset[faceReorderedGraphVertId[meshFaceId]+1] = nFaceInteriorCoeffs;
                }               
            }            
            for(i = 1; i < graphVertOffset.num_elements(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }

            // Allocate the proper amount of space for the class-data 
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numGlobalDirBndCoeffs          = graphVertOffset[firstNonDirGraphVertId];
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int>(nLocBndCondDofs,-1);
            // If required, set up the sign-vector
            if(m_signChange)
            {
                m_localToGlobalSign = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
                m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
                m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD,NekDouble>(nLocBndCondDofs,1.0);
            }

            // Now, all ingredients are ready to set up the actual local to global mapping
            cnt = 0;
            // Loop over all the elements in the domain
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[i]);
                
                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {
                    meshVertId          = (locExpansion->GetGeom3D())->GetVid(j);

                    // Set the global DOF for vertex j of element i
                    m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)] = 
                        graphVertOffset[vertReorderedGraphVertId[meshVertId]];                      
                }      

                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = (locExpansion->GetGeom3D())->GetEorient(j);                        
                    meshEdgeId          = (locExpansion->GetGeom3D())->GetEid(j);
                    
                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                    
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] =  
                            graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+k;
                    }  
                    
                    // Fill the sign vector if required
                    if(m_signChange)
                    {   
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {                 
                            m_localToGlobalSign[cnt+edgeInteriorMap[k]] = (NekDouble) edgeInteriorSign[k];
                        }
                    }                       
                }     

                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    faceOrient          = (locExpansion->GetGeom3D())->GetFaceorient(j);                        
                    meshFaceId          = (locExpansion->GetGeom3D())->GetFid(j);
                    
                    locExpansion->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);
                    
                    // Set the global DOF's for the interior modes of face j
                    for(k = 0; k < nFaceInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+faceInteriorMap[k]] =  
                            graphVertOffset[faceReorderedGraphVertId[meshFaceId]]+k;
                    } 
           
                    if(m_signChange)
                    {
                        for(k = 0; k < nFaceInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+faceInteriorMap[k]] = (NekDouble) faceInteriorSign[k];
                        } 
                    }                    
                }   
              
                cnt += locExpVector[i]->GetNcoeffs();            
            }

            // Set up the mapping for the boundary conditions
            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndCondFaceExp  = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
                
                    for(k = 0; k < bndCondFaceExp->GetNverts(); k++)
                    {
                        meshVertId = (bndCondFaceExp->GetGeom2D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndCondFaceExp->GetVertexMap(k)] = 
                            graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    }

                    for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                    {
                        nEdgeInteriorCoeffs = bndCondFaceExp->GetEdgeNcoeffs(k)-2;
                        edgeOrient          = (bndCondFaceExp->GetGeom2D())->GetEorient(k); 
                        meshEdgeId          = (bndCondFaceExp->GetGeom2D())->GetEid(k); 

                        bndCondFaceExp->GetEdgeInteriorMap(k,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                        for(l = 0; l < nEdgeInteriorCoeffs; ++l)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+edgeInteriorMap[l]] =  
                                graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+l;
                        }  
                    
                        // Fill the sign vector if required 
                        if(m_signChange)
                        {   
                            for(l = 0; l < nEdgeInteriorCoeffs; ++l)
                            {                 
                                m_bndCondCoeffsToGlobalCoeffsSign[cnt+edgeInteriorMap[l]] = (NekDouble) edgeInteriorSign[l];
                            }
                        }  
                    }

                    meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid(); 
                    intDofCnt = 0;
                    for(k = 0; k < bndCondFaceExp->GetNcoeffs(); k++)
                    {
                        if(m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] = 
                                graphVertOffset[faceReorderedGraphVertId[meshFaceId]]+intDofCnt;
                            intDofCnt++;
                        }
                    }
                    cnt += bndCondFaceExp->GetNcoeffs();
                }
            }
    
            globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
            m_numGlobalBndCoeffs = globalId;
            
            cnt=0;
            // Setup interior mapping and the boundary map 
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    if(m_signChange)
                    {
                        m_localToGlobalBndSign[cnt]=m_localToGlobalSign[i];
                    }
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }
            m_numGlobalCoeffs = globalId;  
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
        void LocalToGlobalC0ContMap::CalculateBndSystemBandWidth(const StdRegions::StdExpansionVector &locExpVector)
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locSize = locExpVector[i]->NumBndryCoeffs();
                maxId = -1;
                minId = m_numLocalCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalBndMap[cnt+j] >= m_numGlobalDirBndCoeffs)
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
        
        void LocalToGlobalC0ContMap::CalculateFullSystemBandWidth(const StdRegions::StdExpansionVector &locExpVector)
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locSize = locExpVector[i]->GetNcoeffs();
                maxId = -1;
                minId = m_numLocalCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalMap[cnt+j] >= m_numGlobalDirBndCoeffs)
                    {
                        if(m_localToGlobalMap[cnt+j] > maxId)
                        {
                            maxId = m_localToGlobalMap[cnt+j];
                        }
                        
                        if(m_localToGlobalMap[cnt+j] < minId)
                        {
                            minId = m_localToGlobalMap[cnt+j];
                        }
                    }
                }
                bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

                cnt+=locSize;
            }

            m_fullSystemBandWidth = bwidth;
        }

    }
}

/**
 * $Log: LocalToGlobalC0ContMap.cpp,v $
 * Revision 1.10  2009/04/27 21:34:58  sherwin
 * Modified WriteToField Method
 *
 * Revision 1.9  2009/04/20 16:14:06  sherwin
 * Updates for optimising bandwidth of DG solver and allowing write import on explist
 *
 * Revision 1.8  2009/04/03 15:06:35  sherwin
 * Added LinearAdvection, Generalised GenGlobalLinSys and added in bandwidth routine for DG solver
 *
 * Revision 1.7  2009/04/02 13:06:42  sherwin
 * Modified to take symmetric banded system for HDH solver
 *
 * Revision 1.6  2009/02/27 15:47:13  sherwin
 * Corrected bug in initialising just Neumann boundary conditions
 *
 * Revision 1.5  2008/12/19 15:12:39  pvos
 * Updates for precomputed dirichlet forcing functionality
 *
 * Revision 1.4  2008/11/05 16:15:24  pvos
 * Added bandwith calculation routines
 *
 * Revision 1.3  2008/09/23 18:21:00  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.2  2008/09/17 13:46:40  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.1  2008/09/16 13:36:06  pvos
 * Restructured the LocalToGlobalMap classes
 *
 **/
