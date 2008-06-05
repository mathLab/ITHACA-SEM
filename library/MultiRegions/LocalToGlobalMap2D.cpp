///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap2D.cpp
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
// Description: Local to Global mapping routines in 2D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalMap2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * \page pageConnectivity Connectivity in Nektar++
         *
         * \section Connectivity Connectivity
         * The typical elemental decomposition of the spectral/hp element
         * method requires a global assembly process when considering
         * multi-elemental problems. This global assembly will ensure some
         * level of connectivity between adjacent elements sucht that there is
         * some form of continuity across element boundaries in the global solution. 
         * In this section, we will merely focus on the classical Galerkin 
         * method, where global continuity is typically imposed by making the
         * approximation \f$C^0\f$ continuous.
         * 
         * \subsection Connectivity2D Connectivity in two dimensions
         * As explained in [1], the global assembly process involves
         * the transformation from \a local \a degrees \a of \a freedom to \a global \a degrees
         * \a of \a freedom (DOF). This transformation is typically done by a mapping array
         * which relates the numbering of the local (= elemental) DOF's to the
         * numbering of the global DOF's. To understand how this transformation is set up
         * in Nektar++ one should understand the following:
         * - <b> Starting point </b><BR>
         *   The starting point is the initial numbering of the elemental expansion modes.
         *   This corresponds to the order in which the different local expansion modes are
         *   listed in the coefficient array \a m_coeffs of the elemental (local or standard) expansion. 
         *   The specific order in which the different elemental expansion modes appear is motivated by 
         *   the compatability with the sum-factorisation technique. This also implies that this 
         *   ordering is fixed and should not be changed by the user. Hence, this unchangeable initial 
         *   local numbering will serve as starting input for the connectivity. 
         * - <b> end point </b><BR>
         *   Obviously, we are working towards the numbering of the global DOF's. This global ordering
         *   should:
         *   -  reflect the chosen continuity approach (standard \f$C^0\f$ Galerkin in our case)
         *   -  (optionally) have some optimal ordering (where optimality can be defined in different ways, e.g.
         *      minimal bandwith)
         *
         * All intermittent steps from starting point to end point can basically be chosen freely but they should allow 
         * for an efficient construction of the global numbering system starting from the elemental ordering
         * of the local degrees of freedom. Currently, Nektar++ provides a number of tools and routines
         * in the different sublibraries which can be employed to set up the mapping from local to global DOF's.
         * These tools  will be listed below, but first the connectivity strategies for both modal and nodal expansions
         * will be explained. Note that all explanations below are focussed on quadrilateral elements. However, the general 
         * idea equally holds for triangles.
         *
         * \subsubsection Connectivitystrategies Connectivity strategies
         * For a better understanding of the described strategies, one should first understand
         * how Nektar++ deals with the basic geometric concepts such as edges and (2D) elements. 
         *
         * In Nektar++, A (2D) element is typically defined by a set of edges. In the input .xml files, this should be done as
         * (for a quadrilateral element):
         * \code
         * <Q ID="i"> e0 e1 e2 e3</Q>
         * \endcode
         * where \a e0 to \a e3 correspond to the mesh ID's of the different edges. It is important to know that
         * in Nektar++, the convention is that these edges should be ordered counterclokwise (Note that this order also
         * corresponds to the order in which the edges are passed to the constructors of the 2D geometries). In addition, note that
         * we will refer to edge \a e0 as the edge with local (elemental) edge ID equal to 0, to edge \a e1 as local edge with ID
         * 1, to edge \a e2 as local edge with ID equal 2 and to edge \a e3 as local edge with ID 3.
         * Furthermore, one should note that the local coordinate system is orientated such that the first coordinate axis is
         * aligned with edge 0 and 2 (local edge ID), and the second coordinate axis is aligned with edge 1 and 3. 
         * The direction of these coordinate axis is such that it points in counterclockwise direction for edges 0 and 1, and
         * in clockwise direction for edge 2 and 3.
         *
         * Another important feature in the connectivity strategy is the concept of edge orientation.
         * For a better understanding, consider the input format of an edge as used in the input .xml files which contain
         * the information about the mesh. An edge is defined as:
         * \code
         * <E ID="i"> v0 v1 </E>
         * \endcode
         * where \a v0 and \a v1 are the ID's of the two vertices that define the edge (Note that these
         * vertices are passed in the same order to the constructor of the edge). Now, 
         * the orientation of an edge of a two-dimensional element (i.e. quadrilateral or triangle) is 
         * defined as:
         * - <b> Forward </b> if the vertex with ID \a v0 comes before the vertex with ID \a v1
         *   when considering the vertices the vertices of the element in a counterclockwise direction
         * - <b> Backward </b> otherwise.<BR>
         * 
         *
         * This has the following implications:
         * - The common edge of two adjacent elements has always a forward orientation for one of the
         *   elements it belongs to and a backward orientation for the other.
         * - The orientation of an edge is only relevant when considering two-dimensional elements.
         *   It is a property which is not only inherent to the edge itself, but depends on the
         *   element it belongs to. (This also means that a segment does not have an orientation) 
         *
         * \paragraph Modalexpansions Modal expansions
         * We will follow the basic principles of the connectivity strategy as explained in Section 4.2.1.1 of [1] (such as the
         * hierarchic ordering of the edge modes).
         * However, we do not follow the strategy described to negate the odd modes of an intersecting edge of 
         * two adjacent elements if the local coordinate systems have an opposite direction. The explained strategy 
         * involves checking the direction of the local coordinate systems of the neighbouring elements.
         * However, for a simpler automatic procedure to identify which edges need to have odd mode negated, we 
         * would like to have an approach which can be applied to the elements individually, without information 
         * being coupled between neighbouring elements.<BR>
         * This can be accomplished in the following way:<BR>
         * Note that this approach is based on the earlier observation that the intersecting edge of two elements
         * always has opposite orientation. Proper connectivity can now be guaranteed if:
         * - forward oriented edges always have a counterclockwise local coordinate axis
         * - backward oriented edges always have a clockwise local coordinate axis.
         * 
         * Both the local coordinate axis along an intersecting edge will then point in the same
         * direction.
         * Obviously, these conditions will not be fulfilled by default. But in order to do so,
         * the direction of the local coordinate axis should be reversed in following situations:
         * \code
         * if ((LocalEdgeId == 0)||(LocalEdgeId == 1))
         * {
         *     if( EdgeOrientation == Backward )
         *     {
         *         change orientation of local coordinate axis
         *     }
         * }
         *
         * if ((LocalEdgeId == 2)||(LocalEdgeId == 3))
         * {
         *     if( EdgeOrientation == Forward )
         *     {
         *         change orientation of local coordinate axis
         *     }
         * }
         * \endcode
         * This algorithm above is based on the earlier observation that the local coordinate axis automatically point in
         * counterclockwise direction for edges 0 and 1 and in clockwise direction for the other edges.
         * As explained in [1] the change in local coordinate axis can actually be done by reversing the sign
         * of the odd modes. This is implemented by means of an additional sign vector.
         * \paragraph Nodalexpansions Nodal expansions
         * For the nodal expansions, we will use the connectivity strategy as explained in Section 4.2.1.1 of [1]. 
         * However, we will clarify this strategy from a Nektar++ point of view. As pointed out in [1], the nodal
         * edge modes can be identified with the physical location of the nodal points. In order to ensure proper connectivity
         * between elements the egde modes with the same nodal location should be matched. This will be accomplished if both the 
         * sets of local edge modes along the intersection edge of two elements are numbered in the same direction. And as the 
         * intersecting edge of two elements always has opposite direction, this can be guaranteed if:
         * - the local numbering of the edge modes is counterclockwise for forward oriented edges
         * - the local numbering of the edge modes is clockwise for backward oriented edges.
         *
         * This will ensure that the numbering of the global DOF's on an edge is in the same direction as the tow subsets of
         * local DOF's on the intersecting edge.
         *
         * \subsubsection Implementation Implementation
         * The main entity for the transformation from local DOF's to global DOF's (in 2D) is the LocalToGlobalMap2D class.
         * This class basically is the abstraction of the mapping array <tt>map[e][i]</tt> as introduced in section 4.2.1 of [1]. This mapping
         * array is contained in the class' main data member, LocalToGlobalMap2D::m_locToContMap. Let us recall what this mapping array <BR><BR>
         *      <tt>map[e][i] = globalID</tt><BR><BR>
         * actually represents:
         * - <tt>e</tt> corresponds to the \a e th element 
         * - <tt>i</tt> corresponds to the \a i th expansion mode within element \a e. This index \a i in this map array corresponds to 
         *   the index of the coefficient array m_coeffs.
         * - <tt>globalID</tt> represents the ID of the corresponding global degree of freedom. 
         *
         * However, rather than this two-dimensional structure of the mapping array, LocalToGlobalMap2D::m_locToContMap 
         * stores the mapping array as a one-dimensional array which is the concatenation of the different elemental mapping
         * arrays <tt>map[e]</tt>. This mapping array can then be used to assemble the global system out of the local entries, or to
         * do any other transformation between local and global degrees of freedom (Note that other mapping arrays such as the boundary mapping 
         * <tt>bmap[e][i]</tt> or the sign vector <tt>sign[e][i]</tt> which might be required are also contained in the LocalToGlobalMap2D class).
         *
         * For a better appreciation of the implementation of the connectivity in Nektar++, it might be useful to consider how
         * this mapping array is actually being constructed (or filled). To understand this, first consider the following:
         * - The initial local elemental numbering (referred to as starting point in the beginning of this document) is not suitable
         *   to set up the mapping. In no way does it correspond to the local numbering required for a proper connectivity as elaborated in
         *   Section 4.2.1 of [1]. Hence, this initial ordering asuch cannot be used to implement the connectivity strategies explained above.
         *   As a result, additional routines (see \ref GetVertexMap "here"), which account for some kind of reordering of the local numbering will be required in order to construct
         *   the mapping array properly.
         * - Although the different edge modes can be thought of as to include both the vertex mode, we will make a clear distinction between
         *   them in the implementation. In other words, the vertex modes will be treated separately from the other modes on
         *   an edge as they are conceptually different from an connectivity point of view.
         *   We will refer to these remaining modes as \a interior \a edge \a modes. 
         *
         * The fill-in of the mapping array can than be summarised by the following part of (simplified) code:
         * \code        
         * for(e = 0; e < Number_Of_2D_Elements; e++)
         * {
         *     for(i = 0; i < Number_Of_Vertices_Of_Element_e; i++)
         *     {
         *         offsetValue = ...
         *         map[e][GetVertexMap(i)] = offsetValue;
         *     }
         *
         *     for(i = 0; i < Number_Of_Edges_Of_Element_e; i++)
         *     {
         *         localNumbering = GetEdgeInteriorMap(i);
         *         offsetValue = ...
         *         for(j = 0; j < Number_Of_InteriorEdgeModes_Of_Edge_i; j++)
         *         {
         *             map[e][localNumbering(j)] = offsetValue + j;
         *         }
         *     }
         * }
         * \endcode
         * In this document, we will not cover how the calculate the offsetValue which:
         * - for the vertices, corresponds to the global ID of the specific vertex
         * - for the edges, corresponds to the starting value of the global numbering on the concerned edge
         *
         * However, we would like to focus on the routines GetVertexMap() and GetEdgeInteriorMap(), 2 functions
         * which somehow reorder the initial local numbering in order to be compatible with the connectivity 
         * strategy. 
         *
         * \paragraph GetVertexMap GetVertexMap()
         * The function StdRegions::StdExpansion::GetVertexMap (which is defined for all two-dimensional expansions in the StdRegions::StdExpansion 
         * class tree) is a very straightforward routine which simply does the following
         * - Given the local vertex id (i.e. 0,1,2 or 3 for a quadrilateral element), it returns the position
         *   of the corresponding vertex mode in the elemental coefficient array StdRegions::StdExpansion::m_coeffs. By using this function
         *   as in the code above, it is ensured that the global ID of the vertex is entered in the correct position in the 
         *   mapping array.
         * \paragraph GetEdgeInteriorMap GetEdgeInteriorMap()
         * Like the previous routine, this function is also defined for all two dimensional expanions in the
         * StdRegions::StdExpansion class tree. This is actually the most important function to ensure proper connectivity 
         * between neigbouring elements. It is a function which reorders the numbering of local DOF's according
         * to the connectivity strategy.
         * As input this function takes:
         * - the local edge ID of the edge to be considered
         * - the orientation of this edge.
         *
         * As output, it returns:
         * - the local ordering of the requested (interior) edge modes. This is contained in an array of size \a N-2, 
         *   where N is the number of expansion modes in the relevant direction. The entries in this array represent
         *   the position of the corresponding interior edge mode in the elemental coefficient array StdRegions::StdExpansion::m_coeffs.
         *
         * Rather than the actual values of the local numbering, it is the ordering of local edge modes which is of importance
         * for the connectivity. That is why it is important how the different interior edge modes are sorted in the returned
         * array. This should be such that for both the elements which are connected by the intersecting edge, the local (interior)
         * edge modes are iterated in the same order. This will guarantee a correct global numbering scheme when employing
         * the algorithm shown above. This proper connectivity can be ensured if the function GetEdgeInteriorMap:
         * - for modal expansions: returns the edge interior modes in hierarchical order (i.e. the lowest polynomial order mode first),
         * - for nodal expansions: returns the edge interior modes in:
         *   - counterclockwise order for forward oriented edges
         *   - clockwise order for backward oriented edges.
         * \subsection References References
         * [1] G.E. Karniadakis \& S.J. Sherwin: Spectral/hp Element Methods for 
         * Computational Fluid Dynamics, Oxford Science Publications, 2005
         */
        LocalToGlobalMap2D::LocalToGlobalMap2D()
        {
        }
        
        LocalToGlobalMap2D::LocalToGlobalMap2D(const int loclen, 
                                               const StdRegions::StdExpansionVector &locExpVector, 
                                               const Array<OneD, const MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                               const Array<OneD, const SpatialDomains::BoundaryConditionType> &bndCondTypes)
        {
            int i,j,k;
            int cnt = 0;
            int bndEdgeCnt;
            int meshVertId;
            int meshEdgeId;
            int globalId;
            int nEdgeCoeffs;
            int nEdgeInteriorCoeffs;
            int firstNonDirGraphVertId;
            int nLocBndCondDofs = 0;
            int graphVertId = 0;
            StdRegions::StdExpansion2DSharedPtr locExpansion;  
            LocalRegions::SegExpSharedPtr       bndSegExp;
            LibUtilities::BasisType             bType;  
            StdRegions::EdgeOrientation         edgeOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;  

            // The only unique identifiers of the vertices and edges of the mesh are the
            // vertex id and the mesh id (stored in their corresponding Geometry object).
            // However, setting up a global numbering based on these id's will not lead to 
            // a suitable or optimal numbering. Mainly because:
            // - we want the Dirichlet DOF's to be listed first
            // - we want an optimal global numbering of the remaining DOF's (strategy still
            //   need to be defined but can for example be: minimum bandwith or minimum fill-in
            //   of the resulting global system matrix)
            //
            // That's why the vertices annd egdes will be rearranged. Currently, this is done 
            // in the following way: The vertices and edges of the mesh are considered as vertices 
            // of a graph (in a computer science way)(equivalently, they can also be considered as 
            // boundary degrees of freedom, whereby all boundary modes of a single edge are considered
            // as a single DOF). We then will use algorithms to reorder these 
            // graph-vertices (or boundary DOF's).
            //
            // Two different containers are used to store the graph vertex id's of the different
            // mesh vertices and edges. They are implemented as a STL map such that the graph vertex id
            // can later be retrieved by the unique mesh vertex or edge id's which serve as the key of the map.  
            map<int, int> vertReorderedGraphVertId;
            map<int, int> edgeReorderedGraphVertId; 
           
            // STEP 1: Order the Dirichlet vertices and edges first
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {                               
                    if(bndCondTypes[i]==SpatialDomains::eDirichlet)
                    {          
                        bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
                        meshEdgeId = (bndSegExp->GetGeom())->GetEid();  
                        edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;  
                        for(k = 0; k < 2; k++)
                        {
                            meshVertId = (bndSegExp->GetGeom())->GetVid(k);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                        }
                    }
                    nLocBndCondDofs += bndSegExp->GetNcoeffs();
                }
            }
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
                map<int,int>::iterator mapIt;
                Array<OneD, int>       localVerts;
                Array<OneD, int>       localEdges;               

                m_totLocBndDofs = 0;
                // First we are going to set up a temporary ordering of the mesh vertices and edges
                // in the graph. We will then later use METIS to optimise this ordering
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        m_totLocBndDofs += locExpansion->NumBndryCoeffs();
                        
                        nVerts = locExpansion->GetNverts();
                        // For element i, store the temporary graph vertex id's of all element
                        // edges and verices in these 2 arrays below
                        localVerts = Array<OneD, int>(nVerts,-1);
                        localEdges = Array<OneD, int>(nVerts,-1);
                        vertCnt = 0;
                        edgeCnt = 0;
                        
                        // List the (non-dirichlet) vertices and edges of the mesh as the vertices of the temporary graph
                        //
                        // At the same time, we can already start setting ip the adjacency information of the graph. 
                        // This is required as input of the METIS routines. This adjacency information should list 
                        // all adjacent graph-vertices for all the separate graph-vertices
                        for(j = 0; j < nVerts; ++j) 
                        {   
                            meshVertId = (locExpansion->GetGeom())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {                                
                                if(vertTempGraphVertId.count(meshVertId) == 0)
                                {
                                    vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                                }
                                localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                            }                             
                            
                            meshEdgeId = (locExpansion->GetGeom())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {                            
                                if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                                {
                                    edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                                }
                                localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            }  
                        }

                        // Now loop over all local edges and vertices of this element and define 
                        // that all other edges and verices of this element are adjacent to them.
                        // To do so, we use an STL map where the key is the unique temporary graph-vertex id
                        // and the mapped values is an STL set which contains the temporary graph-vertex id's
                        // of all adajacent graph vertices. By looping over all elements in the mesh and everytime
                        // appending the adjacent grapgh vertices, we will make sure that all the adjacency information
                        // will be complete.
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

                // Now, translate this adjacency information to a format to the format that METIS 
                // understands.
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

                // Fill the vertReorderedGraphVertId and edgeReorderedGraphVertId with the optimal ordering from METIS
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
            // Possibility 2: Minimize the bandwith using a reversed Cuthill-McKee algorithm 
            // (as provided by the Boost Graph Library)
            if(false)
            {
                // the first template parameter (=OutEdgeList) is chosen to be of type std::set
                // as in the set up of the adjacency, a similar edge might be created multiple times.
                // And to prevent the definition of paralelle edges, we use std::set (=boost::setS)
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
                map<int,int>::iterator mapIt;
                Array<OneD, int>       localVerts;
                Array<OneD, int>       localEdges; 
 
                m_totLocBndDofs = 0;
                // First we are going to set up a temporary ordering of the mesh vertices and edges
                // in the graph. We will then later use boost::cuthill_mckee_ordering (reverse Cuthill-McKee algorithm)  
                // to minimize the bandwidth.
                //
                // List the (non-dirichlet) vertices and edges of the mesh as the vertices of the temporary graph.
                // Also define the adjancency between the different vertices of the graph.
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        m_totLocBndDofs += locExpansion->NumBndryCoeffs();

                        nVerts = locExpansion->GetNverts();
                        // For element i, store the temporary graph vertex id's of all element
                        // edges and verices in these 2 arrays below
                        localVerts = Array<OneD, int>(nVerts,-1);
                        localEdges = Array<OneD, int>(nVerts,-1);
                        vertCnt = 0;
                        edgeCnt = 0;
                        
                        for(j = 0; j < nVerts; ++j) 
                        {   
                            meshVertId = (locExpansion->GetGeom())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {                                
                                if(vertTempGraphVertId.count(meshVertId) == 0)
                                {
                                    vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                                }
                                localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                            }                             
                            
                            meshEdgeId = (locExpansion->GetGeom())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {                            
                                if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                                {
                                    edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                                }
                                localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            }  
                        }

                        // Now loop over all local edges and vertices of this element and define 
                        // that all other edges and verices of this element are adjacent to them.
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
            // Possibility 3: Do not use any optomisation at all. Just list the edges and verices in the
            // order they appear when looping over all the elements in the domain
            {
                m_totLocBndDofs = 0;
                for(i = 0; i < locExpVector.size(); ++i)
                { 
                    if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]))
                    {
                        for(j = 0; j < locExpansion->GetNverts(); ++j) // number verts = number edges for 2D geom
                        {   
                            meshVertId = (locExpansion->GetGeom())->GetVid(j);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)        
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                            
                            meshEdgeId = (locExpansion->GetGeom())->GetEid(j);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)        
                            {
                                edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                            }                      
                        }
                        m_totLocBndDofs += locExpansion->NumBndryCoeffs();
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
                    meshEdgeId = (locExpansion->GetGeom())->GetEid(j);
                    meshVertId = (locExpansion->GetGeom())->GetVid(j);
                    graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]+1] = nEdgeCoeffs-2;
                    graphVertOffset[vertReorderedGraphVertId[meshVertId]+1] = 1;
                    
                    bType = locExpansion->GetEdgeBasisType(j);
                    // need a sign vector for modal expansions if nEdgeCoeffs >=4 
                    if( (nEdgeCoeffs >= 4)&&
                        ( (bType == LibUtilities::eModified_A)||
                          (bType == LibUtilities::eModified_B) ) )
                    {
                        m_sign_change = true;
                    }
                }                
            }            
            for(i = 1; i < graphVertOffset.num_elements(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }

            // Allocate the proper amount of space for the class-data 
            m_totLocDofs          = loclen;
            m_numDirichletDofs    = graphVertOffset[firstNonDirGraphVertId];
            m_locToContMap        = Array<OneD, int>(m_totLocDofs,-1);
            m_locToContBndMap     = Array<OneD, int>(m_totLocBndDofs,-1);
            m_locToContBndCondMap = Array<OneD,int>(nLocBndCondDofs,-1);
            // If required, set up the sign-vector
            if(m_sign_change)
            {
                m_sign = Array<OneD, NekDouble>(m_totLocDofs,1.0);
                m_bndSign = Array<OneD, NekDouble>(m_totLocBndDofs,1.0);
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
                    edgeOrient          = (locExpansion->GetGeom())->GetEorient(j);                        
                    meshEdgeId          = (locExpansion->GetGeom())->GetEid(j);
                    meshVertId          = (locExpansion->GetGeom())->GetVid(j);
                    
                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap);

                    // Set the global DOF for vertex j of element i
                    m_locToContMap[cnt+locExpansion->GetVertexMap(j)] = 
                        graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_locToContMap[cnt+edgeInteriorMap[k]] =  
                            graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+k;
                    }  
                    
                    // Fill the sign vector if required
                    if(m_sign_change)
                    {                       
                        if( ((j<2)&&(edgeOrient==StdRegions::eBackwards)) ||
                            ((j>1)&&(edgeOrient==StdRegions::eForwards)) )
                        {
                            for(k = 1; k < nEdgeInteriorCoeffs; k+=2)
                            {
                                m_sign[cnt+edgeInteriorMap[k]] = -1;
                            }
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
                        meshVertId = (bndSegExp->GetGeom())->GetVid(k);
                        m_locToContBndCondMap[cnt+bndSegExp->GetVertexMap(k)] = 
                            graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    }

                    meshEdgeId = (bndSegExp->GetGeom())->GetEid(); 
                    bndEdgeCnt = 0;
                    for(k = 0; k < bndSegExp->GetNcoeffs(); k++)
                    {
                        if(m_locToContBndCondMap[cnt+k] == -1)
                        {
                            m_locToContBndCondMap[cnt+k] = 
                                graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+bndEdgeCnt;
                            bndEdgeCnt++;
                        }
                    }
                    cnt += bndSegExp->GetNcoeffs();
                }
            }
    
            globalId = Vmath::Vmax(m_totLocDofs,&m_locToContMap[0],1)+1;
            m_totGloBndDofs = globalId;
            
            cnt=0;
            // Setup interior mapping and the boundary map 
            for(i = 0; i < m_totLocDofs; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = globalId++;
                }
                else
                {
                    if(m_sign_change)
                    {
                        m_bndSign[cnt]=m_sign[i];
                    }
                    m_locToContBndMap[cnt++]=m_locToContMap[i];
                }
            }
            m_totGloDofs = globalId;  

//             // ----------------------------------------------------------------------------
//             // Calculation of the bandwith ----
//             // The bandwidth here calculated corresponds to what is referred to as half-bandwidth.
//             // If the elements of the matrix are designated as a_ij, it corresponds to
//             // the maximum value of |i-j| for non-zero a_ij.
//             int locSize;
//             int maxId;
//             int minId;
//             int bwidth = -1;
//             cnt=0;
//             for(i = 0; i < locExpVector.size(); ++i)
//             {
//                 locSize = locExpVector[i]->NumBndryCoeffs();
//                 maxId = -1;
//                 minId = m_totLocDofs+1;
//                 for(j = 0; j < locSize; j++)
//                 {
//                     if(m_locToContBndMap[cnt+j] >= m_numDirichletDofs)
//                     {
//                         if(m_locToContBndMap[cnt+j] > maxId)
//                         {
//                             maxId = m_locToContBndMap[cnt+j];
//                         }
                        
//                         if(m_locToContBndMap[cnt+j] < minId)
//                         {
//                             minId = m_locToContBndMap[cnt+j];
//                         }
//                     }
//                 }
//                 bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

//                 cnt+=locSize;
//             }
//             cout<<"MatrixSize : "<<m_totGloBndDofs-m_numDirichletDofs<<endl;
//             cout<<"BandWith   : "<<bwidth<<endl;
//             // ----------------------------------------------------------------------------
        }    
        
        LocalToGlobalMap2D::~LocalToGlobalMap2D()
        {
        } 
    }
}

/**
* $Log: LocalToGlobalMap2D.cpp,v $
* Revision 1.13  2008/05/07 16:05:55  pvos
* Mapping + Manager updates
*
* Revision 1.12  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.11  2008/04/02 22:19:54  pvos
* Update for 2D local to global mapping
*
* Revision 1.10  2008/03/18 14:14:13  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.9  2008/01/23 21:50:52  sherwin
* Update from EdgeComponents to SegGeoms
*
* Revision 1.8  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.7  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.6  2007/07/20 02:04:13  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.5  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
