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
 * All intermittent steps from starting point to end point can
 * basically be chosen freely but they should allow for an
 * efficient construction of the global numbering system
 * starting from the elemental ordering of the local degrees
 * of freedom. Currently, Nektar++ provides a number of tools
 * and routines in the different sublibraries which can be
 * employed to set up the mapping from local to global DOF's.
 * These tools will be listed below, but first the
 * connectivity strategies for both modal and nodal expansions
 * will be explained. Note that all explanations below are
 * focussed on quadrilateral elements. However, the general
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
 * Another important feature in the connectivity strategy is
 * the concept of edge orientation.  For a better
 * understanding, consider the input format of an edge as used
 * in the input .xml files which contain the information about
 * the mesh. An edge is defined as: \code <E ID="i"> v0 v1
 * </E> \endcode where \a v0 and \a v1 are the ID's of the two
 * vertices that define the edge (Note that these vertices are
 * passed in the same order to the constructor of the
 * edge). Now, the orientation of an edge of a two-dimensional
 * element (i.e. quadrilateral or triangle) is defined as: -
 * <b> Forward </b> if the vertex with ID \a v0 comes before
 * the vertex with ID \a v1 when considering the vertices the
 * vertices of the element in a counterclockwise direction -
 * <b> Backward </b> otherwise.<BR>
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
 * We will follow
 * the basic principles of the connectivity strategy as
 * explained in Section 4.2.1.1 of [1] (such as the hierarchic
 * ordering of the edge modes).  However, we do not follow the
 * strategy described to negate the odd modes of an
 * intersecting edge of two adjacent elements if the local
 * coordinate systems have an opposite direction. The
 * explained strategy involves checking the direction of the
 * local coordinate systems of the neighbouring elements.
 * However, for a simpler automatic procedure to identify
 * which edges need to have odd mode negated, we would like to
 * have an approach which can be applied to the elements
 * individually, without information being coupled between
 * neighbouring elements.<BR> This can be accomplished in the
 * following way:<BR> Note that this approach is based on the
 * earlier observation that the intersecting edge of two
 * elements always has opposite orientation. Proper
 * connectivity can now be guaranteed if: - forward oriented
 * edges always have a counterclockwise local coordinate axis
 * - backward oriented edges always have a clockwise local
 * coordinate axis.
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
 * This algorithm above is based on the earlier observation
 * that the local coordinate axis automatically point in
 * counterclockwise direction for edges 0 and 1 and in
 * clockwise direction for the other edges.  As explained in
 * [1] the change in local coordinate axis can actually be
 * done by reversing the sign of the odd modes. This is
 * implemented by means of an additional sign vector.
 * 
 * \paragraph Nodalexpansions Nodal expansions
 * For the nodal
 * expansions, we will use the connectivity strategy as
 * explained in Section 4.2.1.1 of [1].  However, we will
 * clarify this strategy from a Nektar++ point of view. As
 * pointed out in [1], the nodal edge modes can be identified
 * with the physical location of the nodal points. In order to
 * ensure proper connectivity between elements the egde modes
 * with the same nodal location should be matched. This will
 * be accomplished if both the sets of local edge modes along
 * the intersection edge of two elements are numbered in the
 * same direction. And as the intersecting edge of two
 * elements always has opposite direction, this can be
 * guaranteed if: - the local numbering of the edge modes is
 * counterclockwise for forward oriented edges - the local
 * numbering of the edge modes is clockwise for backward
 * oriented edges.
 *
 * This will ensure that the numbering of the global DOF's on
 * an edge is in the same direction as the tow subsets of
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
