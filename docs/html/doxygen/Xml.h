namespace Nektar
{
/**
 * \page pageXML Nektar++ XML File Format
 * The Nektar++ native file format is compliant with XML version 1.0. The root
 * element is \c NEKTAR and has the overall structure as follows
 * @code
 * <NEKTAR>
 *   <GEOMETRY>
 *   ...
 *   </GEOMETRY>
 *   <EXPANSIONS>
 *   ...
 *   </EXPANSIONS>
 *   <CONDITIONS>
 *   ...
 *   </CONDITIONS>
 * </NEKTAR>
 * @endcode
 *
 * \section sectionXMLGeometry XML Geometry definition
 * This section defines the mesh. It specifes a list of vertices, edges (in two
 * or three dimensions) and faces (in three dimensions) and how they connect to
 * create the elemental decomposition of the domain. It also defines a list of
 * composites which are used in the Expansions and Conditions sections of the
 * file to describe the polynomial expansions and impose boundary conditions.
 *
 * The \c GEOMETRY section is structured as
 * @code
 * <GEOMETRY DIM="2" SPACE="2">
 *   <VERTEX>
 *   ...
 *   </VERTEX>
 *   <EDGE>
 *   ...
 *   </EDGE>
 *   <FACE>
 *   ...
 *   </FACE>
 *   <ELEMENT>
 *   ...
 *   </ELEMENT>
 *   <CURVED>
 *   ...
 *   </CURVED>
 *   <COMPOSITE>
 *   ...
 *   </COMPOSITE>
 *   <DOMAIN> ... </DOMAIN>
 * </GEOMETRY>
 * @endcode
 * It has two attributes:
 * - \c DIM: specifies the dimension of the expansion elements.
 * - \c SPACE: specifies the dimension of the space in which the elements exist.
 *
 * These attributes allow, for example, a two-dimensional surface to be
 * embedded in a three-dimensional space. The \c FACES section is only present
 * when \c DIM=3. The \c CURVED section is only present if curved edges or
 * faces are present in the mesh.
 *
 * \subsection subsectionXmlGeometryVertices Vertices
 * Vertices have three coordinates. Each has a unique vertex ID. They are
 * defined in the file with a line of the form
 * @code
 * <V ID="0"> 0.0  0.0  0.0 </V>
 * @endcode
 *
 * \subsection subsectionXmlGeometryEdges Edges
 * Edges are defined by two vertices. Each edge has a unique edge ID. They are
 * defined in the file with a line of the form
 * @code
 * <E ID="0"> 0 1 </E>
 * @endcode
 *
 * \subsection subsectionXmlGeometryFaces Faces
 * Faces are defined by three or more edges. Each face has a unique face ID.
 * They are defined in the file with a line of the form
 * @code
 * <T ID="0"> 0 1 2 </T>
 * <Q ID="1"> 3 4 5 6 </Q>
 * @endcode
 * The choice of tag specified (\c T or \c Q), and thus the number of edges
 * specified depends on the geometry of the face (triangle or quadrilateral).
 *
 * \subsection subsectionXmlGeometryElements Element
 * Elements define the top-level geometric entities in the mesh. Their
 * definition depends upon the dimension of the expansion. For two-dimensional
 * expansions, an element is defined by a sequence of three or four edges. For
 * three-dimensional expansions, the element is defined by a list of faces.
 * Elements are defined in the file with a line of the form
 * @code
 * <T ID="0"> 0 1 2 </T>
 * <H ID="1"> 3 4 5 6 7 8 </H>
 * @endcode
 * Again, the choice of tag specified depends upon the geometry of the element.
 * The element tags are:
 * - \c S: Segment
 * - \c T: Triangle
 * - \c Q: Quadrilateral
 * - \c A: Tetrahedron
 * - \c P: Pyramid
 * - \c R: Prism
 * - \c H: Hexahedron
 *
 * \subsection subsectionXmlGeometryCurved Curved Edges and Faces
 * For mesh elements with curved edges and/or curved faces, a separate entry is
 * used to describe the control points for the curve. Each curve has a unique
 * curve ID and is associated with a predefined edge or face. The total number
 * of points in the curve (including end points) and their distribution is also
 * included as attributes. The control points are listed in order, each
 * specified by three coordinates. Curved edges are defined in the file with a
 * line of the form
 * @code
 * <E ID="3" EDGEID="7" TYPE="PolyEvenlySpaced" NUMPOINTS="3">
 *         0.0  0.0  0.0    0.5  0.5  0.0    1.0  0.0  0.0
 * </E>
 * @endcode
 *
 * \subsection subsectionXmlGeometryComposites Composites
 * Composites define collections of elements, faces or edges. Each has a unique
 * composite ID associated with it. All components of a composite entry must be
 * of the same type. The syntax allows components to be listed individually or
 * using ranges. Examples include
 * @code
 * <C ID="0"> T[0-862] </C>
 * <C ID="1"> E[68,69,70,71] </C>
 * @endcode
 *
 * \subsection subsectionXmlGeometryDomain Domain
 * This tag specifies composites which describe the entire problem domain. It
 * has the form of
 * @code
 * <DOMAIN> C[0] </DOMAIN>
 * @endcode
 *
 * \section sectionXMLExpansions XML Expansions definition
 * This section defines the polynomial expansions used on each of the defined
 * geometric composites. Expansion entries specify the number of modes, the
 * basis type and have the form
 * @code
 * <E COMPOSITE="C[0]" NUMMODES="5" TYPE="MODIFIED" />
 * @endcode
 *
 * \section sectionXMLConditions XML Conditions definition
 * The final section of the file defines parameters and boundary conditions
 * which define the nature of the problem to be solved. These are enclosed in
 * the \c CONDITIONS tag.
 *
 * \subsection subsectionXmlConditionsParameters Parameters
 * Parameters may be required by a particular solver (for instance
 * time-integration parameters or solver-specific parameters), or arbitrary and
 * only used within the context of the session file (e.g. parameters in the
 * definition of an initial condition). All parameters are enclosed in the
 * \c PARAMETERS XML element.
 * @code
 *   <PARAMETERS>
 *   ...
 *   </PARAMETERS>
 * @endcode
 * A parameter may be of integer or real type and may reference other
 * parameters defined previous to it. It is expressed in the file as
 * @code
 * <P> [PARAMETER NAME] = [PARAMETER VALUE] </P>
 * @endcode
 * For example,
 * @code
 * <P> NumSteps = 1000              </P>
 * <P> TimeStep = 0.01              </P>
 * <P> FinTime  = NumSteps*TimeStep </P>
 * @endcode
 *
 * \subsection subsectionXmlConditionsSolverInfo Solver Information
 * These specify properties to define the actions specific to solvers,
 * typically including the equation to solve, the projection type and the
 * method of time integration. The property/value pairs are specified as
 * XML attributes. For example,
 * @code
 *   <SOLVERINFO>
 *     <I PROPERTY="EQTYPE"                VALUE="UnsteadyAdvection"    />
 *     <I PROPERTY="Projection"            VALUE="Continuous"           />
 *     <I PROPERTY="TimeIntegrationMethod" VALUE="ClassicalRungeKutta4" />
 *   </SOLVERINFO>
 * @endcode
 *
 * \subsection subsectionXmlConditionsVariables Variables
 * These define the number (and name) of solution variables. Each variable is
 * prescribed a unique ID. For example a two-dimensional flow simulation may
 * define the velocity variables \f$(u,v)\f$ using
 * @code
 *   <VARIABLES>
 *     <V ID="0"> u </V>
 *     <V ID="1"> v </V>
 *   </VARIABLES>
 * @endcode
 *
 * \subsection subsectionXmlConditionsBoundary Boundary Regions and Conditions
 * Boundary conditions are defined by two XML elements. The first defines the
 * various boundary regions in the domain in terms of composite entities from
 * the \c GEOMETRY section of the file. Each boundary region has a
 * unique ID and are defined as, for example,
 * @code
 *   <BOUNDARYREGIONS>
 *     <B ID="0"> C[2] </B>
 *     <B ID="1"> C[3] </B>
 *   </BOUNDARYREGIONS>
 * @endcode
 * The second defines the actual boundary condition to impose on that composite
 * region for each of the defined solution variables, and has the form,
 * @code
 *   <BOUNDARYCONDITIONS>
 *     <REGION REF="0">
 *       <D VAR="u" VALUE="sin(PI*x)*cos(PI*y)" />
 *       <D VAR="v" VALUE="sin(PI*x)*cos(PI*y)" />
 *     </REGION>
 *   </BOUNDARYCONDITIONS>
 * @endcode
 * Boundary condition specifications may refer to any parameters defined in the
 * session file. The \c REF attribute corresponds to a defined boundary region.
 * The tag used for each variable specifies the type of boundary condition to
 * enforce. These can be either
 * - \c D: Dirichlet
 * - \c N: Neumann
 * - \c R: Robin
 * - \c P: Periodic
 *
 * Time-dependent boundary conditions may be specified through setting the
 * \c USERDEFINEDTYPE attribute. For example,
 * @code
 * <D VAR="u" USERDEFINEDTYPE="TimeDependent" VALUE="sin(PI*(x-t))" />
 * @endcode
 *
 * Periodic boundary conditions reference the corresponding boundary region
 * with which to enforce periodicity.
 *
 * The following example provides an example of three boundary conditions for a
 * two-dimensional flow,
 * @code
 * <BOUNDARYCONDITIONS>
 *   <REGION REF="0">
 *     <D VAR="u" USERDEFINEDTYPE="TimeDependent" VALUE="-cos(x)*sin(y)*exp(-2*t*Kinvis)" />
 *     <D VAR="v" USERDEFINEDTYPE="TimeDependent" VALUE="sin(x)*cos(y)*exp(-2*t*Kinvis)" />
 *     <P VAR="p" VALUE=[2]/>
 *   </REGION>
 *   <REGION REF="1">
 *     <D VAR="u" USERDEFINEDTYPE="TimeDependent" VALUE="-cos(x)*sin(y)*exp(-2*t*Kinvis)" />
 *     <D VAR="v" USERDEFINEDTYPE="TimeDependent" VALUE="sin(x)*cos(y)*exp(-2*t*Kinvis)" />
 *     <N VAR="p" USERDEFINEDTYPE="H" VALUE="0.0"/>
 *   </REGION>
 *   <REGION REF="2">
 *     <D VAR="u" USERDEFINEDTYPE="TimeDependent" VALUE="-cos(x)*sin(y)*exp(-2*t*Kinvis)" />
 *     <D VAR="v" USERDEFINEDTYPE="TimeDependent" VALUE="sin(x)*cos(y)*exp(-2*t*Kinvis)" />
 *     <P VAR="p" VALUE=[0]/>
 *   </REGION>
 *   <REGION REF="3">
 *     <D VAR="u" USERDEFINEDTYPE="TimeDependent" VALUE="-cos(x)*sin(y)*exp(-2*t*Kinvis)" />
 *     <D VAR="v" USERDEFINEDTYPE="TimeDependent" VALUE="sin(x)*cos(y)*exp(-2*t*Kinvis)" />
 *     <D VAR="p" USERDEFINEDTYPE="TimeDependent" VALUE="-0.25*(cos(2*x)+cos(2*y))*exp(-4*t*Kinvis)"/>
 *   </REGION>
 * </BOUNDARYCONDITIONS>
 * @endcode
 * where the boudary regions which are periodic are linked via their region
 * identifier (Region 0 and Region 2).
 *
 *
 * \subsection subsectionXmlConditionsUserDefined User-defined Equations
 * These define spatially-dependent analytic expressions for use in solvers.
 * For example,
 * @code
 *   <USERDEFINEDEQNS>
 *     <F LHS="Vx" VALUE="1" />
 *     <F LHS="Vy" VALUE="1" />
 *   </USERDEFINEDEQNS>
 * @endcode
 *
 * \subsection subsectionXmlConditionsAnalytic Analytic expressions
 * Finally, initial conditions and analytic solutions may be specified for use
 * in, or comparison with, simulations.
 * @code
 *   <EXACTSOLUTION>
 *     <F VAR="u" VALUE="sin(PI*x-advx*t))*cos(PI*(y-advy*t))" />
 *   </EXACTSOLUTION>
 *   <INITIALCONDITIONS>
 *     <F VAR="u" VALUE="sin(PI*x)*cos(PI*y)" />
 *   </INITIALCONDITIONS>
 * @endcode
 Initial conditions can also be specified via a re-start file (typically named as FILENAME.rst).
 A re-start file is basically a solution file (in other words an .fld renamed as .rst) where the solution is specified.
 The expansion order used to generate the .rst file has to be the same as the one you intend to use for your simulation and
 the file has to be located in the same directory where the .xml file and the the executable are.
 An example is:
 * @code
 *   <INITIALCONDITIONS>
 *     <R FILE="Test_KovaFlow_m3.rst" />
 *   </INITIALCONDITIONS>
 * @endcode
 which uses a converged solution to initialise the flow, avoiding the long initial transient state.
 *
 */
}
