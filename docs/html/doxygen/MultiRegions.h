namespace Nektar {
namespace MultiRegions {
/**
 * @namespace Nektar::MultiRegions
 * @brief The namespace associated with the the MultiRegions library.
 * (@ref pageMultiRegions "MultiRegions introduction").
 *
 * This namespace encompasses all those classes related to multiple element
 * formulations, from uncoupled lists of elemental expansions, to coupled
 * fields employing continuous and discontinuous formulations with boundary
 * conditions. For more details see
 * @ref pageMultiRegions "MultiRegions introduction".
 */

/**
 * @page pageMultiRegions MultiRegions sublibrary
 * In the MultiRegions library, all classes and routines are related to
 * the process of assembling a global spectral/hp expansion out of local
 * elemental contributions are bundled together. The most important entities
 * of this library are the base class ExpList and its daughter classes.
 * These classes all are the abstraction of a multi-elemental spectral/hp
 * element expansion. Three different types of multi-elemental expansions
 * can be distinguished:<BR><BR>
 * <b> A collection of local expansions.</b><BR>
 * This collection is just a list of local expansions, without any coupling
 * between the expansions on the different elements, and can be formulated
 * as:
 * \f[ u^{\delta}(\boldsymbol{x})=\sum_{e=1}^{{N_{\mathrm{el}}}}
 * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}) \f]
 * where
 * - \f${N_{\mathrm{el}}}\f$ is the number of elements,
 * - \f$N^{e}_m\f$ is the number of local expansion modes within the
 *   element \f$e\f$,
 * - \f$\phi_n^e(\boldsymbol{x})\f$ is the \f$n^{th}\f$ local expansion
 *   mode within the element \f$e\f$,
 * - \f$\hat{u}_n^e\f$ is the \f$n^{th}\f$ local expansion coefficient
 * - within the element \f$e\f$.
 *
 * These types of expansion are represented by the classes
 * ExpList1D, ExpList2D and ExpList3D, depending on the dimension of the
 * problem.
 *
 * <b>A multi-elemental discontinuous global expansion.</b><BR>
 * The expansions are represented by the classes
 * DisContField1D, DisContField2D and DisContField3D. Objects of these classes
 * should be used when solving partial differential equations using a
 * discontinuous Galerkin approach. These classes enforce a coupling between
 * elements and augment the domain with boundary conditions.
 *
 * All local elemental expansions are now connected to form a global
 * spectral/hp representation. This type of global expansion can be defined as:
 * \f[u^{\delta}(\boldsymbol{x})=\sum_{n=0}^{N_{\mathrm{dof}}-1}\hat{u}_n
 * \Phi_n(\boldsymbol{x})=\sum_{e=1}^{{N_{\mathrm{el}}}}
 * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x})\f]
 * where
 * - \f$N_{\mathrm{dof}}\f$ refers to the number of global modes,
 * - \f$\Phi_n(\boldsymbol{x})\f$ is the \f$n^{th}\f$ global expansion mode,
 * - \f$\hat{u}_n\f$ is the \f$n^{th}\f$ global expansion coefficient.
 *
 * Typically, a mapping array to relate the global degrees of freedom
 * \f$\hat{u}_n\f$ and local degrees of freedom \f$\hat{u}_n^e\f$ is
 * required to assemble the global expansion out of the local
 * contributions. <BR>
 * - In order to solve (second-order) partial differential equations,
 *   information about the boundary conditions should be incorporated in
 *   the expansion. In case of a standard Galerkin implementation, the
 *   Dirichlet boundary conditions can be enforced by lifting a known
 *   solution satisfying these conditions, leaving a homogeneous Dirichlet
 *   problem to be solved. If we denote the unknown solution by
 *   \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet boundary
 *   conditions by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$ then we can
 *   decompose the solution \f$u^{\delta}(\boldsymbol{x})\f$ into the form
 *   \f[u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
 *   u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
 *   \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
 *   \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}-1}
 *   \hat{u}_n^{\mathcal{H}}\Phi_n(\boldsymbol{x}_i).\f]
 *   Implementation-wise, the known solution can be lifted by ordering the
 *   known degrees of freedom \f$\hat{u}_n^{\mathcal{H}}\f$ first in the
 *   global solution array \f$\boldsymbol{\hat{u}}\f$.<BR>
 *
 * <b> A multi-elemental continuous global expansion.</b><BR>
 * The discontinuous case is supplimented with a global continuity condition.
 * In this case a \f$C^0\f$ continuity condition is imposed across the element
 * interfaces and the expansion is therefore globally continuous.
 * - This type of global continuous expansion which incorporates the
 *   boundary conditions are represented by the classes ContField1D,
 *   ContField2D and ContField3D. Objects of these classes should be used
 *   when solving partial differential equations using a standard Galerkin
 *   approach.
 *
 * <b> Additional classes.</b><BR>
 * Furthermore, we have two more sets of classes:
 * - The class LocalToGlobalBaseMap and its daughter classes:
 *   LocalToGlobalC0ContMap and LocalToGlobalDGMap.<BR>
 *   These classes are an abstraction of the mapping from local to global
 *   degrees of freedom and contain one or both of the following mapping
 *   arrays:
 *   - map[\f$e\f$][\f$n\f$]<BR>
 *     This array contains the index of the global degree of freedom
 *     corresponding to the \f$n^{th}\f$ local expansion mode within the
 *     \f$e^{th}\f$ element.
 *   - bmap[\f$e\f$][\f$n\f$]<BR>
 *     This array contains the index of the global degree of freedom
 *     corresponding to the \f$n^{th}\f$ local boundary mode within the
 *     \f$e^{th}\f$ element.
 *   Next to the mapping array, these classes also contain routines to
 *   assemble the global system from the local contributions, and other
 *   routines to transform between local and global level.
 * - The classes GlobalLinSys and GlobalLinSysKey.<BR>
 *   The class GlobalLinSys is an abstraction of the global system matrix
 *   resulting from the global assembly procedure. Depending of the choice
 *   to statically condense the global matrix or not, the relevant blocks
 *   are stored as a member of this class. Given a proper right hand side
 *   vector, this class also contains a routine to solve the resulting matrix
 *   system.<BR>
 *   The class GlobalLinSysKey represents a key which uniquely defines a
 *   global matrix. This key can be used to construct or retrieve the global
 *   matrix associated to a certain key.
 *
 * More information about the implementation of connectivity between elements
 * in Nektar++ can be found \ref pageConnectivity "here".

 * \image html MultiRegions.png
**/
}
}
