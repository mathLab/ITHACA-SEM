/**
 * @mainpage The Nektar++ library
 *
 * Nektar++ is an open source software
 * library currently being developed and designed to provide a toolbox
 * of data structures and algorithms which implement the spectral/hp
 * element method. Nektar++ is the continuation and adaptation of the
 * <A
 * HREF="http://www2.imperial.ac.uk/ssherw/spectralhp/nektar">Nektar</A>
 * flow solver. As opposed to its predecessor which focused on solving
 * fluid dynamics problems, Nektar++ is implemented as a C++
 * object-oriented toolkit which allows developers to implement
 * spectral element solvers for a variety of different engineering
 * problems.
 *
 * The structure of the Nektar++ library, a collection of different
 * sublibraries, is based upon the typical structure of a global
 * spectral/hp approximation, which is characterized by:
 * - <b>The elemental decomposition of the problem</b><BR> As for all finite
 *   element methods, the computational domain is partitioned into a
 *   mesh of many small subdomains or elements. Analogously, the
 *   spectral/hp solution is expanded into a series of local expansions,
 *   each with support on a single element.  This elemental
 *   representation enables the treatment of operations on a local
 *   elemental basis rather than on global level. This not only
 *   simplifies the formulation but also allows many operations to be
 *   performed more efficiently.
 * - <b>The introduction of a standard region</b><BR>
 *   The introduction of a standard region allows the
 *   expansion basis to be defined just once, that is only on the
 *   standard region. All other elements then can be considered as the
 *   image of the standard element under a parametric
 *   mapping. Consequently, the elemental operations of integration and
 *   differentiation can all be executed on the standard element,
 *   subject to a proper treatment of the transformation from local
 *   (world space) to standard (reference space) coordinates. For
 *   curved-sided elements, the mapping from standard element to local
 *   element is generally done using an <em>iso-parametric</em>
 *   representation. In this case, the local geometry is represented
 *   with an expansion of the same form and polynomial order as the
 *   unknown variables.
 *
 * This structure, supplemented with building blocks such as block
 * matrix linear algebra routines and automatic data coordinating
 * objects, can be encapsulated in an efficient object-oriented C++
 * implementation.
 *
 * This conceptual
 * approach of the software leads to a high user-flexibility,
 * including the selection of the preferred expansion basis, its
 * polynomial order and the preferred numerical quadrature.
 *
 * The website of the Nektar++ project can be found on: http://www.nektar.info
 */
