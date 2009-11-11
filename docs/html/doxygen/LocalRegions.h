namespace Nektar {
namespace LocalRegions {
/**
 * @namespace Nektar::LocalRegions
 * @brief The namespace associated with the the LocalRegions library
 * (@ref pageLocalRegions "LocalRegions introduction")
 *
 * This namespace contains classes related to local expansions on elements in
 * the physical space. It specialises the standard element formulations in
 * reference space, as provided in the Nektar#StdRegions namespace with a
 * spatial locality and orientation. For further details see the
 * @ref pageLocalRegions "LocalRegions introduction".
 *
 *
 * @page pageLocalRegions The %LocalRegions library
 * This provides extensions of the spectral element formulation into the world.
 * It provides spatially local forms of the reference space expansions through
 * a one-to-one linear mapping from a standard straight-sided region to the
 * physical space, based on the vertices.
 *
 * @section secLocalRegionsLocalMappings Local Mappings
 * <b>Linear Mappings</b><BR>
 * In one dimension this has the form
 * \f[ x = \chi(\xi) = \frac{1-\xi}{2}x_{e-1} + \frac{1+\xi}{2}x_e 
 *  \quad \xi \Omega^e.\f]
 * In two dimensions, for a quadrilateral, each coordinate is given by
 * \f[ x_i = \chi(\xi_1,\xi_2) = x_i^A\frac{1-\xi_1}{2}\frac{1-\xi_2}{2}
 *                              + x_i^B\frac{1+\xi_1}{2}\frac{1-\xi_2}{2}
 *                              + x_i^D\frac{1-\xi_1}{2}\frac{1+\xi_2}{2}
 *                              + x_i^C\frac{1+\xi_1}{2}\frac{1+\xi_2}{2},
 *                              \quad i=1,2.
 * \f]
 * <b>Curvilinear mappings</b><BR>
 * The mapping can be extended to curved-sided regions through the use of an
 * iso-parametric representation. In contrast to the linear mapping, where only
 * information about the vertices of the element were required, a curvilinear
 * mapping requires information about the shape of each side. This is provided
 * by shape functions, \f$ f^A(\xi_1), f^B(\xi_2), f^C(\xi_1) \f$ and \f$
 * f^D(\xi_2) \f$, in the local coordinate system. For example, the linear
 * blending function is given by
 * \f[x_i = \chi_i(\xi_1,\xi_2) = f^A(\xi_1)\frac{1-\xi_2}{2}
 *                              + f^C(\xi_1)\frac{1+\xi_2}{2}
 *                              + f^B(\xi_2)\frac{1-\xi_1}{2}
 *                              + f^D(\xi_2)\frac{1+\xi_1}{2}\\
 *                              - \frac{1-\xi_1}{2}\frac{1-\xi_2}{2}f^A(-1)
 *                              - \frac{1+\xi_1}{2}\frac{1-\xi_2}{2}f^A(1)
 *                              - \frac{1-\xi_1}{2}\frac{1+\xi_2}{2}f^C(-1)
 *                              - \frac{1+\xi_1}{2}\frac{1+\xi_2}{2}f^C(1)
 * \f]
 *
 * @section secLocalRegionsClasses Classes
 * All local expansions are derived from the top level Expansion base class.
 * Three classes, Expansion1D, Expansion2D and Expansion3D, are derived from
 * this and provided base classes for expansions in one-, two- and threei-
 * dimensions, respectively. The various local expansions are derived from
 * these.
 *
 * <b>One dimension</b><BR>
 * - SegExp - Line expansion, local version of StdRegions#StdSegExp.
 *
 * <b>Two dimensions</b><BR>
 * - TriExp - Triangular expansion.
 * - QuadExp - Quadrilateral expansion.
 *
 * <b>Three dimensions</b><BR>
 * - TetExp - Tetrehedral expansion. (All triangular faces)
 * - HexExp - Hexahedral expansion. (All rectangular faces)
 * - PrismExp - Prism expansion. (Two triangular, three rectangular faces)
 * - PyrExp - Pyramid expansion. (One rectangular, four triangular faces)
 *
 * <b>Other classes</b>
 * - PointExp
 * - LinSys
 * - MatrixKey
 */
}
}
