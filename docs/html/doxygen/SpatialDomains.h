namespace Nektar
{
namespace SpatialDomains
{
/** \page pageSpatialDomains SpatialDomains library
 
 The most important family of classes in the SpatialDomains library is
 the Geometry family, as can also be seen in the figure below. 
 These classes are the representation of a
 (geometric) element in \a physical \a space}. It has been indicated
 before that every local element can be
 considered as an image of the standard element where the corresponding
 one-to-one mapping can be represented as an elemental standard
 spectral/hp expansion. As such, a proper encapsulation should at least
 contain data structures that represent such an expansion in order to
 completely define the geometry of the element. Therefore, we have
 equipped the classes in the Geometry family with the
 following data structures:

* - an object of @link Nektar#StdRegions#StdExpansion StdExpansion @endlink class, and
* - a data structure that contains the metric terms (Jacobian, derivative metrics, \f\ldots\f) of the transformation.

 Note that although the latter data structure is not necessary to
 define the geometry, it contains information inherent to the
 iso-parametric representation of the element that can later be used in
 e.g. the LocalRegions library.  Again, the @link Nektar#StdRegions#StdExpansion StdExpansion @endlink object
 can be defined in the abstract base class Geometry. However,
 for every shape-specific geometry class, it needs to be initialised
 according to the corresponding StdRegions class (e.g. for the
 QuadGeom class, it needs to be initialised as an
 @link Nektar#StdRegions#StdQuadExp StdQuadExp @endlink object).
 
 * \image html SpatialDomains.png

*/

}
}
