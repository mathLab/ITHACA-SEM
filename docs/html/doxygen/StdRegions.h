namespace Nektar
{
namespace StdRegions
{
/** \page pageStdRegions StdRegions library

 The StdRegions library, see the figure below, bundles
 all classes that mimic a spectral/hp element expansion on a standard
 region. 
 Such an expansion, i.e.
 \f[
 u(\boldsymbol{\xi}_i) = \sum_{n\in\mathcal{N}}\phi_n(\boldsymbol{\xi}_i)\hat{u}_n,
 
 \f]
 can be encapsulated in a class that essentially should only contain
 three data structures, respectively representing:
 
* - the coefficient vector \f$\hat{\boldsymbol{u}}\f$,
* - the discrete basis matrix \f$\boldsymbol{B}\f$, and
* - the vector \f$\boldsymbol{u}\f$ which represents the value of the
* - expansion at the quadrature points \f$\boldsymbol{\xi}_i\f$.
 
 All standard expansions, independent of the dimensionality or shape of
 the standard region, can be abstracted in a similar way. Therefore, it
 is possible to define these data structures in an \a abstract base
 class, i.e. the class StdExpansion. This base class can also
 contain the implementation of methods that are identical across all
 shapes. Derived from this base class is another level of abstraction,
 i.e. the abstract classes StdExpansion1D,
 StdExpansion2D and StdExpansion3D.  All other
 shape-specific classes (such as e.g. StdSegExp or
 StdQuadExp) are inherited from these abstract base
 classes. These shape-specific classes are the classes from which
 objects will be instantiated. They also contain the shape-specific
 implementation for operations such as integration or differentiation.
 
 * \image html StdRegions.png
 
 */

}
}
