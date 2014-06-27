/** 

* \page pageStructure Structure
 
A major challenge which arises when one aims to develop a software
package that implements the spectral/hp element method is to implement
the mathematical structure of the method in a digestible and coherent
matter. Obviously, there are many ways to encapsulate the fundamental
concepts related to the spectral/hp element method, depending on e.g.
the intended goal of the developer or the chosen programming
language. We will (without going in too much
detail) give a an overview of how we have chosen to abstract and
implement spectral/hp elements in the \a Nektar++ library. 
However, we want to emphasise that this is not the only possible choice.

Five different sublibraries, employing this characteristic pattern,
are provided in the full \a Nektar++ library:
 
 - the standard elemental region sublibrary (@link Nektar#StdRegions StdRegions @endlink library)
 - the parametric mapping sublibrary (@link Nektar#SpatialDomains SpatialDomains @endlink library)
 - the local elemental region sublibrary (@link Nektar#LocalRegions LocalRegions @endlink library)
 - the global region sublibrary (@link Nektar#MultiRegions MultiRegions @endlink library)
 - the supporting utilities sublibrary (@link Nektar#LibUtilities LibUtilities @endlink library)
 
This structure can also be related to the formulation of a global
spectral/hp element expansion, i.e.

\f[
 
u(\boldsymbol{x})=\overbrace{\sum_{e\in\mathcal{E}}\underbrace{\sum_{n\in\mathcal{N}}\phi^e_n(\boldsymbol{x})\hat{u}^e_n}_{\mbox{\scriptsize{LocalRegions library}}}}^{\mbox{\scriptsize{MultiRegions library}}}=\sum_{e\in\mathcal{E}}\underbrace{\sum_{n\in\mathcal{N}}\phi^{std}_n\overbrace{(\left[\chi^e\right]^{-1}(\boldsymbol{x}))}^{\mbox{\scriptsize{SpatialDomains library}}}\hat{u}^e_n}_{\mbox{\scriptsize{StdRegions library}}}.
\f]
 
A more detailed overview of the \a Nektar++ structure, including an
overview of the most important classes per sublibrary, is depicted in
the figure below.
 
* \image html overview.png "Main structure of \a Nektar++"
 
*/

