#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        // Forward declarations
        class GeomFactors3D;

        /// Shared pointer to GeomFactors3D object.
        typedef boost::shared_ptr<GeomFactors3D>      GeomFactors3DSharedPtr;
        /// Vector of shared pointers to GeomFactors3D objects.
        typedef std::vector< GeomFactors3DSharedPtr > GeomFactors3DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors3DVector::iterator GeomFactors3DVectorIter;

        /// Geometric factors for a 3D expansions.
        class GeomFactors3D : public GeomFactors
        {
        public:
            /// One dimensional geometric factors based on one-, two- or three-
            /// dimensional coordinate description.
            SPATIAL_DOMAINS_EXPORT GeomFactors3D(const GeomType gtype,
                          const int coordim,
                          const Array<OneD, const StdRegions
                                            ::StdExpansion3DSharedPtr> &Coords,
                          const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                          const bool QuadMetrics = false,
                          const bool LaplMetrics  = false);

        private:
            // Set up 3D Jacobian
            void SetUpJacGmat3D(
                            const Array<OneD, Array<OneD, NekDouble> > d1,
                            const Array<OneD, Array<OneD, NekDouble> > d2,
                            const Array<OneD, Array<OneD, NekDouble> > d3);

            /// Set up quadrature metrics for 2D expansions.
            virtual void v_SetUpQuadratureMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);
        };
    }
}

#endif
