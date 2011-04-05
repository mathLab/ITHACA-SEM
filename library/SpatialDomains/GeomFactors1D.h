#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS1D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS1D_H

#include <SpatialDomains/GeomFactors.h>
#include <StdRegions/StdExpansion2D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        // Forward declarations
        class GeomFactors1D;

        /// Shared pointer to a GeomFactors1D object.
        typedef boost::shared_ptr<GeomFactors1D>      GeomFactors1DSharedPtr;
        /// Vector of shared pointers to GeomFactors1D objects.
        typedef std::vector< GeomFactors1DSharedPtr > GeomFactors1DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors1DVector::iterator         GeomFactors1DVectorIter;

        /// Geometric factors for a 1D expansion.
        class GeomFactors1D : public GeomFactors
        {
        public:
            /// One dimensional geometric factors based on one-, two- or three-
            /// dimensional coordinate description.
            SPATIAL_DOMAINS_EXPORT GeomFactors1D(const GeomType gtype,
                          const int coordim,
                          const Array<OneD, const StdRegions
                                            ::StdExpansion1DSharedPtr> &Coords,
                          const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                          const bool QuadMetrics = false,
                          const bool LaplMetrics  = false);

            /// Copy constructor
            SPATIAL_DOMAINS_EXPORT GeomFactors1D(const GeomFactors1D& S);

            /// Destructor
            SPATIAL_DOMAINS_EXPORT virtual ~GeomFactors1D();

        private:
            /// Set up 1D Jacobian
            void SetUpJacGmat1D();

            /// Set up edge normals (given 2D element)
            virtual void v_ComputeNormals(
                        const GeometrySharedPtr &geom,
                        const int edge,
                        const LibUtilities::PointsKey &to_key);

            /**
             * @todo Implement tangents.
             * @todo Implement 1D quadrature metrics.
             * @todo Implement 1D laplacian metrics.
             */
        };
    }
}

#endif
