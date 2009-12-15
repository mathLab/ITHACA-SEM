#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS2D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS2D_H

#include <SpatialDomains/GeomFactors.h>
#include <StdRegions/StdExpansion3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        // Forward declaration
        class GeomFactors2D;

        /// Shared pointer to GeomFactors2D object.
        typedef boost::shared_ptr<GeomFactors2D>      GeomFactors2DSharedPtr;
        /// Vector of shared pointers to GeomFactors2D objects.
        typedef std::vector< GeomFactors2DSharedPtr > GeomFactors2DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors2DVector::iterator GeomFactors2DVectorIter;

        /// Geometric factors for a 2D expansions.
        class GeomFactors2D : public GeomFactors
        {
        public:
            /// One dimensional geometric factors based on one-, two- or three-
            /// dimensional coordinate description.
            GeomFactors2D(const GeomType gtype,
                          const int coordim,
                          const Array<OneD, const StdRegions
                                            ::StdExpansion2DSharedPtr> &Coords,
                          const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                          const bool QuadMetrics = false,
                          const bool LaplMetrics  = false);

            /// Copy constructor.
            GeomFactors2D(const GeomFactors2D& S);

            /// Destructor.
            virtual ~GeomFactors2D();

        private:
            /// Set up 1D Jacobian.
            void SetUpJacGmat2D();

            /// Compute vectors in a principle direction.
            void ComputePrincipleDirection(
                        Array<OneD,Array<OneD,NekDouble> > &output);

            /// Compute outward edge normals along an edge.
            virtual void v_ComputeEdgeNormals(
                        const int edge,
                        const LibUtilities::PointsKey &to_key,
                        Array<OneD, Array<OneD, NekDouble> > &output) const;

            /// Set up surface normals
            virtual void v_ComputeSurfaceNormals();

            /// Set up the tangent vectors
            virtual void v_ComputeTangents();

            /// Set up quadrature metrics for 2D expansions.
            virtual void v_SetUpQuadratureMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

            /// Set up Laplacian metrics for 2D expansions.
            virtual void v_SetUpLaplacianMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

        };
    }
}

#endif
