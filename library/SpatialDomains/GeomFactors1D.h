///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors1D.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// 
// Description: 1D geometric factors class.
//
///////////////////////////////////////////////////////////////////////////////

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

            /// Set up edge tangents (given 1D element)
            virtual void v_ComputeEdgeTangents(
            	    	const GeometrySharedPtr &geom,
            	    	const int edge,
            	    	const LibUtilities::PointsKey &to_key);            

            /**
             * @todo Implement 1D quadrature metrics.
             * @todo Implement 1D laplacian metrics.
             */
        };
    }
}

#endif
