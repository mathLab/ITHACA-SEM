///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors2D.h
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
// Description: 2D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS2D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS2D_H

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion2D.h>  // for StdExpansion2DSharedPtr

namespace Nektar
{
    namespace SpatialDomains
    {
        class GeomFactors2D : public GeomFactors
        {
        public:
            /// One dimensional geometric factors based on one-, two- or three-
            /// dimensional coordinate description.
            SPATIAL_DOMAINS_EXPORT GeomFactors2D(const GeomType gtype,
                                                 const int coordim,
                                                 const Array<OneD, const StdRegions
                                                 ::StdExpansion2DSharedPtr> &Coords,
                                                 const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                 &tbasis,
                                                 const bool QuadMetrics = false,
                                                 const bool LaplMetrics  = false,
                                                 const bool CheckJacPositive = true);

            /// Copy constructor.
            SPATIAL_DOMAINS_EXPORT GeomFactors2D(const GeomFactors2D& S);

            /// Destructor.
            SPATIAL_DOMAINS_EXPORT virtual ~GeomFactors2D();

        private:
            /// Set up 1D Jacobian.
            void SetUpJacGmat2D(bool CheckJacPositive = true);

            /// Compute vectors in a principle direction.
            void ComputePrincipleDirection(
                        Array<OneD,Array<OneD,NekDouble> > &output);
            
            /// Set up surface normals
            virtual void v_ComputeSurfaceNormals();
            
            /// Set up the tangent vectors
            virtual void v_ComputeTangents();

            /// Set up quadrature metrics for 2D expansions.
            virtual void v_SetUpQuadratureMetrics(LibUtilities::ShapeType shape,
                const Array<OneD, const LibUtilities::BasisSharedPtr>
                    &tbasis);

            /// Set up Laplacian metrics for 2D expansions.
            virtual void v_SetUpLaplacianMetrics(
                LibUtilities::ShapeType shape,
                const Array<OneD, const LibUtilities::BasisSharedPtr>
                    &tbasis);
        };

        /// Shared pointer to GeomFactors2D object.
        typedef boost::shared_ptr<GeomFactors2D>      GeomFactors2DSharedPtr;
        /// Vector of shared pointers to GeomFactors2D objects.
        typedef std::vector< GeomFactors2DSharedPtr > GeomFactors2DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors2DVector::iterator GeomFactors2DVectorIter;
    }
}

#endif
