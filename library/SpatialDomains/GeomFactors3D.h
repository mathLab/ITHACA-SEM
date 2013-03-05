///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors3D.cpp
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
// Description: Implementation of 3D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion3D.h>

namespace Nektar
{

    namespace SpatialDomains
    {
        // Forward declarations

        /// Geometric factors for a 3D expansions.
        class GeomFactors3D : public GeomFactors
        {
        public:
            /// One dimensional geometric factors based on one-, two- or three-
            /// dimensional coordinate description.
            SPATIAL_DOMAINS_EXPORT GeomFactors3D(
                const GeomType gtype,
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
                LibUtilities::ShapeType shape,
                const Array<OneD, const LibUtilities::BasisSharedPtr>
                    &tbasis);
        };

        /// Shared pointer to GeomFactors3D object.
        typedef boost::shared_ptr<GeomFactors3D>      GeomFactors3DSharedPtr;
        /// Vector of shared pointers to GeomFactors3D objects.
        typedef std::vector< GeomFactors3DSharedPtr > GeomFactors3DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors3DVector::iterator GeomFactors3DVectorIter;

    }
}

#endif
