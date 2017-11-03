///////////////////////////////////////////////////////////////////////////////
//
// File SpatialDomainsDeclarations.hpp
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
// Description: Class definition in SpatialDomains required in StdRegions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SPATIALDOMDEF_H
#define SPATIALDOMDEF_H

#include <memory>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry;
        class EdgeComponent;
        class SegGeom;
        class Geometry;
        class Geometry1D;
        class Geometry2D;
        class Geometry3D;
        class GeomFactors;

        static std::shared_ptr<GeomFactors> NullGeomFactorsSharedPtr;
        static std::shared_ptr<Geometry>    NullGeometrySharedPtr;
        static std::shared_ptr<Geometry1D>  NullGeometry1DSharedPtr;
        static std::shared_ptr<Geometry2D>  NullGeometry2DSharedPtr;
        static std::shared_ptr<Geometry3D>  NullGeometry3DSharedPtr;
    } // end of namespace
} // end of namespace

#endif //SPATIALDOMDEF_H

