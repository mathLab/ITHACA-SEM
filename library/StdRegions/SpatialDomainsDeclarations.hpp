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

namespace Nektar
{
    namespace SpatialDomains
    {
        
        class Geometry;
        class EdgeComponent;
        class SegGeom;
        class GeomFactors;
        
    } // end of namespace
} // end of namespace

#endif //SPATIALDOMDEF_H

/**
 * $Log: SpatialDomainsDeclarations.hpp,v $
 * Revision 1.4  2007/03/25 15:48:22  sherwin
 * UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
 *
 * Revision 1.3  2007/03/20 09:12:46  kirby
 * update of geofac and metric info; fix style issues
 *
 * Revision 1.2  2007/03/14 21:24:09  sherwin
 * Update for working version of MultiRegions up to ExpList1D
 *
 * Revision 1.1  2006/05/04 18:58:30  kirby
 * *** empty log message ***
 *
 * Revision 1.2  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.1  2006/02/26 23:37:29  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/
