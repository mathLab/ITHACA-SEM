////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/Geometry1D.h,v $
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry1D: public Geometry
        {
        public:
            Geometry1D();
            Geometry1D(int coordim);
            ~Geometry1D();

        protected:

        private:
        };

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

//
// $Log: Geometry1D.h,v $
// Revision 1.2  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.19  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.18  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.17  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.16  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
