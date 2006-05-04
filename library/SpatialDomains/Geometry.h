////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.h
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
//  Description:  This file contains the base class specification for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY_H

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/GeoFac.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry
        {
        public:
            Geometry();
            Geometry(int coordim);

            virtual ~Geometry();

            inline StdRegions::GeomType GetGtype()
            {
                return m_xgeofac->GetGtype();
            }

            inline const double* GetJac()
            {
                return m_xgeofac->GetJac();
            }

            inline const double** GetGmat()
            {
                return m_xgeofac->GetGmat();
            }

            inline const int GetCoordim()
            {
                return m_coordim;
            }

            inline void SetXGeoFac(GeoFac *gfac)
            {
                m_xgeofac = gfac;
            }

            inline GeoFac* GetXGeoFac()
            {
                return m_xgeofac;
            }


        protected:
            int m_coordim;     // coordinate dimension
            GeoFac* m_xgeofac;
            GeomState m_state; // enum identifier to determine if quad points are filled

        private:
            virtual void v_GenXGeoFac(void)
            {
	      ErrorUtil::Error(ErrorUtil::efatal,__FILE__, __LINE__,
                    "This function is only valid for shape type geometries");
            }

        };
    }; //end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY_H

//
// $Log: Geometry.h,v $
// Revision 1.21  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.20  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.19  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.18  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.17  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.16  2006/02/26 21:19:42  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.15  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//



