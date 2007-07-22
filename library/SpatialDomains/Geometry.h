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

#include "pchSpatialDomains.h"

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/GeomFactors.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry; // Forward declaration for typedef.
        typedef boost::shared_ptr<Geometry> GeometrySharedPtr;
        typedef std::vector< GeometrySharedPtr > GeometryVector;
        typedef std::vector< GeometrySharedPtr >::iterator GeometryVectorIter;

        class Geometry
        {
            public:
                Geometry();
                Geometry(int coordim);

                virtual ~Geometry();

                inline GeomType GetGtype()
                {
                    return m_geomfactors->GetGtype();
                }

                inline const ConstArray<OneD,NekDouble> &GetJac()
                {
                    return m_geomfactors->GetJac();
                }

                inline const ConstArray<TwoD, NekDouble>& GetGmat()
                {
                    return m_geomfactors->GetGmat();
                }

                inline const int GetCoordim() const 
                {
                    return m_coordim;
                }

                inline GeomFactorsSharedPtr GetGeomFactors(void)
                {
                    if(!(m_geomfactors.get()))
                    {
                        GenGeomFactors();
                    }

                    return m_geomfactors;
                }


                // Wrappers around virtual Functions
                void GenGeomFactors(void)
                {
                    return v_GenGeomFactors();
                }

            protected:

                static GeomFactorsSharedPtr ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor);

                int                  m_coordim;     // coordinate dimension
                GeomFactorsSharedPtr m_geomfactors;
                GeomState            m_state;       // enum identifier to determine if quad points are filled
                static GeomFactorsVector m_RegGeomFactorsManager;

            private:
                virtual void v_GenGeomFactors(void)
                {
                    NEKERROR(ErrorUtil::efatal,
                        "This function is only valid for shape type geometries");
                }

        };
    }; //end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY_H

//
// $Log: Geometry.h,v $
// Revision 1.12  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.11  2007/07/10 22:21:00  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.10  2007/07/10 17:06:31  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.9  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.8  2007/05/25 17:52:02  jfrazier
// Updated to use new Array classes.
//
// Revision 1.7  2007/03/20 09:17:39  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.6  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.5  2007/01/18 20:59:28  sherwin
// Before new configuration
//
// Revision 1.4  2006/08/18 19:37:17  jfrazier
// *** empty log message ***
//
// Revision 1.3  2006/08/17 22:55:00  jfrazier
// Continued adding code to process composites in the mesh2d.
//
// Revision 1.2  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
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



