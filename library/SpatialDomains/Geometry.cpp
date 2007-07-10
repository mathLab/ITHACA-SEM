////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.cpp
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
//  Description:  This file contains the base class implementation for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        Geometry::Geometry():
            m_coordim(0), 
            m_state(eNotFilled)	  
        {
        }

        Geometry::Geometry(const int coordim):
            m_coordim(coordim),
            m_state(eNotFilled)	  
        {
        }

        Geometry::~Geometry()
        {
        }
        
        RegGeomFactorsMap Geometry::m_RegGeomFactorsManager;
        GeomFactorsSharedPtr Geometry::ValidateRegGeomFactor(GeomFactorsSharedPtr geomFactor)
        {
            GeomFactorsSharedPtr returnval;
            if (geomFactor->GetGtype() == eRegular)
            {
                const GeomFactorsKey &geomFacKey = geomFactor->GetGeomFactorsKey();
                RegGeomFactorsMap::iterator iter = m_RegGeomFactorsManager.find(geomFacKey);

                if (iter != m_RegGeomFactorsManager.end())
                {
                    returnval = iter->second;
                }
                else
                {
                    m_RegGeomFactorsManager[geomFacKey] = geomFactor;
                    returnval = geomFactor;
                }
            }

            return returnval;
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: Geometry.cpp,v $
// Revision 1.4  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.3  2006/08/24 18:50:00  jfrazier
// Completed error checking on permissable composite item combinations.
//
// Revision 1.2  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.14  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/03/13 11:17:03  sherwin
//
// First compiing version of Demos in SpatialDomains and LocalRegions. However they do not currently seem to execute properly
//
// Revision 1.12  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.11  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
