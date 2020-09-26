///////////////////////////////////////////////////////////////////////////////
//
// File Expansion0D.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header file for Expansion0D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION0D_H
#define EXPANSION0D_H

#include <LocalRegions/Expansion.h>
#include <StdRegions/StdExpansion0D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/Expansion1D.h>
#include <SpatialDomains/Geometry0D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        class Expansion0D;
        typedef std::shared_ptr<Expansion0D>      Expansion0DSharedPtr;
        typedef std::vector< Expansion0DSharedPtr > Expansion0DVector;

        class Expansion0D: virtual public Expansion,
            virtual public StdRegions::StdExpansion0D
        {
        public:
            LOCAL_REGIONS_EXPORT Expansion0D(
                SpatialDomains::Geometry0DSharedPtr pGeom);
            
            LOCAL_REGIONS_EXPORT virtual ~Expansion0D() {}
                        
            inline SpatialDomains::Geometry0DSharedPtr GetGeom0D() const;

        protected:
			
        private:
        };
        

        inline SpatialDomains::Geometry0DSharedPtr
            Expansion0D::GetGeom0D() const
        {
            return std::dynamic_pointer_cast<SpatialDomains::
                Geometry0D>(m_geom);
			
        }
    } //end of namespace
} //end of namespace

#endif
