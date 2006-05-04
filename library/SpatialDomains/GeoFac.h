////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/GeoFac.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOFAC_H
#define NEKTAR_SPATIALDOMAINS_GEOFAC_H

#include <SpatialDomains/SpatialDomains.hpp>

#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdExpansion3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class GeoFac
        {
        public:
            GeoFac(const StdRegions::GeomType gtype, const int expdim, 
                const int coordim);

            /** \brief One dimensional geometric factors based on one,
            two or three dimensional coordinate description
            **/
            GeoFac(const StdRegions::GeomType gtype, const int coordim,
                const StdRegions::StdExpansion1D **Coords);

            /**  \brief Two dimensional geometric factors based on two
            or three dimensional coordinate description
            **/
            GeoFac(const StdRegions::GeomType gtype, const int coordim,
                const StdRegions::StdExpansion2D **Coords);

            /**  \brief Three dimensional geometric factors and Jacobian
            **/
            GeoFac(const StdRegions::GeomType gtype, 
                const StdRegions::StdExpansion3D **Coords);

            ~GeoFac();

            inline int GetMemCount()
            {
                return m_mem_count;
            }

            inline void AddMemCount()
            {
                m_mem_count += 1;
            }

            inline void SubMemCount()
            {
                m_mem_count -= 1;
            }

            inline StdRegions::GeomType GetGtype()
            {
                return m_gtype;
            }

            inline const double ** GetGmat() const
            {
                return (const double**) m_gmat;
            }

            inline const double *GetJac() const 
            {
                return m_jac;
            }

            inline void ResetGmat(double *ndata, int nq, int expdim, 
                int coordim)
            {
                if(m_gmat[0])
                {
                    delete[] m_gmat[0];
                }

                m_gmat[0] = ndata;
                for(int i = 1; i < expdim*coordim; ++i)
                {
                    m_gmat[i] = m_gmat[i-1]+nq;
                }
            }

            inline void ResetJac(double *ndata)
            {
                if(m_jac)
                { 
                    delete[] m_jac;
                }

                m_jac = ndata;
            }

        protected:
            double ** m_gmat;
            double *  m_jac;

        private:
            StdRegions::GeomType m_gtype;
            int m_mem_count;
        };

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOFAC_H

//
// $Log: GeoFac.h,v $
// Revision 1.16  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.15  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.14  2006/03/13 18:20:03  sherwin
//
// Fixed error in ResetGmat:
//
// Revision 1.13  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.12  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.11  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

