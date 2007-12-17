////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/GeomFactors.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H

#include <SpatialDomains/SpatialDomains.hpp>

#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpansion2D.h>
#ifdef HIGH_D_FUNCTIONS
#include <StdRegions/StdExpansion3D.h>
#endif

namespace Nektar
{
    namespace SpatialDomains
    {
        class GeomFactors;
        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs);

        typedef boost::shared_ptr<GeomFactors>      GeomFactorsSharedPtr;
        typedef std::vector< GeomFactorsSharedPtr > GeomFactorsVector;
        typedef GeomFactorsVector::iterator GeomFactorsVectorIter;

        class GeomFactors
        {
        public:

            friend bool operator==(const GeomFactors &lhs, const GeomFactors &rhs);

            GeomFactors(void);

            GeomFactors(const GeomType gtype, const int expdim, const int coordim);

            /** \brief One dimensional geometric factors based on one,
            two or three dimensional coordinate description
            **/
            GeomFactors(const GeomType gtype, const int coordim,
                const ConstArray<OneD, StdRegions::StdExpansion1DSharedPtr> &Coords);

            /**  \brief Two dimensional geometric factors based on two
            or three dimensional coordinate description
            **/
            GeomFactors(const GeomType gtype, const int coordim,
                        const ConstArray<OneD,StdRegions::StdExpansion2DSharedPtr> &Coords);

#ifdef HIGH_D_FUNCTIONS
            /**  \brief Three dimensional geometric factors and Jacobian
            **/
            GeomFactors(const GeomType gtype, 
                const StdRegions::StdExpansion3D **Coords);
#endif

            ~GeomFactors();

            inline GeomType GetGtype()
            {
                return m_gtype;
            }

            inline const ConstArray<TwoD,NekDouble> &GetGmat() const
            {
                return m_gmat;
            }

            inline const ConstArray<OneD,NekDouble> &GetJac() const 
            {
                return m_jac;
            }

            inline const ConstArray<TwoD,NekDouble> &GetNormals() const
            {
                return m_normals;
            }

            
            inline void  ResetNormals(const ConstArray<TwoD,NekDouble> &newnorm)
            {
                m_normals = Array<TwoD,NekDouble>(newnorm.GetRows(),newnorm.GetColumns(),newnorm.data());
            }

            inline void ResetGmat(const ConstArray<OneD,NekDouble> &ndata, 
                                  const int nq, const int expdim, 
                                  const int coordim)
            {
                m_gmat = Array<TwoD,NekDouble>(expdim*coordim,nq,ndata.data());
            }


            inline void ResetGmat(const ConstArray<TwoD,NekDouble> &ndata)
            {
                m_gmat = Array<TwoD,NekDouble>(ndata.GetRows(),ndata.GetColumns(),ndata.data());
            }

            inline void ResetJac(int nq, const ConstArray<OneD,NekDouble> &ndata)
            {
                m_jac = Array<OneD, NekDouble>(nq, ndata.data());
            }

        protected:
            Array<TwoD,NekDouble> m_gmat;
            Array<OneD,NekDouble> m_jac;
            Array<TwoD,NekDouble> m_normals;

            GeomType m_gtype;
            int m_expdim;
            int m_coordim;
        };
    } //end of namespace
} //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GeomFactors_H

//
// $Log: GeomFactors.h,v $
// Revision 1.14  2007/12/03 21:30:43  sherwin
// Added normal details
//
// Revision 1.13  2007/07/22 23:04:23  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.12  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.11  2007/07/10 22:20:59  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.10  2007/07/10 17:06:31  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.9  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.8  2007/05/28 08:35:26  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.7  2007/05/25 17:52:02  jfrazier
// Updated to use new Array classes.
//
// Revision 1.6  2007/05/17 18:45:25  jfrazier
// Minor changes to accommodate Array class.
//
// Revision 1.5  2007/04/08 03:34:48  jfrazier
// Updated to compile with SharedArray.  This has not been converted to SharedArray, just made to work with others that have been converted.
//
// Revision 1.4  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.3  2007/03/29 19:24:00  bnelson
// Replaced boost::shared_array with SharedArray
//
// Revision 1.2  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.1  2007/03/20 09:17:39  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.7  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.6  2007/03/02 12:01:59  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.5  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.4  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.3  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.2  2006/05/29 17:05:17  sherwin
// Updated to use shared_ptr around Geom types - added typedef
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
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

