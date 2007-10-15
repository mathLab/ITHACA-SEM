///////////////////////////////////////////////////////////////////////////////
//
// File StdBasis.cpp
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
// Description: Basis definition 
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/LibUtilities.h>

#include <math.h>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>


namespace Nektar
{
    namespace LibUtilities 
    {

        bool operator<(const BasisKey &lhs, const BasisKey &rhs)
        {
            PointsKey lhsPointsKey = lhs.GetPointsKey();
            PointsKey rhsPointsKey = rhs.GetPointsKey();

            if (lhsPointsKey  < rhsPointsKey) 
            {
                return true;
            }
            if (lhsPointsKey != rhsPointsKey) 
            {
                return false;
            }

            if (lhs.m_nummodes < rhs.m_nummodes)
            {
                return true;
            }
            if (lhs.m_nummodes > rhs.m_nummodes) 
            {
                return false;
            }

            return (lhs.m_basistype < rhs.m_basistype);
        }

        bool operator>(const BasisKey &lhs, const BasisKey &rhs)
        {
            return (rhs < lhs);
        }

        bool BasisKey::opLess::operator()(const BasisKey &lhs, const BasisKey &rhs) const
        {
            return (lhs.m_basistype < rhs.m_basistype);
        }

        std::ostream& operator<<(std::ostream& os, const BasisKey& rhs)
        {
            os << "NumModes: " << rhs.GetNumModes() << " BasisType: " << BasisTypeMap[rhs.GetBasisType()];
            os << " " << rhs.GetPointsKey() << std::endl;

            return os;
        }

        boost::shared_ptr<Basis> Basis::Create(const BasisKey &bkey)
        {
            boost::shared_ptr<Basis> returnval(new Basis(bkey));

            returnval->Initialize();

            return returnval;
        }

        void Basis::Initialize()
        {
            ASSERTL0(GetNumModes()>0, "Cannot call Basis initialisation with zero or negative order");
            ASSERTL0(GetTotNumPoints()>0, "Cannot call Basis initialisation with zero or negative numbers of points");

            // Allocate Memory
            int size = GetTotNumModes()*GetTotNumPoints();
            m_bdata  = Array<OneD, NekDouble>(size);
            m_dbdata = Array<OneD, NekDouble>(size);

            GenBasis();
        };

        // Method used to generate appropriate basis
        void Basis::GenBasis()
        {
            int i,p,q;
            NekDouble scal;
            Array<OneD, NekDouble> modeSharedArray;
            NekDouble *mode;
            ConstArray<OneD, NekDouble> z;
            ConstArray<OneD, NekDouble> w;
            const NekDouble *D;

            boost::shared_ptr< Points<NekDouble> > pointsptr = PointsManager()[GetPointsKey()];
            pointsptr->GetZW(z,w);

            D = &(pointsptr->GetD()->GetPtr())[0];

            int numModes = GetNumModes();
            int numPoints = GetNumPoints();

            switch(GetBasisType())
            {
            case eOrtho_A:
            case eLegendre:
                mode = m_bdata.data();

                for (p=0; p<numModes; ++p, mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z.data(), mode, NULL, p, 0.0, 0.0);

                    // normalise 
                    scal = sqrt(0.5*(2.0*p+1.0));
                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= scal;
                    }
                }

                // define derivative basis
                Blas::Dgemm('n','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints,m_bdata.data(),numPoints,0.0,
                    m_dbdata.data(),numPoints);
                break;


            case eOrtho_B:
                {
                    const NekDouble *one_m_z_pow;//,fac;

                    // bdata should be of size order*(order+1)/2*zorder;

                    mode = m_bdata.data();


                    for(i = 0; i < numPoints; ++i)
                        mode[i] = 1.0;

                    mode += numPoints;

                    for (q = 1; q < numModes; ++q, mode += numPoints)
                    {
                        Polylib::jacobfd(numPoints, z.data(), mode, NULL, q, 1.0, 0.0);
                    }

                    one_m_z_pow = m_bdata.data();

                    for(p = 1; p < numModes; ++p)
                    {

                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] = 0.5*(1-z[i])*one_m_z_pow[i];
                        }

                        one_m_z_pow = mode;
                        mode += numPoints;

                        for(q = 1; q < numModes-p; ++q, mode+=numPoints){
                            Polylib::jacobfd(numPoints, z.data(), mode, NULL, q, 2.0*p+1.0, 0.0);

                            for(i = 0; i < numPoints; ++i)
                            {
                                mode[i] *= one_m_z_pow[i];
                            }
                        }
                    }

                    // normalise (recalling factor of 2 for weight (1-b)/2) 
                    for(p = 0, mode=m_bdata.data(); p < numModes; ++p)
                    {
                        for(q = 0; q < numModes-p; ++q,mode+=numPoints)
                        {
                            scal = sqrt((p+q+1.0));
                            for(i = 0; i < numPoints; ++i)
                            {
                                mode[i] *= scal;
                            }
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('n','n',numPoints,numModes*(numModes+1)/2,numPoints,1.0,D,numPoints,
                        m_bdata.data(),numPoints,0.0,m_dbdata.data(),numPoints);
                }
                break; 

            // This is tilde psi_pqr in Spen's book, page 105
            // The 4-dimensional array is laid out in memory such that 
            // 1) Eta_z values are the changing the fastest, then r, q, and finally p.
            // 2) r index increases by the stride of numPoints.
            case eOrtho_C:  
                {
                    int P = numModes - 1, Q = numModes - 1, R = numModes - 1;
                    NekDouble *mode = m_bdata.data();

                    for( int p = 0, m = 0; p <= P; ++p ) {
                        for( int q = 0; q <= Q - p; ++q ) {
                            for( int r = 0; r <= R - p - q; ++r, mode += numPoints ) {
                                Polylib::jacobfd(numPoints, z.data(), mode, NULL, r, 2*p + 2*q + 1.0, 0.0);
                                for( int k = 0; k < numPoints; ++k ) {
                                    mode[k] *= 0.5*pow(1.0 - z[k], p+q);
                                }
                            }
                        }
                    }


                    // Define derivative basis
                    Blas::Dgemm('n','n',numPoints,numModes*(numModes+1)*
                        (numModes+2)/6,numPoints,1.0, D, numPoints,
                        m_bdata.data(),numPoints,0.0,m_dbdata.data(),numPoints);
                }
                break;

            case eModified_A:

                for(i = 0; i < numPoints; ++i)
                {
                    m_bdata[i] = 0.5*(1-z[i]);
                    m_bdata[numPoints + i] = 0.5*(1+z[i]);
                }

                mode = m_bdata.data() + 2*numPoints;

                for(p = 2; p < numModes; ++p, mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z.data(), mode, NULL, p-2,1.0,1.0);

                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= m_bdata[i]*m_bdata[numPoints+i];
                    }
                }

                // define derivative basis 
                Blas::Dgemm('n','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints,m_bdata.data(),numPoints,0.0,m_dbdata.data(),
                    numPoints);
                break;

            case eModified_B: case eModified_C:
                {

                    const NekDouble *one_m_z_pow, *one_p_z;

                    // bdata should be of size order*(order+1)/2*zorder

                    // first fow 
                    for(i = 0; i < numPoints; ++i)
                    {
                        m_bdata[0*numPoints + i] = 0.5*(1-z[i]);
                        m_bdata[1*numPoints + i] = 0.5*(1+z[i]);
                    }

                    mode = m_bdata.data() + 2*numPoints;

                    for(q = 2; q < numModes; ++q, mode+=numPoints)
                    {
                        Polylib::jacobfd(numPoints, z.data(), mode, NULL, q-2,1.0,1.0);

                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] *= m_bdata[i]*m_bdata[numPoints+i];
                        }
                    }

                    // second row
                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] = 0.5*(1-z[i]);
                    }

                    mode += numPoints;

                    for(q = 2; q < numModes; ++q, mode+=numPoints)
                    {
                        Polylib::jacobfd(numPoints, z.data(), mode, NULL, q-2,1.0,1.0);

                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] *= m_bdata[i]*m_bdata[numPoints+i];
                        }
                    }

                    // third and higher rows 
                    one_m_z_pow = m_bdata.data();
                    one_p_z     = m_bdata.data()+numPoints;

                    for(p = 2; p < numModes; ++p)
                    {
                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] = m_bdata[i]*one_m_z_pow[i];
                        }

                        one_m_z_pow  = mode;
                        mode        += numPoints;

                        for(q = 1; q < numModes-p; ++q, mode+=numPoints)
                        {
                            Polylib::jacobfd(numPoints,z.data(),mode,NULL,q-1,2*p+1,1.0);

                            for(i = 0; i <  numPoints; ++i)
                            {
                                mode[i] *= one_m_z_pow[i]*one_p_z[i];
                            }
                        }
                    }

                    Blas::Dgemm('n','n',numPoints,numModes*(numModes+1)/2,
                        numPoints,1.0,D,numPoints,
                        m_bdata.data(),numPoints,0.0,m_dbdata.data(),numPoints);
                }
                break;

            case eGLL_Lagrange: 
                {
                    mode = m_bdata.data();
                    boost::shared_ptr< Points<NekDouble> > pointsptr = PointsManager()[PointsKey(numModes,eGaussLobattoLegendre)];
                    const ConstArray<OneD, NekDouble>& zp(pointsptr->GetZ());

                    for (p=0; p<numModes; ++p,mode += numPoints)
                    {
                        for(q = 0; q < numPoints; ++q)
                        {
                            mode[q] = Polylib::hglj(p,z[q],zp.data(),numModes,0.0,0.0);
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('n','n',numPoints,numModes,numPoints,1.0,D,
                        numPoints, m_bdata.data(),numPoints,0.0,
                        m_dbdata.data(),numPoints);

                }//end scope
                break;
            case eFourier:

                ASSERTL0(numModes%2==0, "Fourier modes should be a factor of 2");

                for(i = 0; i < numPoints; ++i)
                {
                    m_bdata[i] = 1.0;
                    m_bdata[numPoints+i] = 0.0; 

                    m_dbdata[i] = m_dbdata[numPoints+i] = 0.0; 
                }

                for (p=1; p < numModes/2; ++p)
                {
                    for(i = 0; i < numPoints; ++i)
                    {
                        m_bdata[ 2*p   *numPoints+i] = cos(p*M_PI*z[i]);
                        m_bdata[(2*p+1)*numPoints+i] = sin(p*M_PI*z[i]);

                        m_dbdata[ 2*p   *numPoints+i] = -p*M_PI*sin(p*M_PI*z[i]);
                        m_dbdata[(2*p+1)*numPoints+i] =  p*M_PI*cos(p*M_PI*z[i]);
                    }
                }

                break;

            case eChebyshev:

                mode = m_bdata.data();

                for (p=0,scal = 1; p<numModes; ++p,mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z.data(), mode, NULL, p, -0.5, -0.5);

                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= scal;
                    }

                    scal *= 4*(p+1)*(p+1)/(2*p+2)/(2*p+1);    
                }

                // define derivative basis 
                Blas::Dgemm('n','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints, m_bdata.data(),numPoints,0.0,m_dbdata.data(),
                    numPoints);
                break;

            default:
                ASSERTL0(false, "Basis Type not known or "
                    "not implemented at this time.");
            }
        }

        /** \brief Determine if polynomial basis can be eactly integrated
        *  with itself
        */
        bool BasisKey::ExactIprodInt(void) const 
        {
            bool returnval = false;

            switch(GetPointsType())
            {
            case eGaussGaussLegendre:
            case eGaussRadauMLegendre:
            case eGaussRadauPLegendre:
            case eGaussLobattoLegendre:
                returnval = (GetNumPoints() >= GetNumModes());
                break;
            }

            return returnval;
        }


        /** \brief Determine if basis has collocation property,
        *  i.e. GLL_Lagrange with Lobatto integration of appropriate order.
        */
        bool BasisKey::Collocation() const 
        {
            return ( m_basistype == eGLL_Lagrange &&
                GetPointsType() == eGaussLobattoLegendre &&
                GetNumModes() == GetNumPoints());
        }


        // BasisKey compared to BasisKey
        bool operator  == (const BasisKey& x, const BasisKey& y)
        {
            return (x.GetPointsKey() == y.GetPointsKey() &&
                x.m_basistype == y.m_basistype &&
                x.GetNumModes() == y.GetNumModes());
        }


        // BasisKey* compared to BasisKey
        bool operator  == (const BasisKey* x, const BasisKey& y)
        {
            return (*x == y);
        }

        // \brief BasisKey compared to BasisKey*
        bool operator  == (const BasisKey& x, const BasisKey *y)
        {
            return (x == *y);
        }

        // \brief BasisKey compared to BasisKey
        bool operator  != (const BasisKey& x, const BasisKey& y)
        {
            return (!(x == y));
        }

        // BasisKey* compared to BasisKey
        bool operator  != (const BasisKey* x, const BasisKey& y)
        {
            return (!(*x == y));
        }


        // BasisKey compared to BasisKey*
        bool operator  != (const BasisKey& x, const BasisKey* y)
        {
            return (!(x == *y));
        }

    } // end of namespace stdregion
} // end of namespace stdregion

/** 
* $Log: Basis.cpp,v $
* Revision 1.23  2007/10/03 03:00:13  bnelson
* Added precompiled headers.
*
* Revision 1.22  2007/09/27 12:52:03  pvos
* Column major Blas calls corrections
*
* Revision 1.21  2007/09/25 14:24:40  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.20  2007/07/22 23:03:26  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.19  2007/07/20 00:28:24  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.18  2007/05/15 03:37:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.17  2007/04/29 00:31:56  jfrazier
* Updated to use multi_arrays.
*
* Revision 1.16  2007/04/08 03:31:18  jfrazier
* Modified to use SharedArray.
*
* Revision 1.15  2007/04/04 02:10:30  bnelson
* Made comparison operators const correct.
*
* Revision 1.14  2007/04/03 03:58:24  bnelson
* Moved Lapack.hpp, Blas.hpp, Transf77.hpp to LinearAlgebra
*
* Revision 1.13  2007/02/26 15:52:30  sherwin
* Working version for Fourier points calling from StdRegions. Fourier interpolations not working yet
*
* Revision 1.12  2007/02/17 04:06:47  jfrazier
* Added greater-than operator for easy comparison when using keys that use a basis and provided a method to access the basiskey.
*
* Revision 1.11  2007/02/06 17:12:27  jfrazier
* Fixed a problem with global initialization in libraries.
*
* Revision 1.10  2007/02/01 23:28:41  jfrazier
* Basis is not working, but not fully tested.
*
* Revision 1.9  2007/02/01 15:30:07  jfrazier
* Reformatted.
*
Revision 1.8  2007/01/31 23:27:56  kirby
*** empty log message ***

* Revision 1.7  2007/01/31 18:37:32  kirby
*
* fully compiling but not fully tested update to Foundations
* Basis now should work (Compiles but not tested yet)
*
* Revision 1.6  2007/01/29 19:47:29  jfrazier
* Cleanup and restructuring to accommodate private default and copy constructors.
*
* Revision 1.5  2007/01/26 04:00:54  jfrazier
* Cleanup after initial creation.
*
* Revision 1.4  2007/01/24 23:43:01  kirby
* *** empty log message ***
*
* Revision 1.3  2007/01/20 22:21:17  kirby
* *** empty log message ***
*
* Revision 1.2  2007/01/17 11:35:52  pvos
* updating doxygen documentation
*
* Revision 1.1  2006/06/06 18:45:14  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:30  kirby
* *** empty log message ***
*
* Revision 1.30  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.29  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.28  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.27  2006/03/01 17:07:33  sherwin
*
* Added new location of polylib in header definitions
*
* Revision 1.26  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.25  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.24  2006/02/19 13:26:13  sherwin
*
* Coding standard revisions so that libraries compile
*
**/ 

