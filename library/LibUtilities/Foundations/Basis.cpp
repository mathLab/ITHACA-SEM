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


#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Blas.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        bool operator<(const BasisKey &lhs, const BasisKey &rhs)
        {
            PointsKey lhsPointsKey = lhs.GetPointsKey();
            PointsKey rhsPointsKey = rhs.GetPointsKey();

            if (lhsPointsKey < rhsPointsKey) return true;
            if (lhsPointsKey != rhsPointsKey) return false;

            if (lhs.m_nummodes < rhs.m_nummodes) return true;
            if (lhs.m_nummodes > rhs.m_nummodes) return false;

            return (lhs.m_basistype < rhs.m_basistype);
        }

        bool BasisKey::opLess::operator()(const BasisKey &lhs, const BasisKey &rhs)
        {
#pragma message("Address this operator")
            return (lhs<rhs);
        }

        std::ostream& operator<<(std::ostream& os, const BasisKey& rhs)
        {
            return os;
        }

        
        void Basis::Initialize()
        {
            ASSERTL0(GetNumModes()>0, "Cannot call Basis initialisation with zero or negative order");
            ASSERTL0(GetTotNumPoints()>0, "Cannot call Basis initialisation with zero or negative numbers of points");

            // Allocate Memory
            int size = GetTotNumModes()*GetTotNumPoints();
            m_bdata  = new double [size];
            m_dbdata = new double [size];

            GenBasis();
        };

        // Method used to generate appropriate basis
        void Basis::GenBasis()
        {
            int i,p,q;
            double scal,*mode;
            const double *z, *w, *D;

            std::cout << "I AM HERE" << std::endl;


            boost::shared_ptr<Points<double> > pointsptr = PointsManager()[GetPointsKey()];
            pointsptr->GetZW(z,w);
            D = &(pointsptr->GetD()->GetPtr())[0];

            int numModes = GetNumModes();
            int numPoints = GetNumPoints();

            switch(GetBasisType())
            {
            case eOrtho_A:
            case eLegendre:
                mode = m_bdata;

                for (p=0; p<numModes; ++p, mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z, mode, NULL, p, 0.0, 0.0);

                    // normalise 
                    scal = sqrt(0.5*(2.0*p+1.0));
                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= scal;
                    }
                }

                // define derivative basis
                Blas::Dgemm('t','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints,m_bdata,numPoints,0.0,
                    m_dbdata,numPoints);
                break;


            case eOrtho_B:
                {
                    double *one_m_z_pow;//,fac;

                    // bdata should be of size order*(order+1)/2*zorder;

                    mode = m_bdata;

                    for(i = 0; i < numPoints; ++i)
                        mode[i] = 1.0;

                    mode += numPoints;

                    for (q = 1; q < numModes; ++q, mode += numPoints)
                    {
                        Polylib::jacobfd(numPoints, z, mode, NULL, q, 1.0, 0.0);
                    }

                    one_m_z_pow = m_bdata;

                    for(p = 1; p < numModes; ++p)
                    {

                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] = 0.5*(1-z[i])*one_m_z_pow[i];
                        }

                        one_m_z_pow = mode;
                        mode += numPoints;

                        for(q = 1; q < numModes-p; ++q, mode+=numPoints){
                            Polylib::jacobfd(numPoints, z, mode, NULL, q, 2.0*p+1.0, 0.0);

                            for(i = 0; i < numPoints; ++i)
                            {
                                mode[i] *= one_m_z_pow[i];
                            }
                        }
                    }

                    // normalise (recalling factor of 2 for weight (1-b)/2) 
                    for(p = 0, mode=m_bdata; p < numModes; ++p)
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
                    Blas::Dgemm('t','n',numPoints,numModes*(numModes+1)/2,numPoints,1.0,D,numPoints,
                        m_bdata,numPoints,0.0,m_dbdata,numPoints);
                }
                break; 

            case eOrtho_C:
                {

                    double *one_m_z_pow,scal;//,fac;

                    // bdata should be of size _order*(order+1)*(order+2)/6*zorder;

                    mode = m_bdata;
                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[0] = 1.0;
                    }

                    mode += numPoints;
                    for (q = 1; q < numModes; ++q, mode += numPoints)
                    {
                        Polylib::jacobfd(numPoints, z, mode, NULL, q, 2.0, 0.0);
                    }

                    one_m_z_pow = m_bdata;
                    for(p = 1; p < numModes; ++p)
                    {
                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] = 0.5*(1-z[i])*one_m_z_pow[i];
                        }

                        one_m_z_pow = mode;
                        mode += numPoints;

                        for(q = 1; q < numModes-p; ++q,mode+=numPoints)
                        {
                            Polylib::jacobfd(numPoints, z, mode, NULL, q, 2.0*p+2.0, 0.0);

                            for(i = 0; i < numPoints; ++i)
                            {
                                mode[i] *= one_m_z_pow[i];
                            }
                        }
                    }

                    ASSERTL2(false, "Normalisation might need fixing");

                    for(p = 0, mode=m_bdata; p < numModes; ++p)
                    {
                        for(q = 0; q < numModes-p; ++q,mode+=numPoints)
                        {
                            scal = sqrt((p+q+1.5));
                            for(i = 0; i < numPoints; ++i)
                            {
                                mode[i] *= scal;
                            }
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('t','n',numPoints,numModes*(numModes+1)*
                        (numModes+2)/6,numPoints,1.0,D,numPoints,
                        m_bdata,numPoints,0.0,m_dbdata,numPoints);
                }       
                break;

            case eModified_A:

                for(i = 0; i < numPoints; ++i)
                {
                    m_bdata[i] = 0.5*(1-z[i]);
                    m_bdata[numPoints + i] = 0.5*(1+z[i]);
                }

                mode = m_bdata + 2*numPoints;

                for(p = 2; p < numModes; ++p, mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z, mode, NULL, p-2,1.0,1.0);

                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= m_bdata[i]*m_bdata[numPoints+i];
                    }
                }

                // define derivative basis 
                Blas::Dgemm('t','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints,m_bdata,numPoints,0.0,m_dbdata,
                    numPoints);
                break;
            case eModified_B: case eModified_C:
                {

                    double *one_m_z_pow,*one_p_z;

                    // bdata should be of size order*(order+1)/2*zorder

                    // first fow 
                    for(i = 0; i < numPoints; ++i)
                    {
                        m_bdata[0*numPoints + i] = 0.5*(1-z[i]);
                        m_bdata[1*numPoints + i] = 0.5*(1+z[i]);
                    }

                    mode = m_bdata + 2*numPoints;

                    for(q = 2; q < numModes; ++q, mode+=numPoints)
                    {
                        Polylib::jacobfd(numPoints, z, mode, NULL, q-2,1.0,1.0);

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
                        Polylib::jacobfd(numPoints, z, mode, NULL, q-2,1.0,1.0);

                        for(i = 0; i < numPoints; ++i)
                        {
                            mode[i] *= m_bdata[i]*m_bdata[numPoints+i];
                        }
                    }

                    // third and higher rows 
                    one_m_z_pow = m_bdata;
                    one_p_z     = m_bdata+numPoints;

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
                            Polylib::jacobfd(numPoints,z,mode,NULL, q-1,2*p+1,1.0);

                            for(i = 0; i <  numPoints; ++i)
                            {
                                mode[i] *= one_m_z_pow[i]*one_p_z[i];
                            }
                        }
                    }

                    Blas::Dgemm('t','n',numPoints,numModes*(numModes+1)/2,
                        numPoints,1.0,D,numPoints,
                        m_bdata,numPoints,0.0,m_dbdata,numPoints);
                }
                break;

            case eGLL_Lagrange: 
                {
                    const double *zp,*wp;

                    mode = m_bdata;
                    (PointsManager()[PointsKey(numModes,eGaussLobattoLegendre)])->GetZW(zp,wp);

                    for (p=0; p<numModes; ++p,mode += numPoints)
                    {
                        for(q = 0; q < numPoints; ++q)
                        {
                            mode[q] = Polylib::hglj(p,z[q],zp,numModes,0.0,0.0);
                        }
                    }

                    // define derivative basis 
                    Blas::Dgemm('t','n',numPoints,numModes,numPoints,1.0,D,
                        numPoints, m_bdata,numPoints,0.0,
                        m_dbdata,numPoints);

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

                mode = m_bdata;

                for (p=0,scal = 1; p<numModes; ++p,mode += numPoints)
                {
                    Polylib::jacobfd(numPoints, z, mode, NULL, p, -0.5, -0.5);

                    for(i = 0; i < numPoints; ++i)
                    {
                        mode[i] *= scal;
                    }

                    scal *= 4*(p+1)*(p+1)/(2*p+2)/(2*p+1);	
                }

                // define derivative basis 
                Blas::Dgemm('t','n',numPoints,numModes,numPoints,1.0,D,
                    numPoints, m_bdata,numPoints,0.0,m_dbdata,
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
            return (m_basistype == eGLL_Lagrange &&
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




