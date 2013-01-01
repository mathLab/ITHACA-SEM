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
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <boost/math/special_functions/gamma.hpp>

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

        Basis::Basis(const BasisKey &bkey):
            m_basisKey(bkey),
            m_points(PointsManager()[bkey.GetPointsKey()]),
            m_bdata(bkey.GetTotNumModes()*bkey.GetTotNumPoints()),
            m_dbdata(bkey.GetTotNumModes()*bkey.GetTotNumPoints())
        {
            m_InterpManager.RegisterGlobalCreator(boost::bind(&Basis::CalculateInterpMatrix,this,_1));
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

            GenBasis();
        };

        /** \brief Calculate the interpolation Matrix for coefficient from
        *  one base (m_basisKey) to another (tbasis0)
        */
        boost::shared_ptr< NekMatrix<NekDouble> > Basis::CalculateInterpMatrix(const BasisKey &tbasis0)
        {
            int dim = m_basisKey.GetNumModes();
            const PointsKey pkey(dim,LibUtilities::eGaussLobattoLegendre);
            BasisKey fbkey(m_basisKey.GetBasisType(),dim,pkey);
            BasisKey tbkey(tbasis0.GetBasisType(),dim,pkey);

            // "Constructur" of the basis
            BasisSharedPtr fbasis = BasisManager()[fbkey];
            BasisSharedPtr tbasis = BasisManager()[tbkey];

            // Get B Matrices
            Array<OneD, NekDouble> fB_data = fbasis->GetBdata();
            Array<OneD, NekDouble> tB_data = tbasis->GetBdata();

            // Convert to a NekMatrix
            NekMatrix<NekDouble> fB(dim,dim,fB_data);
            NekMatrix<NekDouble> tB(dim,dim,tB_data);

            // Invert the "to" matrix: tu = tB^(-1)*fB fu = ftB fu
            tB.Invert();

            // Compute transformation matrix
            Array<OneD, NekDouble> zero1D(dim*dim,0.0);
            boost::shared_ptr< NekMatrix<NekDouble> > ftB(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(dim,dim,zero1D));
            (*ftB) = tB*fB;

            return ftB;
        }

        // Method used to generate appropriate basis
             /** The following expansions are generated depending on the
             * enum type defined in \a m_basisKey.m_basistype:
             *
             * NOTE: This definition does not follow the order in the
             * Karniadakis \& Sherwin book since this leads to a more
             * compact hierarchical pattern for implementation
             * purposes. The order of these modes dictates the
             * ordering of the expansion coefficients.
             *
             * In the following m_numModes = P
             *
             * \a eModified_A:
             *
             * m_bdata[i + j*m_numpoints] =
             * \f$ \phi^a_i(z_j) = \left \{
             * \begin{array}{ll} \left ( \frac{1-z_j}{2}\right ) & i = 0 \\
             * \\
             * \left ( \frac{1+z_j}{2}\right ) & i = 1 \\
             * \\
             * \left ( \frac{1-z_j}{2}\right )\left ( \frac{1+z_j}{2}\right )
             *  P^{1,1}_{i-2}(z_j) & 2\leq i < P\\
             *  \end{array} \right . \f$
             *
             * \a eModified_B:
             *
             * m_bdata[n(i,j) + k*m_numpoints] =
             * \f$ \phi^b_{ij}(z_k) = \left \{ \begin{array}{lll}
             * \phi^a_j(z_k) & i = 0, &   0\leq j < P  \\
             * \\
             * \left ( \frac{1-z_k}{2}\right )^{i}  & 1 \leq i < P,&   j = 0 \\
             * \\
             * \left ( \frac{1-z_k}{2}\right )^{i} \left ( \frac{1+z_k}{2}\right )
             * P^{2i-1,1}_{j-1}(z_k) & 1 \leq i < P,\ &  1\leq j < P-i\ \\
             * \end{array}  \right . , \f$
             *
             * where \f$ n(i,j) \f$ is a consecutive ordering of the
             * triangular indices \f$ 0 \leq i, i+j < P \f$ where \a j
             * runs fastest.
             *
             *
             * \a eModified_C:
             *
             * m_bdata[n(i,j,k) + l*m_numpoints] =
             * \f$ \phi^c_{ij,k}(z_l) = \phi^b_{i+j,k}(z_l) =
             *  \left \{ \begin{array}{llll}
             * \phi^b_{j,k}(z_l) & i = 0, &   0\leq j < P  &  0\leq k < P-j\\
             * \\
             * \left ( \frac{1-z_l}{2}\right )^{i+j}  & 1\leq i < P,\
             * &  0\leq j <  P-i,\  & k = 0 \\
             * \\
             * \left ( \frac{1-z_l}{2}\right )^{i+j}
             * \left ( \frac{1+z_l}{2}\right )
             * P^{2i+2j-1,1}_{k-1}(z_k) & 1\leq i < P,&  0\leq j < P-i&
             * 1\leq k < P-i-j \\
             * \\
             * \end{array}  \right . , \f$
             *
             * where \f$ n(i,j,k) \f$ is a consecutive ordering of the
             * triangular indices \f$ 0 \leq i, i+j, i+j+k < P \f$ where \a k
             * runs fastest, then \a j and finally \a i.
             *
             */
        void Basis::GenBasis()
        {
            int i,p,q;
            NekDouble scal;
            Array<OneD, NekDouble> modeSharedArray;
            NekDouble *mode;
            Array<OneD, const NekDouble> z;
            Array<OneD, const NekDouble> w;
            const NekDouble *D;

            m_points->GetZW(z,w);

            D = &(m_points->GetD()->GetPtr())[0];
            int numModes = GetNumModes();
            int numPoints = GetNumPoints();

            switch(GetBasisType())
            {

            /** \brief Orthogonal basis A

            \f$\tilde \psi_p^a (\eta_1) = L_p(\eta_1) = P_p^{0,0}(\eta_1)\f$

           */
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
                Blas::Dgemm('n','n',numPoints,numModes,numPoints,1.0,D,numPoints,
                            m_bdata.data(),numPoints,0.0,m_dbdata.data(),numPoints);
                break;

            /** \brief Orthogonal basis B

            \f$\tilde \psi_{pq}^b(\eta_2) = \left ( {1 - \eta_2} \over 2 \right)^p P_q^{2p+1,0}(\eta_2)\f$ \\

           */

            // This is tilde psi_pq in Spen's book, page 105
            // The 3-dimensional array is laid out in memory such that
            // 1) Eta_y values are the changing the fastest, then q and p.
            // 2) q index increases by the stride of numPoints.
            case eOrtho_B:
                {
                     NekDouble *mode = m_bdata.data();

                     for( int p = 0; p < numModes; ++p )
                     {
                         for( int q = 0; q < numModes - p; ++q,  mode += numPoints  )
                         {
                             Polylib::jacobfd(numPoints, z.data(), mode, NULL, q, 2*p + 1.0, 0.0);
                             for( int j = 0; j < numPoints; ++j )
                             {
                                 mode[j] *= sqrt(p+q+1.0)*pow(0.5*(1.0 - z[j]), p);
                             }
                         }
                     }

                     // define derivative basis
                     Blas::Dgemm('n','n',numPoints,numModes*(numModes+1)/2,numPoints,1.0,D,numPoints,
                                  m_bdata.data(),numPoints,0.0,m_dbdata.data(),numPoints);
                }
                break;

            /** \brief Orthogonal basis C

            \f$\tilde \psi_{pqr}^c = \left ( {1 - \eta_3} \over 2 \right)^{p+q} P_r^{2p+2q+2, 0}(\eta_3)\f$ \\

           */

            // This is tilde psi_pqr in Spen's book, page 105
            // The 4-dimensional array is laid out in memory such that
            // 1) Eta_z values are the changing the fastest, then r, q, and finally p.
            // 2) r index increases by the stride of numPoints.
            case eOrtho_C:
                {
                    int P = numModes - 1, Q = numModes - 1, R = numModes - 1;
                    NekDouble *mode = m_bdata.data();

                    for( int p = 0; p <= P; ++p )
                    {
                        for( int q = 0; q <= Q - p; ++q )
                        {
                            for( int r = 0; r <= R - p - q; ++r, mode += numPoints )
                            {
                                Polylib::jacobfd(numPoints, z.data(), mode, NULL, r, 2*p + 2*q + 2.0, 0.0);
                                for( int k = 0; k < numPoints; ++k )
                                {
                                    // Note factor of 0.5 is part of normalisation
                                    mode[k] *= pow(0.5*(1.0 - z[k]), p+q);

                                    // finish normalisation
                                    mode[k] *= sqrt(r+p+q+1.5);
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

                // Note the following packing deviates from the
                // definition in the Book by Karniadakis in that we
                // put the vertex degrees of freedom at the lower
                // index range to follow a more hierarchic structure.

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

            case eModified_B:
                {

                    // Note the following packing deviates from the
                    // definition in the Book by Karniadakis in two
                    // ways. 1) We put the vertex degrees of freedom
                    // at the lower index range to follow a more
                    // hierarchic structure. 2) We do not duplicate
                    // the singular vertex definition so that only a
                    // triangular number (i.e. (modes)*(modes+1)/2) of
                    // modes are required consistent with the
                    // orthogonal basis.

                    // In the current structure the q index runs
                    // faster than the p index so that the matrix has
                    // a more compact structure

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
                            Polylib::jacobfd(numPoints,z.data(),mode,NULL,q-1,2*p-1,1.0);

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


            case eModified_C:
                {
                    // Note the following packing deviates from the
                    // definition in the Book by Karniadakis in two
                    // ways. 1) We put the vertex degrees of freedom
                    // at the lower index range to follow a more
                    // hierarchic structure. 2) We do not duplicate
                    // the singular vertex definition (or the
                    // duplicated face information in the book ) so
                    // that only a tetrahedral number
                    // (i.e. (modes)*(modes+1)*(modes+2)/6) of modes
                    // are required consistent with the orthogonal
                    // basis.

                    // In the current structure the r index runs
                    // fastest rollowed by q and than the p index so
                    // that the matrix has a more compact structure

                    // Note that eModified_C is a re-organisation/
                    // duplication of eModified_B so will get a
                    // temporary Modified_B expansion and copy the
                    // correct components.

                    // Generate Modified_B basis;
                    BasisKey ModBKey(eModified_B,m_basisKey.GetNumModes(),
                                    m_basisKey.GetPointsKey());
                    BasisSharedPtr  ModB = BasisManager()[ModBKey];

                    Array<OneD, const NekDouble> ModB_data = ModB->GetBdata();

                    // Copy Modified_B basis into first
                    // (numModes*(numModes+1)/2)*numPoints entires of
                    // bdata.  This fills in the complete (r,p) face.

                    // Set up \phi^c_{p,q,r} = \phi^b_{p+q,r}

                    int N;
                    int B_offset = 0;
                    int offset = 0;
                    for(p = 0; p < numModes; ++p)
                    {
                        N = numPoints*(numModes-p)*(numModes-p+1)/2;
                        Vmath::Vcopy(N, &ModB_data[0]+B_offset,1,&m_bdata[0] + offset,1);
                        B_offset += numPoints*(numModes-p);
                        offset   += N;
                    }

                    // set up derivative of basis.
                    Blas::Dgemm('n','n',numPoints,
                                numModes*(numModes+1)*(numModes+2)/6,
                                numPoints,1.0,D,numPoints,
                                m_bdata.data(),numPoints,0.0,
                                m_dbdata.data(),numPoints);
                }
                break;

            case eGLL_Lagrange:
                {
                    mode = m_bdata.data();
                    boost::shared_ptr< Points<NekDouble> > m_points = PointsManager()[PointsKey(numModes, eGaussLobattoLegendre)];
                    const Array<OneD, const NekDouble>& zp(m_points->GetZ());

                    for (p=0; p<numModes; ++p, mode += numPoints)
                    {
                        for(q = 0; q < numPoints; ++q)
                        {
                            mode[q] = Polylib::hglj(p, z[q], zp.data(), numModes, 0.0, 0.0);
                        }
                    }

                    // define derivative basis
                    Blas::Dgemm('n', 'n', numPoints, numModes, numPoints, 1.0,
                                D, numPoints, m_bdata.data(), numPoints, 0.0,
                                m_dbdata.data(), numPoints);

                }//end scope
                break;
            case eGauss_Lagrange:
                {
                    mode = m_bdata.data();
                    boost::shared_ptr< Points<NekDouble> > m_points = PointsManager()[PointsKey(numModes, eGaussGaussLegendre)];
                    const Array<OneD, const NekDouble>& zp(m_points->GetZ());
					
                    for (p=0; p<numModes; ++p,mode += numPoints)
                    {
                        for(q = 0; q < numPoints; ++q)
                        {
                            mode[q] = Polylib::hgj(p, z[q], zp.data(), numModes, 0.0, 0.0);
                        }
                    }
					
                    // define derivative basis
                    Blas::Dgemm('n', 'n', numPoints, numModes, numPoints, 1.0,
                                D, numPoints, m_bdata.data(), numPoints, 0.0,
								m_dbdata.data(), numPoints);
					
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
                        m_bdata[(2*p+1)*numPoints+i] = -sin(p*M_PI*z[i]);

                        m_dbdata[ 2*p   *numPoints+i] = -p*M_PI*sin(p*M_PI*z[i]);
                        m_dbdata[(2*p+1)*numPoints+i] = -p*M_PI*cos(p*M_PI*z[i]);
                    }
                }

                break;
					
			
            // Fourier Single Mode (1st mode)
            case eFourierSingleMode:
					
                for(i = 0; i < numPoints; ++i)
                {
                    m_bdata[i] = cos(M_PI*z[i]);
                    m_bdata[numPoints+i] = -sin(M_PI*z[i]);
						
                    m_dbdata[i] = -M_PI*sin(M_PI*z[i]);
                    m_dbdata[numPoints+i] = -M_PI*cos(M_PI*z[i]);
                }
					
                for (p=1; p < numModes/2; ++p)
                {
                    for(i = 0; i < numPoints; ++i)
                    {
                        m_bdata[ 2*p   *numPoints+i] = 0.;
                        m_bdata[(2*p+1)*numPoints+i] = 0.;
							
                        m_dbdata[ 2*p   *numPoints+i] = 0.;
                        m_dbdata[(2*p+1)*numPoints+i] = 0.;
                    }
                }
                break;
					
					//Fourier Real Half Mode
					case eFourierHalfModeRe:
					
					m_bdata[0] = cos(M_PI*z[0]);
					m_dbdata[0] = -M_PI*sin(M_PI*z[0]);
					break;
					
					//Fourier Imaginary Half Mode
				    case eFourierHalfModeIm:
					
					m_bdata[0] = -sin(M_PI*z[0]);
					m_dbdata[0] = -M_PI*cos(M_PI*z[0]);
					break;

					
					
            case eChebyshev:
                {
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

                    // Define derivative basis
                    Blas::Dgemm('n', 'n', numPoints, numModes, numPoints, 1.0,
                                D, numPoints, m_bdata.data(), numPoints, 0.0,
                                m_dbdata.data(), numPoints);
                }
                break;

            case eMonomial:
                {
                    int P = numModes - 1;
                    NekDouble *mode = m_bdata.data();

                    for( int p = 0; p <= P; ++p, mode += numPoints )
                    {
                        for( int i = 0; i < numPoints; ++i )
                        {
                            mode[i] = pow(z[i], p);
                        }
                    }

                    // define derivative basis
                    Blas::Dgemm('n', 'n', numPoints, numModes, numPoints, 1.0,
                                D, numPoints, m_bdata.data(), numPoints, 0.0, 
                                m_dbdata.data(),numPoints);
                }//end scope
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
            case eGaussKronrodLegendre:
            case eGaussRadauKronrodMLegendre:
            case eGaussLobattoKronrodLegendre:
                returnval = (GetNumPoints() >= GetNumModes());
                break;

            default:
                break;
            }

            return returnval;
        }

        /** \brief Determine if basis has collocation property,
         *  i.e. GLL_Lagrange with Lobatto integration of appropriate order,
         *  Gauss_Lagrange with Gauss integration of appropriate order.
         */
        bool BasisKey::Collocation() const
        {
            return ((m_basistype     == eGLL_Lagrange         &&
                     GetPointsType() == eGaussLobattoLegendre &&
                     GetNumModes()   == GetNumPoints())       || 
                    (m_basistype     == eGauss_Lagrange       &&
                     GetPointsType() == eGaussGaussLegendre   &&
                     GetNumModes()   == GetNumPoints()));
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

