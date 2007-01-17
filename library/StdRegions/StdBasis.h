///////////////////////////////////////////////////////////////////////////////
//
// File StdBasis.h
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
// Description: Header file of Basis definition 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDBASIS_H
#define STDBASIS_H

#include <math.h>

#include <StdRegions/PolyManager.h>
#include <loki/Factory.h>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
  namespace StdRegions 
  {
  

    class PolyManager;

    class BasisKey
    {
    
    public:
      BasisKey()
      {
	    m_basisorder = 0;
	    m_pointsorder = 0;
      }
    
      BasisKey(const BasisType btype, const int border, const PointsType ptstype,
	       const int ptsorder,  const double alpha, const double beta):
	m_basistype(btype),
	m_basisorder(border),
	m_pointstype(ptstype),
	m_pointsorder(ptsorder),
	m_alpha(alpha),
	m_beta (beta)
      {
      };
    
      BasisKey(const BasisKey &B):
	m_basistype   (B.m_basistype),
	m_basisorder  (B.m_basisorder),
	m_pointstype  (B.m_pointstype),
	m_pointsorder (B.m_pointsorder),
	m_alpha       (B.m_alpha),
	m_beta        (B.m_beta)
	{
	}
      
      ~BasisKey()
       {
       };

      /** \brief return order of basis */
      inline int GetBasisOrder() const
      {
	return m_basisorder;
      }
    
      /** \brief return points order at which  basis is defined */
      inline int GetPointsOrder() const
      {
		return m_pointsorder;
      }

      /** \brief return type of expansion basis */
      inline BasisType GetBasisType() const
      {
	return m_basistype;
      }

      /** \brief return type of quadrature */
      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }    
    
      /** \brief return alpha weight of Guassian quadrature */
      inline double GetAlpha() const
      {
	return m_alpha;
      }

      /** \brief return beta weight of Guassian quadrature */
      inline double GetBeta() const
      {
	return m_beta;
      }

      /** \brief Check to see if the quadrature of expansions x is the
       *  same as the calling basis
       */
      inline int SamePoints(const BasisKey &x) const
      {
	if((x.m_pointsorder == m_pointsorder)&&(x.m_pointstype == m_pointstype)
	   && (x.m_alpha == m_alpha)&&(x.m_beta == m_beta))
	{
	  return true;
	}
	else
	{
	  return false;
	}
      }
    
      /** \brief Check to see if basis expansions x is the same as the
       *  calling basis
       */
      inline int SameExp(const BasisKey &x) const
      {
	if((x.m_basisorder == m_basisorder)&&(x.m_basistype == m_basistype))
	{
	  return true;
	}
	else
	{
	  return false;
	}
      }
    
      /** \brief determine if basis definition has exact integration for
       *  inner product
       */
      int  ExactIprodInt() const;

      /** \brief Determine if basis has collocation properties,
       *  i.e. GLL_Lagrange with appropriate quadrature
       */
      int  Collocation() const;

      void GetInterpVec(const double zi, double *I) const;    

      //Overloaded Operators
      friend bool operator  == (const BasisKey& x, const BasisKey& y);
      friend bool operator  == (const BasisKey* x, const BasisKey& y);
      friend bool operator  == (const BasisKey& x, const BasisKey *y);
      friend bool operator  != (const BasisKey& x, const BasisKey& y);
      friend bool operator  != (const BasisKey* x, const BasisKey& y);
      friend bool operator  != (const BasisKey& x, const BasisKey *y);
    
    protected:
      int        m_basisorder;  //< Expansion Order
      BasisType  m_basistype;   //< Expansion Type 
      int        m_pointsorder; //< Number of discrete evaluation points 
      PointsType m_pointstype;  //< Quadrature/points type
      double     m_alpha;       //< Quadrature alpha
      double     m_beta;        //< Quadrature beta
      
      /** All BasisKeys  share the same PolyManager */
      typedef Loki::SingletonHolder<PolyManager> PolyManagerSingleton;
      
    private:
    };
  
    ///////////////////////////////////////////////////////////////////////////
    class Basis: public BasisKey
    {
    
    public:    
      
      Basis()
      {
	    ASSERTL0(false,
			 "Error in Basis.C: Default Constructor instantiation "
			 "not allowed."); 
      }

      
      Basis(BasisType btype,int order, PointsType ptype, 
	    const int nq, const double alpha, 
	    const double beta):
	BasisKey(btype,order,ptype,nq,alpha,beta)
      {
	Basis::Initialize();
      };
    
      Basis(const BasisKey &B):
	BasisKey(B)
      {
	Basis::Initialize();
      };


      //Copy Constructor
      Basis(const Basis &B)
      {
	m_basistype   = B.m_basistype;
	m_basisorder  = B.m_basisorder;
	m_pointstype  = B.m_pointstype;
	m_pointsorder = B.m_pointsorder;
	m_alpha       = B.m_alpha;
	m_beta        = B.m_beta;

	int size = BasisMem();

	// Allocate Memory
	m_bdata  = new double [size];
	m_dbdata = new double [size];
	
	for(int i=0;i<size;i++)
	{
	  m_bdata[i]  = B.m_bdata[i];
	  m_dbdata[i] = B.m_dbdata[i];
	}
      };
      
      // default destructor()
      ~Basis()
      {
	if(m_bdata)
	{
	  delete [] m_bdata;  
	}

	if(m_dbdata)
	{
	  delete [] m_dbdata;  
	}
      
	m_bdata  = NULL;
	m_dbdata = NULL;
      };
    
      void ResetBasisOrder(int order)
      {
	m_basisorder = order;
	
	if(m_bdata)
	{
	  delete[] m_bdata;  //delete old space
	  m_bdata = NULL;
	}

	if(m_dbdata)
	{
	  delete[] m_dbdata;  //delete old space
	  m_bdata = NULL;
	}
      
	m_bdata  = new double[BasisMem()]; //allocate new space
	m_dbdata = new double[BasisMem()]; //allocate new space
	
	GenBasis();
      }

      /** \brief return basis definition array m_bdata */
      inline const double * GetBdata() const
      {
	return m_bdata;
      }

      /** \brief return basis definition array m_dbdata */
      inline const double * GetDbdata() const
      {
	return m_dbdata;
      }

    protected:
      double    *m_bdata;       /**< Basis definition */
      double    *m_dbdata;      /**< Derivative Basis definition */

      /** \brief wrapper around function BasisMem calling private data */
      int BasisMem()
      {
	switch(m_basistype)
	{
	case eOrtho_B: case eModified_B: case eModified_C:
	  return  m_basisorder*(m_basisorder+1)/2*m_pointsorder;
	  break;
	case eOrtho_C:
	  return  m_basisorder*(m_basisorder+1)*
	         (m_basisorder+2)/6*m_pointsorder;
	  break;
	default:
	  return m_basisorder*m_pointsorder;
	}
    }  
      
      void Initialize()
      {
	// Allocate Memory
	int size = BasisMem();
	m_bdata  = new double [size];
	m_dbdata = new double [size];
	
	GenBasis();
      };


      /** \brief Generate a basis array stored in \a m_bdata[m][i] of type 
       *  \a init_ptype of order \a init_order at \a init_zorder points
       *  distributed according to \a init_z[i].
       *
       *  \b init_pointstype can take the following arguments: 
       *
       *  \li \b eOrtho_A, \b eLegendre: Orthogonal cardinal function A 
       *  (Legendre polynomials) \f$\widetilde{\psi}^a_p(z) = L_p(z)\f$\n\n
       *  \f$ \mbox{\bf m\_bdata}[p][i] = \widetilde{\psi}^a_p(z_i) 
       *  = P^{0,0}_p(z_i) \f$ 
       *  normalised so that 
       *  \f$ (\widetilde{\psi}^a_{p},\widetilde{\psi}^a_{p})=1\f$
       *
       *  \li \b eOrtho_B: Orthogonal cardinal function B 
       *  \f$\widetilde{\psi}^b_{pq}(z)\f$\n\n  
       *  \f$ \mbox{\bf m\_bdata}[m(p,q)][i]=\widetilde{\psi}^b_{pq}(z_i) = 
       *  [(1-z_i)]^p P^{2p+1,0}_q(z_i) \f$
       *  normalised so that 
       *  \f$(\widetilde{\psi}^b_{pq},\widetilde{\psi}^b_{pq})=1\f$\n\n
       *  where \f$ m(p,q)=\frac{p}{2}(2\ \mbox{init\_order}-p+1)+q\f$ 
       *  and  \f$0 \leq p, p+q < \mbox{init\_order}, \f$
       * 
       *  \li \b eOrtho_C: Orthogonal cadinal function C:
       *  \f$ \widetilde{\psi}^b_{pqr}(z)\f$ \n
       *  \f$ \widetilde{\psi}^c_{pqr}(z_i)
       *  =  [(1-z_i)]^{p+q} P^{2p+2q+1,0}_q(z_i) \f$
       *  normalised so that 
       *  \f$(\widetilde{\psi}^c_{pqr},\widetilde{\psi}^c_{pqr})=1\f$\n
       *  Due symmetry we note that \n
       *  \f$ \widetilde{\psi}^c_{pqr}(z_i) = 
       *  \widetilde{\psi}^{c\star}_{(p+q)r}(z_i) = 
       *  \widetilde{\psi}^{c\star}_{sr}(z_i) \f$ \n
       *  and so we only need calculate \n\n
       *  \f$ \mbox{\bf m\_bdata}[m(s,r)][i]=
       *  \widetilde{\psi}^{c\star}_{sr}(z_i)	\
       *  = [(1-z)]^{s} P^{2s+2,0}_q(z)\f$
       *  (Which is very similar to the cardinal function B)\n\n
       *  where \f$ s = p+q, m(s,r)=\frac{s}{2}
       *  (2\ \mbox{init\_order}-s+1)+r\f$		    
       *  and  \f$0 \leq s, s+r < \mbox{init\_order}, \f$
       *
       *  \li \b eModified_A:  Modified cardial functions A, 
       *  \f$ \psi^a_p(z)\f$\n\n
       *  \f$ \mbox{\bf m\_bdata}[p][i] = \left \{ \begin{array}{ll}
       *  (1-z_i)/2    & p=0          \\          
       *  (1+z_i)/2    & p=1          \\         
       *  (1-z_i)/2 (1+z_i)/2 & p=2   \\
       *  (1-z_i)/2 (1+z_i)/2 P_1^{1,1}(z_i)  & p=3\\
       *  \vdots &   \vdots \\
       *  (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & p=j \end{array} \right . 
       *  \f$\n\n
       *  Note this has been changed from %Nektar where (1+z)/2 was
       *  stored first in list rather than (1-z)/2. The new format is
       *  more consistent with a left to right convention of
       *  interpreting the basis index \a p.\n
       *  The vertex modes are ordered in a hierarchical fashion rather
       *  than the nodal format as used in the book.
       *
       *  \li \b eModified_B, \b eModified_C Modified cardinal functions B and C
       *  \f$ \psi^b_{pq}(z)\f$, \f$ \psi^c_{pqr}(z)\f$\n
       *  \b Note: eModified_B and eModified_C can use the same definition
       *  since \f$ \psi^c_{pqr}(z) =  \psi^b_{(p+q)r}(z) \f$\n 
       *  \f$ \mbox{\bf m\_bdata[m(p,q)][i]} = 
       *  \left \{ \begin{array}{ll}
       *  (1-z_i)/2            & m = 0,  (p=0,q=0) \\
       *  (1+z_i)/2            & m = 1,  (p=0,q=1) \\
       *  (1-z_i)/2 (1+z_i)/2  & m = 2,  (p=0,q=2) \\
       *  \vdots       &       \vdots \\
       *  (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & m = j, (p=0,q=j) \\  
       *  & \\
       *  (1-z_i)/2            & m = \mbox{order},  (p=1,q=0) \\
       *  (1-z_i)/2 (1+z_i)/2  & m = \mbox{order}+1,  (p=1,q=2) \\
       *  \vdots       &       \vdots    \\
       *  (1-z_i)/2 (1+z_i)/2 P_{j-2}^{1,1}(z_i) & m = \mbox{order}+j,
       *  (p=1,q=j) \\
       *  & \\
       *  \left [(1-z_i)/2\right]^2  & m=\mbox{2 order-1}, (p=2,q=0) \\
       *  \left [(1-z_i)/2\right]^2 (1+z_i)/2 & m=\mbox{2 order}, (p=2,q=1) \\
       *  \vdots	           &         \vdots              \\
       *  \left [(1-z_i)/2\right]^2 (1+z_i)/2 P_{j-1}^{3,1}(z_i) &  
       *   m=\mbox{2 order}+j-1,  (p=2,q=j) \\
       *   & \\
       *   \left [(1-z_i)/2\right]^{k+1} & m(p,q), (p=3,q=0) \\
       *  \left [(1-z_i)/2\right]^{k+1}(1+z_i)/2 & m(p,q), (p=3.q=1) \\
       *  \vdots	              &         \vdots         \\
       *  \left [(1-z_i)/2\right]^{k+1} (1+z_i)/2 P_{j-1}^{2k+1,1}(z_i) & 
       *  m(p,q), (p=k,q=j) 
       *  \end{array} \right . \f$ \n \n
       *  Note: this has been changed from %Nektar where (1+z)/2 was
       *  stored first in list rather than (1-z)/2. The new format is
       *  more consistent with a left to right convention. The modes
       *  used in edge 0 are now also interlaced with the interior modes
       *  which is convenient for matrix operations in the  
       *  sumfactorizations\n\n
       *  The vertex modes are ordered in a hierarchical fashion
       *  rather than the nodal format as used in the book.\n\n
       *  Finally the second row which is almost identical to the first 
       *  row but does not include the mode (1+z)/2. This is to make the block 
       *  consistent with the size of the orthonormal basis and means that
       *  similar sumfactorisation routines can be called for both basis althoug
       *  the top vertex mode will remain a special case in the Modified_B
       *  bases.
       *
       *  \li \b eFourier: Fourier bases stored in real and then imaginary 
       *  parts, in the interval [-1,1] i.e.\n \n
       *  \f$ \mbox{\bf m\_bdata}[p][i] = \left \{ \begin{array}{ll}
       *  1 & p=0 \\ 
       *  \sin((m+1) \pi z) & p=1 \\ 
       *  \cos(\pi z_i) & p=2 \\
       *  \sin(\pi z_i) & p=3 \\
       *  \cos(2 \pi z_i) & p=4 \\
       *  \sin(2 \pi z_i) & p=5 \\
       *  \vdots & \vdots \\
       *  \cos(m \pi z) & p = 2m \\ 
       *  \sin(m \pi z) & p = 2m+1  \end{array} \right . \f$\n \n
       *  Note: Should alway be called with a factor of 2 modes
       *
       *  \li \b eChebychev: Chebychev polynomials \f$ T_p(z) \f$\n\n
       *  \f$\mbox{\bf m\_bdata}[p][i] = T_p(z_i) = 
       *  2^{2p} (p!)^2/(2n!) P^{-1/2,-1/2}_p(z_i) \f$\n
       */
      void GenBasis();

    private:
      
  };
    
  } // end of namespace
} // end of namespace

#endif //STDBASIS_H


