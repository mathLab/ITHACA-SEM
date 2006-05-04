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
#include <LibUtilities/ErrorUtil.hpp>

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

      /// \brief return order of basis
      inline int GetBasisOrder() const
      {
	return m_basisorder;
      }
    
      /// \brief return points order at which  basis is defined 
      inline int GetPointsOrder() const
      {
	return m_pointsorder;
      }

      /// \brief return type of expansion basis
      inline BasisType GetBasisType() const
      {
	return m_basistype;
      }

      /// \brief return type of quadrature 
      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }    
    
      /// \brief return alpha weight of Guassian quadrature 
      inline double GetAlpha() const
      {
	return m_alpha;
      }

      /// \brief return beta weight of Guassian quadrature 
      inline double GetBeta() const
      {
	return m_beta;
      }

      /// \brief Check to see if the quadrature of expansions x is the
      /// same as the calling basis
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
    
      /// \brief Check to see if basis expansions x is the same as the
      /// calling basis
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
    
      /// \brief determine if basis definition has exact integration for
      /// inner product
      int  ExactIprodInt() const;

      /// \brief Determine if basis has collocation properties,
      /// i.e. GLL_Lagrange with appropriate quadrature
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
      
      /// All BasisKeys  share the same PolyManager
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

      /// \brief return basis definition array m_bdata
      inline const double * GetBdata() const
      {
	return m_bdata;
      }

      /// \brief return basis definition array m_dbdata
      inline const double * GetDbdata() const
      {
	return m_dbdata;
      }

    protected:
      double    *m_bdata;       //!< Basis definition
      double    *m_dbdata;      //!< Derivative Basis definition

      /// wrapper around function BasisMem calling private data
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

      // Method used to generate appropriate basis
      void GenBasis();

    private:
      
  };
    
  } // end of namespace
} // end of namespace

#endif //STDBASIS_H


