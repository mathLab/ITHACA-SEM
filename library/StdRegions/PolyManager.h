///////////////////////////////////////////////////////////////////////////////
//
// File PolyManager.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Header for Polynomial routines definition and manager
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POLYMANGER_H
#define POLYMANGER_H

#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdBasis.h>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
  namespace StdRegions
  {

    class BasisKey;
    class PolyManager;
    class ZeroWeightsKey;
    class DerivMatrixKey;
    class InterpMatrixKey;
    class ZeroWeights;
    class DerivMatrix;
    class InterpMatrix;

    /* *********************
       ZeroWeights Class
       ********************* */

    class ZeroWeightsKey
    {
    public:
      ZeroWeightsKey(PointsType ptype, int np, double alpha, double beta):
	m_pointstype  (ptype),
	m_pointsorder (np),
	m_alpha       (alpha),
	m_beta        (beta)
      {
      }

      ~ZeroWeightsKey()
      {
      }

      inline int GetPointsOrder() const
      {
	return m_pointsorder;
      }

      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }

      inline double GetAlpha() const
      {
	return m_alpha;
      }

      inline double GetBeta() const
      {
	return m_beta;
      }

      // Overloaded Operators

      friend bool operator  == (const ZeroWeightsKey& x, const ZeroWeights& y);
      friend bool operator  != (const ZeroWeightsKey& x, const ZeroWeights& y);

      friend bool operator  == (const ZeroWeights& x, const ZeroWeightsKey& y);
      friend bool operator  != (const ZeroWeights& x, const ZeroWeightsKey& y);

      friend bool operator  == (const ZeroWeightsKey& x, const ZeroWeights* y);
      friend bool operator  != (const ZeroWeightsKey& x, const ZeroWeights* y);

      friend bool operator  == (const ZeroWeights* x, const ZeroWeightsKey& y);
      friend bool operator  != (const ZeroWeights* x, const ZeroWeightsKey& y);

    private:
      double m_alpha, m_beta;
      int    m_pointsorder;
      PointsType m_pointstype;
    };


    class ZeroWeights
    {
    public:

      ZeroWeights(const PointsType ptype, const int np,
          const double alpha, const double beta)
      {
	SetZW(ptype,np,alpha,beta);
      }

      ZeroWeights(const ZeroWeightsKey &Z)
      {
	SetZW(Z.GetPointsType(),Z.GetPointsOrder(),Z.GetAlpha(),Z.GetBeta());
      }

      ZeroWeights(const ZeroWeights &Z)
      {
	m_pointsorder = Z.m_pointsorder;
	m_alpha = Z.m_alpha;
	m_beta = Z.m_beta;

	if(m_pointsorder > 0)
	{
	  m_z = new double[m_pointsorder];
	  m_w = new double[m_pointsorder];
	  
	  for(int i=0;i<m_pointsorder;i++)
	  {
	    m_z[i] = Z.m_z[i];
	    m_w[i] = Z.m_w[i];
	  }
	}
	else
	{
	  m_z = (double *)NULL;
	  m_w = (double *)NULL;
	}
      }

      // default destructor()
      ~ZeroWeights()
      {
	if(m_pointsorder+1)//use initialised value as trip
	{
	  if(m_z)
	  {
	    delete [] m_z;
	  }

	  if(m_w)
	  {
	    delete [] m_w;
	  }
	}
	m_z = (double*)NULL;
	m_w = (double*)NULL;
      }

      inline int GetPointsOrder() const
      {
	return m_pointsorder;
      }

      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }

      inline double GetAlpha() const
      {
	return m_alpha;
      }

      inline double GetBeta() const
      {
	return m_beta;
      }

      inline double *GetZ() const
      {
	return m_z;
      };

      inline double *GetW() const
      {
	return m_w;
      };

      // Overloaded Operators
      friend bool operator  == (const ZeroWeightsKey& x, const ZeroWeights& y);
      friend bool operator  != (const ZeroWeightsKey& x, const ZeroWeights& y);

      friend bool operator  == (const ZeroWeights& x, const ZeroWeightsKey& y);
      friend bool operator  != (const ZeroWeights& x, const ZeroWeightsKey& y);

      friend bool operator  == (const ZeroWeightsKey& x, const ZeroWeights* y);
      friend bool operator  != (const ZeroWeightsKey& x, const ZeroWeights* y);

      friend bool operator  == (const ZeroWeights* x, const ZeroWeightsKey& y);
      friend bool operator  != (const ZeroWeights* x, const ZeroWeightsKey& y);

    private:
      double *m_z;
      double *m_w;
      double  m_alpha, m_beta;
      int     m_pointsorder;
      PointsType m_pointstype;

      void SetZW(const PointsType ptype, const int np, const double alpha,
         const double beta);
    };


    /**********************
     DerivMatrix Class
     ********************* */

    class DerivMatrixKey
    {
    public:

      DerivMatrixKey(const PointsType ptype, const int np, const double alpha,
             const double beta):
	m_pointstype (ptype),
	m_pointsorder (np),
	m_alpha (alpha),
	m_beta  (beta)
      {
      }

      // default destructor()
      ~DerivMatrixKey()
      {
      }

      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }

      inline int GetPointsOrder() const
      {
	return m_pointsorder;
      }

      inline double GetAlpha() const
      {
	return m_alpha;
      }

      inline double GetBeta() const
      {
	return m_beta;
      }

      // Overloaded Operators

      friend bool operator  == (const DerivMatrixKey& x, const DerivMatrix& y);
      friend bool operator  != (const DerivMatrixKey& x, const DerivMatrix& y);

      friend bool operator  == (const DerivMatrix& x, const DerivMatrixKey& y);
      friend bool operator  != (const DerivMatrix& x, const DerivMatrixKey& y);

      friend bool operator  == (const DerivMatrixKey& x, const DerivMatrix* y);
      friend bool operator  != (const DerivMatrixKey& x, const DerivMatrix* y);

      friend bool operator  == (const DerivMatrix* x, const DerivMatrixKey& y);
      friend bool operator  != (const DerivMatrix* x, const DerivMatrixKey& y);

    private:
      double m_alpha;
      double m_beta;
      int    m_pointsorder;
      PointsType m_pointstype;
    };


    class DerivMatrix
    {
    public:

      DerivMatrix(PolyManager *P, const PointsType ptype, const int np,
          const double alpha, const double beta)
      {
	SetD(P,ptype,np,alpha,beta);
      }

      DerivMatrix(PolyManager *P, const DerivMatrixKey &D)
      {
	SetD(P,D.GetPointsType(),D.GetPointsOrder(),D.GetAlpha(),D.GetBeta());
      }

      DerivMatrix(const DerivMatrix &D)
      {
	m_pointsorder   = D.m_pointsorder;
	m_alpha = D.m_alpha;
	m_beta  = D.m_beta;
	
	if(m_pointsorder > 0)
	{
	  m_dm = new double[m_pointsorder*m_pointsorder];
	  for(int i=0;i<m_pointsorder*m_pointsorder;i++)
	  {
	    m_dm[i] = D.m_dm[i];
	  }
	}
      }

      // default destructor()
      ~DerivMatrix()
      {
	if((m_pointsorder+1) && m_dm )//use initialised value as trip
        {
	  delete [] m_dm;
	}

	m_dm = (double*) NULL;
      }

      inline PointsType GetPointsType() const
      {
	return m_pointstype;
      }

      inline int GetPointsOrder() const
      {
	return m_pointsorder;
      }

      inline double GetAlpha() const
      {
	return m_alpha;
      }

      inline double GetBeta() const
      {
	return m_beta;
      }

      inline double *GetD() const
      {
	return m_dm;
      };

      // Overloaded Operators

      friend bool operator  == (const DerivMatrixKey& x, const DerivMatrix& y);
      friend bool operator  != (const DerivMatrixKey& x, const DerivMatrix& y);

      friend bool operator  == (const DerivMatrix& x, const DerivMatrixKey& y);
      friend bool operator  != (const DerivMatrix& x, const DerivMatrixKey& y);

      friend bool operator  == (const DerivMatrixKey& x, const DerivMatrix* y);
      friend bool operator  != (const DerivMatrixKey& x, const DerivMatrix* y);

      friend bool operator  == (const DerivMatrix* x, const DerivMatrixKey& y);
      friend bool operator  != (const DerivMatrix* x, const DerivMatrixKey& y);

    private:
      double *m_dm;
      double m_alpha, m_beta;
      int    m_pointsorder;
      PointsType m_pointstype;

      void SetD(PolyManager *P, const PointsType ptype, const int np, const
        double alpha, const double beta);
    };


    /* *********************
       InterpMatrix Class
       ********************* */

    class InterpMatrixKey
    {

    public:
      InterpMatrixKey(const PointsType orig_ptype, const int orig_np,
              const double orig_alpha,    const double orig_beta,
              const PointsType dest_ptype, const int dest_np,
              const double dest_alpha,    const double dest_beta):
	m_orig_np    (orig_np),
	m_dest_np    (dest_np),
	m_orig_ptype (orig_ptype),
	m_dest_ptype (dest_ptype),
	m_orig_alpha (orig_alpha),
	m_dest_alpha (dest_alpha),
	m_orig_beta  (orig_beta),
	m_dest_beta  (dest_beta)
      {
      }

      // default destructor()
      ~InterpMatrixKey()
      {
      }

      inline int GetOrigNp() const
      {
	return m_orig_np;
      }

      inline PointsType GetOrigPtype() const
      {
	return m_orig_ptype;
      }

      inline double GetOrigAlpha() const
      {
	return m_orig_alpha;
      }

      inline double GetOrigBeta() const
      {
	return m_orig_beta;
      }

      inline int GetDestNp() const
      {
	return m_dest_np;
      }

      inline PointsType GetDestPtype() const
      {
	return m_dest_ptype;
      }

      inline double GetDestAlpha() const
      {
	return m_dest_alpha;
      }

      inline double GetDestBeta() const
      {
	return m_dest_beta;
      }

      // Overloaded Operators

      friend bool operator  == (const InterpMatrixKey& x, const InterpMatrix& y);
      friend bool operator  != (const InterpMatrixKey& x, const InterpMatrix& y);

      friend bool operator  == (const InterpMatrix& x, const InterpMatrixKey& y);
      friend bool operator  != (const InterpMatrix& x, const InterpMatrixKey& y);

      friend bool operator  == (const InterpMatrixKey& x, const InterpMatrix* y);
      friend bool operator  != (const InterpMatrixKey& x, const InterpMatrix* y);

      friend bool operator  == (const InterpMatrix* x, const InterpMatrixKey& y);
      friend bool operator  != (const InterpMatrix* x, const InterpMatrixKey& y);


    private:
      double    m_orig_alpha;
      double    m_orig_beta;
      double    m_dest_alpha;
      double    m_dest_beta;
      int       m_orig_np;
      int       m_dest_np;
      PointsType m_orig_ptype;
      PointsType m_dest_ptype;
    };


    class InterpMatrix
    {
    public:

      InterpMatrix(PolyManager *P, const InterpMatrixKey &I)
      {
	SetI(P,I.GetOrigPtype(),I.GetOrigNp(),I.GetOrigAlpha(),I.GetOrigBeta(),
	     I.GetDestPtype(),I.GetDestNp(),I.GetDestAlpha(),I.GetDestBeta());
      }

      InterpMatrix(PolyManager *P,const PointsType orig_ptype, const int orig_np,
           const double orig_alpha, const double orig_beta,
           PointsType dest_ePtype, const int dest_np,
           const double dest_alpha, const double dest_beta)
      {
	SetI(P, orig_ptype, orig_np, orig_alpha, orig_beta,
	     dest_ePtype, dest_np, dest_alpha, dest_beta);
      }

      InterpMatrix(const InterpMatrix &I)
      {
	m_orig_np    = I.m_orig_np;
	m_dest_np    = I.m_dest_np;
	m_orig_ptype = I.m_orig_ptype;
	m_dest_ptype = I.m_dest_ptype;
	m_orig_alpha = I.m_orig_alpha;
	m_dest_alpha = I.m_dest_alpha;
	m_orig_beta  = I.m_orig_beta;
	m_dest_beta  = I.m_dest_beta;

	if((m_orig_np > 0)&&(m_dest_np > 0))
	{
	  m_im = new double[m_orig_np*m_dest_np];
	  
	  for(int i=0;i<m_orig_np*m_dest_np;i++)
	  {
	    m_im[i] = I.m_im[i];
	  }
	}
      }

      // default destructor()
      ~InterpMatrix()
      {
	if((m_orig_np+1) && m_im)//use initialised value as trip
	{
	  delete [] m_im;
	}
	m_im = (double*) NULL;
      }

      inline int GetOrigNp() const
      {
	return m_orig_np;
      }

      inline PointsType GetOrigPtype() const
      {
	return m_orig_ptype;
      }

      inline double GetOrigAlpha() const
      {
	return m_orig_alpha;
      }

      inline double GetOrigBeta() const
      {
	return m_orig_beta;
      }

      inline int GetDestNp() const
      {
	return m_dest_np;
      }

      inline PointsType GetDestPtype() const
      {
	return m_dest_ptype;
      }

      inline double GetDestAlpha() const
      {
	return m_dest_alpha;
      }

      inline double GetDestBeta() const
      {
	return m_dest_beta;
      }

      double * GetInterp() const
      {
	return m_im;
      };

      // Overloaded Operators

      friend bool operator  == (const InterpMatrixKey& x, const InterpMatrix& y);
      friend bool operator  != (const InterpMatrixKey& x, const InterpMatrix& y);

      friend bool operator  == (const InterpMatrix& x, const InterpMatrixKey& y);
      friend bool operator  != (const InterpMatrix& x, const InterpMatrixKey& y);

      friend bool operator  == (const InterpMatrixKey& x, const InterpMatrix* y);
      friend bool operator  != (const InterpMatrixKey& x, const InterpMatrix* y);

      friend bool operator  == (const InterpMatrix* x, const InterpMatrixKey& y);
      friend bool operator  != (const InterpMatrix* x, const InterpMatrixKey& y);

    private:
      double *m_im;
      double  m_orig_alpha;
      double  m_orig_beta;
      double  m_dest_alpha;
      double  m_dest_beta;
      int     m_orig_np;
      int     m_dest_np;
      PointsType m_orig_ptype;
      PointsType m_dest_ptype;

      void SetI(PolyManager *P, const PointsType orig_ptype, const int orig_np,
        const double orig_alpha, const double orig_beta,
        const PointsType dest_ptype, const int dest_np,
        const double dest_alpha,const double dest_beta);

    };

    class PolyManager
    {

    public:
      PolyManager()
      {
	m_zwvals = new std::vector<ZeroWeights*>[SIZE_PointsType];
	m_dmvals = new std::vector<DerivMatrix*>[SIZE_PointsType];
	m_imvals = new std::vector<InterpMatrix*>[SIZE_PointsType];

	for(int i=0;i<SIZE_PointsType;i++)
	{
	  m_zwvalscur.push_back(m_zwvals[i].begin());
	  m_dmvalscur.push_back(m_dmvals[i].begin());
	  m_imvalscur.push_back(m_imvals[i].begin());
	}
      }

      ~PolyManager()
      {
	std::vector<ZeroWeights*>::iterator defzw;
	std::vector<DerivMatrix*>::iterator defder;
	std::vector<InterpMatrix*>::iterator definterp;

	for(int i=0;i<SIZE_PointsType;i++)
	{
	  for(defzw = m_zwvals[i].begin(); defzw != m_zwvals[i].end(); ++defzw)
	  {
	    delete defzw[0];
	  }

	  for(defder = m_dmvals[i].begin(); defder != m_dmvals[i].end();
	      ++defder)
	  {
	    delete defder[0];
	  }

	  for(definterp = m_imvals[i].begin(); definterp != m_imvals[i].end();
	      ++definterp)
	  {
	    delete definterp[0];
	  }
	}
	
	delete[] m_zwvals;
	delete[] m_dmvals;
	delete[] m_imvals;
      };
      
      void GetZW(const PointsType ptype, const int  np, const double * &z,
		 const double * &w, const double alpha, const double beta);
      void GetZW(const BasisKey *B, const double * &z,  const double * &w);
      void ResetZW(void);

      
      void GetD(const PointsType ptype, const int  np, const double * &D,
		const double alpha, const double beta);
      void GetD(const BasisKey *B, const double * &D);
      void ResetD(void);

      void GetI(const PointsType orig_ptype, const int orig_np,
        const double orig_alpha, const double orig_beta,
        const PointsType dest_ptype, const int dest_np,
        const double dest_alpha, const double dest_beta,
        const double * &I);

      void GetI(const BasisKey *orig_B,
        const BasisKey *dest_B, const double * &I);
      void ResetI(void);

      void GetInterpVec(const double zi, const PointsType ptype,
            const double *z, const int np, const double alpha,
            const double beta, double *I);
    private:

      std::vector<ZeroWeights*> *m_zwvals;
      std::vector<std::vector<ZeroWeights*>::iterator> m_zwvalscur;

      std::vector<DerivMatrix*> *m_dmvals;
      std::vector<std::vector<DerivMatrix*>::iterator> m_dmvalscur;

      std::vector<InterpMatrix*> *m_imvals;
      std::vector<std::vector<InterpMatrix*>::iterator> m_imvalscur;

      void ZWLookup(const PointsType ptype, const  int np,
            const double* &z, const double* &w,
            const double alpha, const double beta,
            std::vector<ZeroWeights*> &val);
      void DLookup (const PointsType ptype, const  int np,
            const double* &D,const double alpha,
            const double beta, std::vector<DerivMatrix*> &val);
      void ILookup (const PointsType orig_ptype, const int orig_np,
            const double orig_alpha, const double orig_beta,
            const PointsType dest_ptype, const int dest_np,
            const double dest_alpha, const double dest_beta,
            const double* &I, std::vector<InterpMatrix*> &val);
    };

  } // end of namespace
} // end of namespace

#endif //POLYMANGER_H

/**
 * $Log: PolyManager.h,v $
 * Revision 1.1  2006/05/04 18:58:30  kirby
 * *** empty log message ***
 *
 * Revision 1.32  2006/03/21 09:21:31  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.31  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.30  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.29  2006/03/01 17:07:33  sherwin
 *
 * Added new location of polylib in header definitions
 *
 * Revision 1.28  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 * Revision 1.27  2006/02/26 23:37:29  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 * Revision 1.26  2006/02/19 13:26:13  sherwin
 *
 * Coding standard revisions so that libraries compile
 *
 * Revision 1.25  2006/02/15 08:06:35  sherwin
 *
 * Put files into coding standard (although they do not compile)
 *
 *
 **/
