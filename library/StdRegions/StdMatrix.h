///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrix.h
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
// Description: Standard Matrix container details
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDMATRIX_H
#define STDMATRIX_H

//Standard Matrix
#include <math.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/PolyManager.h>
#include <StdRegions/StdExpansion.h>

#include <stdio.h>

namespace Nektar
{
    namespace StdRegions
    {
	
	class StdMatContainer
	{
	public:
	    
	    StdMatContainer(double * M):
		m_factored (false),
		m_lda      (0),
		m_bwidth   (0),
		m_ldiag    (0),
		m_packed_matrix ((double *)NULL)
	    {
		m_matrix = M;
	    }
	    
	    StdMatContainer(StdExpansion * E, double * M):
		m_factored (false),
		m_lda      (0),
		m_bwidth   (0),
		m_ldiag    (0),
		m_packed_matrix ((double *)NULL)
	    {
		m_numbases = E->GetNumBases();
		m_basisorder = new int[m_numbases];
		m_base = new const Basis*[m_numbases];
		
		for(int i=0;i<m_numbases;i++)
		{
		    m_basisorder[i] = E->GetBasisOrder(i);
		    m_base[i]     = E->GetBasis(i);
		}
		m_matrix = M;
	    }

	    
	    ~StdMatContainer()
	    {
		 
		if(m_numbases)
		{
		    if(m_basisorder)
		    {
			delete[] m_basisorder;
			m_basisorder = NULL;
		    }
		    
		    if(m_base)
		    {
			delete[] m_base;
			m_base = NULL;
		    }
		}
		
		// It is assumed that StdMatContainer always owns the vector
		// attached to _matrix, and hence it is responsible for
		// deleting the space
		
		if(m_matrix)
		{
		    delete[] m_matrix;
		}
		
		if(m_packed_matrix)
		{
		    delete[] m_packed_matrix;
		    m_packed_matrix = NULL;
		}
	    }

      inline int GetNumBases() const
      {
	return m_numbases;
      }

      inline int GetBasisOrder(int index) const
      {
	return(m_basisorder[index]);
      }

      inline const Basis * GetBasis(int index) const
      {
	return m_base[index];
      }

      inline double * GetMatrix() const
      {
	return m_matrix;
      }

      /** \brief Forward Multiplication
       *
       *  Multiply full matrix so that \a out = \a _matrix x \a in 
       */
      void Mxv(const double *in, double *out);

      // Matrix factorisation initialisation routine

      //  Inversion routines

      /** \brief declare matrix memory */
      void SetMemPackedMatrix();

      /** \brief factorise matrix depending upon definition stored in 
       *  m_matform and store in m_packed_matrix
       *
       *  Options for m_matfrom are: 
       *
       *  - m_matform = eSymmetric_Positive implies Cholesky factorization 
       *  using lapack::dpptrf.
       *
       *  - m_matform = eSymmetric implies Factorisation using Bunch-Kaufman
       *  pivoting using lapack::dsptrf.
       *
       *  - m_matform = eSymmetric_Positive_Banded implies lower diagonal 
       *  banded cholesky factorisation using lapack::dpbtrf.
       *
       *  - m_matform = eGeneral_Banded implies factoring using 
       *  lapack::dgbtrf.
       *
       *  - m_matform = eGeneral_Full implies factoring using lapack::dgetrf.
       */
      void Factor();

      /** \brief evaluate \f$u \leftarrow m\_packed\_matrix^{-1} u\f$ for
       *  \a nrhs solves 
       *
       *  Options for mat_from are:
       *
       *  - m_matform = eSymmetric_Positive implies Cholesky back solve using
       *  lapack::dpptrs.
       *
       *  - m_matform = eSymmetric implies Back solve from using lapack::dsptrs.
       *
       *  - m_matform = eSymmetric_Positive_Banded implies lower diagonal 
       *  banded cholesky backsole using lapack::dpbtrs.
       *
       *  - m_matform = eGeneral_Banded implies LU back solve using 
       *  lapack::dgbtrs.
       *
       *  - m_matform = eGeneral_Full implies LU back solve using 
       *  lapack::dgetrs.
       */
      void Solve(double *u, const int nrhs);
      void SolveTranspose(double *u, const int nrhs);

      /** \brief Fill m_packed_matrix from m_matrix 
       *
       *  \b Note: Assume input matrix is stored in RowMajor 'C' type format.
       *  The packed matrix however is stored in Column Major Fortran type 
       *  format so that the appropriate Lapack routines can be called. This is
       *  only important for the General Matrix forms.
       */
      void FillPackedMatrix();

      // Spectral analysis routines

      /** \brief  evaluate the eigenvalues of the factored Matrix
       *
       *  Determine the eigen specturm of the Packed matrix by  returning its
       *  read values in \a er and its imaginary values in \a ei. 
       *
       *  if evecs is given then the eigenvectors are also returned. 
       *  This option is currently only setup for a General_Full matrix.
       */
      void   EigenValues(double *er, double *ei, double *evecs);

      /** \brief Write the eigenvalues of matrix m_packed_matrix to 
       *  file \a file
       */
      void   EigenValues(const char file[]);

      /** \brief return the maximum eigenvalue */
      double MaxEigenValue();

      /** \brief return the L2 condition number */
      double L2ConditionNo();

      /** \brief return the null space of the matrix
       *  
       *  Return the null space of the Packed Matrix 
       *  where the null space is defined by the number of eigenvalues
       *  with magnitude smaller than tol.
       */
      int    NullSpaceDim (double);


      // Private data access and setting

      /** \brief return private data for m_matform */
      MatrixForm GetMatForm()
      {
	return m_matform;
      }

      /** \brief set private data: m_matform */
      void SetMatForm(MatrixForm matform)
      {
	m_matform = matform;
      }

      /** \brief return private data for m_lda */
      inline int GetLda ()
      {
	return m_lda;
      }

      /** \brief set private data: m_lda */
      void SetLda    (const int val)
      {
	m_lda = val;
      }

      /** \brief return private data for m_bwidth */
      int GetBwidth    ()
      {
	return m_bwidth;
      }

      /** \brief set private data: m_bwidth */
      void SetBwidth (int val)
      {
	m_bwidth = val;
      }

      /** \brief return private data for m_ldiag */
      int GetLdiag    ()
      {
	return m_ldiag;
      }

      /** \brief set private data: m_ldiag */
      void SetLdiag (const int val)
      {
	m_ldiag = val;
      }

      /** \brief dump the matrix to file \a out */
      void DumpMatrix           (FILE *out);

      /** \brief show the matrix structure in file \a out */
      void ShowMatrixStructure (FILE *out);
      inline void ShowMatrixStructure (void)
      {
	  FILE *out = stdout; 
	  ShowMatrixStructure(out);
      }
      
      // Overloaded Operators
      friend bool operator  == (const StdExpansion &x, const StdMatContainer *y);
      friend bool operator  != (const StdExpansion &x, const StdMatContainer *y);

      friend bool operator  == (const StdMatContainer *x, const StdExpansion &y);
      friend bool operator  != (const StdMatContainer *x, const StdExpansion &y);

    private:
      int m_numbases;
      int * m_basisorder;     /**< Order per expansion basis */
      const Basis **m_base;   /**< Bases needed for the expansion */
      double * m_matrix;      /**< Matrix pointer */

      MatrixForm m_matform;    /**< enum list of type of matrix form */
      bool    m_factored;      /**< Flag to identify if matrix is factored */
      int    *m_ipiv;          /**< Pivoting array */
      int     m_lda;           /**< leading diagonal of matrix */
      /** Bandwdith for positive banded matrix Upper sug diagonals plus one
       *  for general banded matrix
       */
      int     m_bwidth;        
      /** Low sub diagonals for general banded matrix */
      int     m_ldiag;         
      double* m_packed_matrix; /**< Inverse/Factorisation of Matrix pointer */
    };

    class StdMatrix
    {

    public:

      StdMatrix()
      {
      };

      ~StdMatrix()
      {
	std::vector<StdMatContainer*>::iterator def;

	for(def = m_local_mass.begin(); def != m_local_mass.end(); ++def)
	{
	  delete def[0];
	}

	for(def = m_nbasis_trans.begin(); def != m_nbasis_trans.end(); ++def)
	{
	delete def[0];
	}

	for(def = m_local_lap.begin(); def != m_local_lap.end(); ++def)
	{
	  delete def[0];
	}
      }

      /** \brief Get the Local mass matrix \f$ \bf M\f$ given a standard
       *  expansion E
       */
      StdMatContainer * GetLocalMass  (StdExpansion * E);

      StdMatContainer * GetNBasisTrans(StdExpansion * E);

      /** \brief Get the Local  Laplacian matrix \f$ L\f$ given a standard
       *  expansion E
       */
      StdMatContainer * GetLocalLap   (StdExpansion * E);

    private:

      std::vector<StdMatContainer*> m_local_mass;
      std::vector<StdMatContainer*>::iterator m_local_mass_cur;

      std::vector<StdMatContainer*> m_nbasis_trans;
      std::vector<StdMatContainer*>::iterator m_nbasis_trans_cur;

      std::vector<StdMatContainer*> m_local_lap;
      std::vector<StdMatContainer*>::iterator m_local_lap_cur;

    };

  } // end of namespace
} // end of namespace


#endif //STDMATRIX_H

/**
 * $Log: StdMatrix.h,v $
 * Revision 1.4  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.3  2006/06/02 18:48:40  sherwin
 * Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 * Revision 1.2  2006/06/01 13:43:19  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.20  2006/03/12 21:59:48  sherwin
 *
 * compiling version of LocalRegions
 *
 * Revision 1.19  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.18  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.17  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/



