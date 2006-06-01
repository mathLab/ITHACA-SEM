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

      ///Forward Multiply
      void Mxv(const double *in, double *out);

      /// Matrix factorisation initialisation routine
      void SetInvInfo(StdExpansion *E, const MatrixType  Mform);

      //  Inversion routines

      /// declare matrix memory
      void SetMemPackedMatrix();

      /// factorise matrix depending upon definition stored in _mat_form and
      /// store in _packed_matrix
      void Factor();

      /// evaulate \f$u \leftarrow A^{-1} u\f$ for \a nrhs solves
      void Solve(double *u, const int nrhs);

      /// Fill _packed_matrix from _matrix
      void FillPackedMatrix();

      // Spectral analysis routines

      /// evaluate the eigenvlaues of factored Matrix
      void   EigenValues(double *er, double *ei, double *evecs);
      /// Write Eigenvalues to file
      void   EigenValues(const char file[]);
      /// Determine Maxximum eigenvalue
      double MaxEigenValue();
      /// L2 condition Number
      double L2ConditionNo();
      /// Null Space of the Matrix
      int    NullSpaceDim (double);


      // Private data access and setting

      /// return private data for _mat_form
      MatrixForm GetMatForm()
      {
	return m_matform;
      }

      /// set private data: _mat_form
      void SetMatForm(MatrixForm matform)
      {
	m_matform = matform;
      }

      /// return private data for m_lda
      inline int GetLda ()
      {
	return m_lda;
      }

      /// set private data: lda
      void SetLda    (const int val)
      {
	m_lda = val;
      }

      /// return private data for lda
      int GetBwidth    ()
      {
	return m_bwidth;
      }

      /// set private data: bwidth
      void SetBwidth (int val)
      {
	m_bwidth = val;
      }

      /// return private data for _ldiag
      int GetLdiag    ()
      {
	return m_ldiag;
      }

      /// set private data: bwidth
      void SetLdiag (const int val)
      {
	m_ldiag = val;
      }

      /// dump the matrix to file out
      void DumpMatrix           (FILE *out);

      /// show the matrix structure in file out
      void ShowMatrixStructure (FILE *out);

      // Overloaded Operators
      friend bool operator  == (const StdExpansion &x, const StdMatContainer *y);
      friend bool operator  != (const StdExpansion &x, const StdMatContainer *y);

      friend bool operator  == (const StdMatContainer *x, const StdExpansion &y);
      friend bool operator  != (const StdMatContainer *x, const StdExpansion &y);

    private:
      int m_numbases;
      int * m_basisorder;     ///< Order per expansion basis
      const Basis **m_base;   ///< Bases needed for the expansion
      double * m_matrix;      ///< Matrix pointer

      MatrixForm m_matform;    ///< enum list of type of matrix form
      bool    m_factored;      ///< Flag to identify if matrix is factored
      int    *m_ipiv;          ///< Pivoting array
      int     m_lda;           ///< leading diagonal of matrix
      int     m_bwidth;        ///< Bandwdith for positive banded matrix \n
                               ///< Upper sug diagonals plus one for general
                               ///< banded matrix
      int     m_ldiag;         ///< Low sub diagonals for general banded matrix
      double* m_packed_matrix; ///< Inverse/Factorisation of Matrix pointer
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

      StdMatContainer * GetLocalMass  (StdExpansion * E);
      StdMatContainer * GetNBasisTrans(StdExpansion * E);
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



