///////////////////////////////////////////////////////////////////////////////
//
// File: NistSparseBlas.hpp
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
// Description: wrapper of functions around NIST sparse BLAS routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSEBLAS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSEBLAS_HPP


namespace SparseBlas
{
    extern "C"
    {
        // -- BLAS Level 1:
        // -- BLAS level 2:

        // -- BLAS level 3:

        ///////////////////////////////////////////////////////////////////////////////////
        //    dcsrmm -- compressed sparse row format matrix-matrix multiply
        //
        //    C <- alpha A B + beta C
        //
        //    Arguments:
        //
        //    int transa   Indicates how to operate with the sparse matrix
        //   		   0 : operate with matrix
        //   		   1 : operate with transpose matrix
        //
        //    int m	   Number of - rows in matrix A
        //                           - rows in matrix C
        //
        //    int n	   Number of - columns in matrix B
        //                           - columns in matrix C
        //
        //    int k	   Number of - columns in matrix A
        //                           - rows in matrix B
        //   
        //    double alpha Scalar parameter
        //       
        //    double beta  Scalar parameter
        //      
        //    int descra[] Descriptor argument.  Nine element integer array
        //   		   descra[0] matrix structure
        //   			   0 : general
        //   			   1 : symmetric
        //   		           2 : Hermitian
        //   			   3 : Triangular
        //   			   4 : Skew(Anti)-Symmetric
        //   			   5 : Diagonal
        //   		   descra[1] upper/lower triangular indicator
        //   			   1 : lower
        //   			   2 : upper
        //   		   descra[2] main diagonal type
        //   			   0 : non-unit
        //   			   1 : unit
        //   		   descra[3] Array base 
        //   			   0 : C/C++ compatible
        //   			   1 : Fortran compatible
        //   		   descra[4] repeated indices?
        //   			   0 : unknown
        //   			   1 : no repeated indices
        //    
        //
        //    double *val  scalar array of length nnz containing matrix entries
        //
        //    int *indx    integer array of length nnz containing column indices
        //  
        //    int *pntrb   integer array of length k such that pntrb(j)-pntrb(1)
        //                 points to location in val of the first nonzero element in row j
        //
        //    int *pntre   integer array of length k such that pntre(j)-pntrb(1)
        //                 points to location in val of the last nonzero element in row j
        //
        //    double *b	   rectangular array with first dimension ldb
        //
        //    double *c	   rectangular array with first dimension ldc
        //
        //    double *work scratch array of length lwork.  lwork should be at least
        //   		   max(m,n)
        ///////////////////////////////////////////////////////////////////////////////////

        void  dcsrmm(const int transa, const int m, const int n, const int k,
                     const double alpha, const int descra[], const double val[],
                     const int indx[], const int pntrb[], const int pntre[],
                     const double b[], const int ldb,
                     const double beta, double c[], const int ldc,
                     double work[], const int lwork);


        ///////////////////////////////////////////////////////////////////////////////////
        //   dbsrmm -- block sparse row format matrix-matrix multiply
        //  
        //   C <- alpha A B + beta C
        //  
        //   Arguments:
        //  
        //   int transa	Indicates how to operate with the sparse matrix
        //  		0 : operate with matrix
        //  		1 : operate with transpose matrix
        //  
        //   int mb	Number of block rows in matrix A
        //  
        //   int n	Number of columns in matrix c
        //  
        //   int kb	Number of block columns in matrix A
        //  
        //   double alpha Scalar parameter
        //  
        //   double beta	Scalar parameter
        //  
        //   int descra[]	Descriptor argument.  Nine element integer array
        //  		descra[0] matrix structure
        //  			0 : general
        //  			1 : symmetric
        //  			2 : Hermitian
        //  			3 : Triangular
        //  			4 : Skew(Anti)-Symmetric
        //  			5 : Diagonal
        //  		descra[1] upper/lower triangular indicator
        //  			1 : lower
        //  			2 : upper
        //  		descra[2] main diagonal type
        //  			0 : non-unit
        //  			1 : unit
        //  		descra[3] Array base 
        //  			0 : C/C++ compatible
        //  			1 : Fortran compatible
        //  		descra[4] repeated indices?
        //  			0 : unknown
        //  			1 : no repeated indices
        //  
        //  
        //   double *val	scalar array of length nnz containing matrix entries
        //  
        //   int *bindx	integer array of length bnnz consisting of the block column
        //   		indices of the entries of A.
        //  
        //   int *bpntrb	integer array of length mb such that bpntrb(i)-bpntrb(1)
        //                points to location in bindx of the first block entry of 
        //  		the j-th row of A.
        //  
        //   int *bpntre	integer array of length mb such that bpntre(i)-bpntrb(1)
        //                points to location in bindx of the last block entry of
        //  		the j-th row of A.
        //  
        //   int lb	dimension of blocks
        //  
        //   double *b	rectangular array with first dimension ldb
        //  
        //   double *c	rectangular array with first dimension ldc
        //  
        //   double *work	scratch array of length lwork.  lwork should be at least
        //  		max(m,n)
        ///////////////////////////////////////////////////////////////////////////////////

        void  dbsrmm(const int transa, const int mb, const int n, const int kb,
                     const double alpha, const int descra[], const double val[],
                     const int bindx[], const int bpntrb[], const int bpntre[],
                     const int lb, const double b[], const int ldb,
                     const double beta, double c[], const int ldc,
                     double work[], const int lwork);


        void CSR_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);


        void CSRsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);


        void BCO_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base);

        void BCO_VecMult_CATB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base);

        void BSR_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base);

        void BSRsymm_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base);

        void BSR_VecMult_CATB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base);


    } // extern "C"

    static inline void Dbsrmm(const int transa, const int mb, const int n, const int kb,
                              const double alpha, const int descra[], const double val[],
                              const int bindx[], const int bpntrb[], const int bpntre[],
                              const int lb, const double b[], const int ldb,
                              const double beta, double c[], const int ldc,
                              double work[], const int lwork) 
    {
        dbsrmm (transa,mb,n,kb,alpha,descra,val,bindx,
                         bpntrb,bpntre,lb,b,ldb,beta,c,ldc,work,lwork);
    }


    static inline void Dcsrmm(const int transa, const int m, const int n, const int k,
                              const double alpha, const int descra[], const double val[],
                              const int indx[], const int pntrb[], const int pntre[],
                              const double b[], const int ldb,
                              const double beta, double c[], const int ldc,
                              double work[], const int lwork) 
    {
        dcsrmm (transa,m,n,k,alpha,descra,val,indx,
                         pntrb,pntre,b,ldb,beta,c,ldc,work,lwork);
    }

    static inline void Dcsrmv(const int nrows, const int ncols,
                              const double *mat, const int *indx, const int *ptr,
                              double *b, double*c, double* work)
    {
        static int descra[] = {0,0,0,0,0,0,0,0,0};

        dcsrmm (0,nrows,1,ncols,1.0,descra,mat,indx,ptr,ptr+1,b,nrows,0.0,c,nrows,work,nrows);
    }
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSEBLAS_HPP

