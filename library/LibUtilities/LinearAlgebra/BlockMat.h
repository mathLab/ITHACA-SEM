// #ifndef NEKTAR_LIB_UTILITIES_BLOCK_MAT_H
// #define NEKTAR_LIB_UTILITIES_BLOCK_MAT_H
//
// #include <LibUtilities/Vmath.hpp>
// #include <LibUtilities/Lapack.hpp>
// #include <LibUtilities/ErrorUtil.hpp>
//
// #include <stdlib.h>
// #include <iostream>
// #include <vector>
//
// #include <assert.h>
//
// namespace blockmat
// {
//
//     enum MatStorage
//     {
//         RowMajor,  ///< Row major matrix storage
//         ColMajor   ///< column major matrix storage
//     };
//
//     class BlockSubMat
//     {
//         public:
//             /// default constructor  - no matrix
//             BlockSubMat(const int id, const int rows, const int cols)
//             {
//                 _id   = id;
//                 _rows = rows;
//                 _cols = cols;
//                 _mat  = new double [_rows*_cols];
//                 Vmath::Zero(_rows*_cols,_mat,1);
//             }
//
//     /// default constructor - matrix
//     BlockSubMat(const int id, const int rows, const int cols,
//           const double *mat)
//     {
//         _id   = id;
//         _rows = rows;
//         _cols = cols;
//         _mat  = new double [_rows*_cols];
//         Vmath::Vcopy(_rows*_cols,mat,1,_mat,1);
//     }
//
//     // copy constructor
//     BlockSubMat(const BlockSubMat& a)
//     {
//         _id   = a._id;
//         _rows = a._rows;
//         _cols = a._cols;
//         _mat  = new double [_rows*_cols];
//         Vmath::Vcopy(_rows*_cols,a._mat,1,_mat,1);
//     }
//
//     /// default destructor
//     ~BlockSubMat()
//     {
//         if(_rows*_cols)
//         {
//             if(_mat) delete[] _mat;
//         }
//     }
//
//     void AddVal(const int id1, const int id2, const double val)
//     {
//         _mat[id1*_cols+id2] += val;
//     }
//
//     void AddMat(const int rows, const int cols, const double *mat)
//     {
//
//         if((_rows != rows) || (_cols != cols))
//         {
//             ErrorUtil::Error(ErrorUtil::ewarning,"BlVeck::AddMat",
//                 "Adding sub-matrices of different dimension");
//         }
//
//         Vmath::Vadd(_rows*_cols,mat,1,_mat,1,_mat,1);
//     }
//
//
//     /// y = alpha*A*v + beta y
//     /// formA,  provides the storage format of matrix A
//     void geMxv(MatStorage formA, const double& alpha, const double *v,
//            const double& beta, double *y)
//     {
//
//         switch(formA)
//         {
//             case RowMajor:
//                 Blas::Dgemv('T',_cols,_rows,alpha,_mat,_cols,v,1,beta,y,1);
//                 break;
//
//             case ColMajor:
//                 Blas::Dgemv('N',_cols,_rows,alpha,_mat,_cols,v,1,beta,y,1);
//                 break;
//
//             default:
//                 ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxv",
//                     "unknown storage format for matrix A");
//
//             break;
//         }
//     }
//
//     void Mxvpy(const double * v, double * y)
//     {
//         Blas::Dgemv('T',_cols,_rows,1.0,_mat,_cols,v,1,1.0,y,1);
//     }
//
//     void Mtxvpy(const double * v, double * y)
//     {
//         Blas::Dgemv('N',_cols,_rows,1.0,_mat,_cols,v,1,1.0,y,1);
//     }
//
//
//     void scal(const double alpha)
//     {
//         Blas::Dscal(_rows*_cols,alpha,_mat,1);
//     }
//
//     void axpy(const double alpha, const BlockSubMat& A)
//     {
//         Blas::Daxpy(A._rows*A._cols,alpha,A._mat,1,_mat,1);
//     }
//
//     /// alpha*A*B + beta (this)
//     /// formA, formB provide the storage format of matrices A and B
//     void geMxM(MatStorage formA, MatStorage formB, const double& alpha,
//            const BlockSubMat& A, const BlockSubMat& B, const double& beta)
//     {
//
//         switch(formB)
//         {
//             case RowMajor:
//                 switch(formA)
//                 {
//                     case RowMajor:
//                         Blas::Cdgemm(A._rows,B._cols,A._cols,alpha,A._mat,
//                             A._cols,B._mat,B._cols,beta,_mat,_cols);
//                     break;
//
//                     case ColMajor:
//                         Blas::Dgemm('N','T',B._cols,A._cols,A._rows,alpha,B._mat,
//                                 B._cols,A._mat,A._cols,beta,_mat,_cols);
//                         break;
//                     default:
//                         ErrorUtil::Error(ErrorUtil::ewarning, "BlockMat::geMxM",
//                                 "unknown storage format for matrix A");
//                         break;
//                 }
//                 break;
//
//
//             case ColMajor:
//                 switch(formA)
//                 {
//                     case RowMajor:
//                         Blas::Dgemm('T','N',B._rows,A._rows,A._cols,alpha,B._mat,
//                                 B._cols,A._mat,A._cols,beta,_mat,_cols);
//                         break;
//                     case ColMajor:
//                         Blas::Dgemm('N','N',A._rows,B._cols,A._cols,alpha,A._mat,
//                                 A._cols,B._mat,B._cols,beta,_mat,_cols);
//                         break;
//                     default:
//                         ErrorUtil::Error(ErrorUtil::ewarning, "BlockMat::geMxM",
//                                 "unknown storage format for matrix A");
//                         break;
//                 }
//                 break;
//
//             default:
//                 ErrorUtil::Error(ErrorUtil::ewarning, "BlockMat::geMxM",
//                     "unknown storage format for matrix B");
//             break;
//         }
//
//     }
//
//     // A * B + (this)
//     void MxMpM(const BlockSubMat& a, const BlockSubMat& b)
//     {
//         Blas::Cdgemm(a._rows,b._cols,a._cols,1.0,a._mat,
//             a._cols,b._mat,b._cols,1.0,_mat,_cols);
//     }
//
//     // A^T * B + (this)
//     void MtxMpM(const BlockSubMat& a, const BlockSubMat& b)
//     {
//         Blas::Dgemm('N','T',b._cols,a._cols,a._rows,1.0,b._mat,
//             b._cols,a._mat,a._cols,1.0,_mat,_cols);
//     }
//
//     // A * B^T + (this)
//     void MxMtpM(const BlockSubMat& a, const BlockSubMat& b)
//     {
//         Blas::Dgemm('T','N',b._rows,a._rows,a._cols,1.0,b._mat,
//             b._cols,a._mat,a._cols,1.0,_mat,_cols);
//     }
//
//     void PrintBlock()
//     {
//         for(int i = 0; i < _rows; ++i)
//         {
//             for(int j = 0; j < _cols; ++j)
//             {
//                 std::cout << _mat[i*_cols+j] << " ";
//             }
//             std::cout << std::endl;
//         }
//     }
//
//     int     get_rows() const { return _rows; }
//     int     get_cols() const { return _cols; }
//     int     get_id()   const { return _id;   }
//     double *get_mat() { return _mat; }
//     double  get_val(const int i, const int j) { return _mat[i*_cols+j]; }
//
//     // copy operator;
//     BlockSubMat& operator = (const BlockSubMat& a)
//     {
//         if(this != &a)
//         {
//             Vmath::Vcopy(_rows*_cols,a._mat,1,_mat,1);
//         }
//       return *this;
//     }
//     // (this) = (this) + a
//     BlockSubMat& operator += (const BlockSubMat& a)
//     {
//         Vmath::Vadd(_rows*_cols,a._mat,1,_mat,1,_mat,1);
//         return *this;
//     }
//     // (this) = (this) - a
//     BlockSubMat& operator -= (const BlockSubMat& a)
//     {
//         Vmath::Vsub(_rows*_cols,_mat,1,a._mat,1,_mat,1);
//         return *this;
//     }
//
//     void neg()
//     {
//         Vmath::Neg(_rows*_cols,_mat,1);
//     }
//
//     void zero()
//     {
//         Vmath::Zero(_rows*_cols,_mat,1);
//     }
//
//     void invert(void)
//     {
//         int    *ipiv, info;
//         double *Wk;
//
//         if(_rows != _cols)
//         {
//             ErrorUtil::Error(ErrorUtil::efatal,"BlockSubMat::invert()","Matrix not square");
//         }
//
//         ipiv = new int   [_rows];
//         Wk   = new double[_rows*_cols];
//
//         // invert C as a general matrix
//         Lapack::Dgetrf(_rows,_cols,_mat,_cols,ipiv,info);          assert(info==0);
//         Lapack::Dgetri(_rows,_mat,_cols,ipiv,Wk,_rows*_cols,info); assert(info==0);
//
//         delete[] Wk;
//         delete[] ipiv;
//     }
//
//     void rescale(const double *scale)
//     {
//         if(_rows != _cols)
//         {
//             ErrorUtil::Error(ErrorUtil::efatal,"BlockSubMat::rescale",
//                 "_rows not equal to _cols");
//         }
//
//         for(int i = 0; i < _rows; ++i)
//         {
//             // scale rows
//             Vmath::Vmul(_cols,scale,1,_mat+i*_cols,1,_mat+i*_cols,1);
//             // scale cols
//             Vmath::Vmul(_rows,scale,1,_mat+i,_cols,_mat+i,_cols);
//         }
//     }
//
//     private:
//         int _id;       // id of this submatrix
//         int _rows;     // number of rows in matrix;
//         int _cols;     // number of columns in matrix;
//         double *_mat;  // linear matrix stored in row major format
//
//     };
//
//     // local Class for differential matrix structure
//     class BlockVec
//     {
//         private:
//             friend class BlockMat;
//             friend class BlockSubMat;
//
//         public:
//             /// default constructor
//             BlockVec()
//             {
//                 _max_blk = 0;
//             }
//
//             /// default destructor
//             ~BlockVec()
//             {
//                 for(int i = 0; i < _Block.size(); ++i)
//                 {
//                     delete _Block[i];
//                 }
//                 delete[] _entries;
//             }
//
//             /// default destructor
//             void reset()
//             {
//                 for(int i = 0; i < _Block.size(); ++i)
//                 {
//                     delete _Block[i];
//                 }
//
//                 _Block.erase(_Block.begin(),_Block.end());
//                 Vmath::Fill(_max_blk,-1,_entries,1);
//             }
//
//             void SetMem(int nblk)
//             {
//                 _max_blk = nblk;
//                 _entries = new int[_max_blk];
//                 Vmath::Fill(_max_blk,-1,_entries,1);
//             }
//
//             void GenBlock(const int id, const int rows, const int cols,
//                 const double *mat);
//
//             void GenBlock(const int id, const int rows, const int cols,
//                 const int id1, const int id2, const double val);
//
//             void PrintBlocks(const int *offset)
//             {
//                 int i,j,k,id,cols;
//                 double *mat;
//
//                 if(_Block.size())
//                 {
//                     for(i = 0; i < _Block[0]->get_rows(); ++i)
//                     {
//                         for(j = 0; j < _max_blk; ++j)
//                         {
//                             cols = offset[j+1]-offset[j];
//                             if((id=_entries[j])+1)
//                             {
//                                 mat = _Block[id]->get_mat();
//                                 for(k = 0; k < cols; ++k)
//                                 {
//                                     std::cout << mat[i*cols+k] << " ";
//                                 }
//                             }
//                             else
//                             {
//                                 for(k = 0; k < cols; ++k)
//                                 {
//                                     std::cout <<  "0 ";
//                                 }
//                             }
//                         }
//                         std::cout << std::endl;
//                     }
//                 }
//                 else
//                 {
//                     std::cout << "-- Empty row --" << std::endl;
//                 }
//             }
//
//             BlockVec& add    (const BlockVec& a,    const BlockVec& b);
//             BlockVec& sub    (const BlockVec& a,    const BlockVec& b);
//             BlockVec& MxMpM  (const BlockSubMat& a, const BlockVec& b);
//             BlockVec& MtxMpM (const BlockSubMat& a, const BlockVec& b);
//             BlockVec& geMxM  (MatStorage formA, MatStorage formB,
//                     const double& alpha, const BlockSubMat& A,
//                     const BlockVec& B, const double& beta);
//             // (this) = (this) + alpha * A
//             BlockVec& axpy  (const double alpha, const BlockVec& A);
//
//             double* Mxvpy  (const int *offset, const double *v, double *y);
//             double* Mtxvpy (const int *offset, const double *v, double *y);
//             double* geMxv  (MatStorage form, const double& alpha, const double *v,
//                     const double& beta, double *y, const int *offset);
//
//
//
//         private:
//             int _max_blk;
//             int *_entries;
//
//             std::vector <BlockSubMat*> _Block; // Block entry in Row
//     };
//
//     // local Class Block Matrix
//     class BlockMat
//     {
//         public:
//             /// default constructor
//             BlockMat(int rowblk, int colblk)
//             {
//                 _max_rowblk = rowblk;
//                 _max_colblk = colblk;
//
//                 _offset = (int *)NULL;
//
//                 _Row = new BlockVec [_max_rowblk];
//                 for(int i =0; i < _max_rowblk; ++i)
//                 {
//                     _Row[i].SetMem(_max_colblk);
//                 }
//             }
//
//             /// default destructor
//             ~BlockMat()
//             {
//                 delete[] _Row;
//
//                 if(_offset)
//                 {
//                     delete[] _offset;
//                 }
//             }
//
//             void Reset_Rows()
//             {
//                 int i,j;
//
//                 for(i=0; i < _max_rowblk; ++i)
//                 {
//                     for(j=0;j < _Row[i]._Block.size(); ++j)
//                     {
//                         _Row[i].reset();
//                     }
//                 }
//             }
//
//             void GenBlock(const int row_id, const int col_id, const int rows,
//                 const int cols, const double *mat)
//             {
//                 _Row[row_id].GenBlock(col_id,rows,cols,mat);
//             }
//
//             void GenBlock(const int row_id, const int col_id, const int rows,
//                 const int cols, const int id1, const int id2,
//                 const double val)
//             {
//                 _Row[row_id].GenBlock(col_id,rows,cols,id1,id2,val);
//             }
//
//             int get_tot_row_entries()
//             {
//                 if(!_offset)
//                 {
//                     setup_offset();
//                 }
//                 return _offset[_max_rowblk];
//             }
//
//             void setup_offset()
//             {
//                 // set up _offset if not already done at init
//
//                 if(!_offset) _offset = new int [_max_colblk+1];
//                 Vmath::Zero(_max_colblk+1,_offset,1);
//
//                 for(int i = 0; i < _max_rowblk; ++i)
//                 for(int j = 0; j < _Row[i]._Block.size(); ++j)
//                 update_offset(_Row[i]._Block[j]->get_id(),_Row[i]._Block[j]->get_cols());
//             }
//
//             void update_offset(int id,int cols)
//             {
//                 int i;
//                 // check and update _offset - assumes
//                 if(_offset[id+1] == _offset[id])
//                 {
//                     // assume cols not set
//                     for(i = id+1; i <= _max_colblk; ++i)
//                     {
//                         _offset[i] += cols;
//                     }
//                 }
//             }
//
//             BlockMat& invert_diag();
//
//             void PrintBlocks()
//             {
//
//                 if(!_offset) setup_offset();
//
//                 for(int i = 0; i < _max_rowblk; ++i)
//                 _Row[i].PrintBlocks(_offset);
//             }
//
//             BlockMat& add    (const BlockMat& a, const BlockMat& b);
//             BlockMat& sub    (const BlockMat& a, const BlockMat& b);
//             BlockMat& geMxM  (MatStorage formA, MatStorage formB, const double& alpha,
//                     const BlockMat& A, const BlockMat& b,   const double& beta);
//             BlockMat& MxM    (const BlockMat& a, const BlockMat& b);
//             BlockMat& MtxM   (const BlockMat& a, const BlockMat& b);
//             BlockMat& MxMt   (const BlockMat& a, const BlockMat& b);
//
//             //(this) = A*(this) where A is block diagonal
//             BlockMat& diagMxy (const BlockMat& A);
//
//             // (this) = (this) + alpha * A
//             BlockMat& axpy  (const double alpha, const BlockMat& A);
//
//             double* geMxv  (MatStorage form, const double& alpha, const double *v,
//                     const double& beta,  double *y);
//             double* Mxvpy  (double* v, double* y);
//             double* Mtxvpy (double* v, double* y);
//
//             int  get_mrowblk(){ return _max_rowblk;}
//             int  get_mcolblk(){ return _max_colblk; }
//             int* get_offset(){ if(!_offset) setup_offset(); return _offset;}
//
//             int get_rows(const int i, const int j)
//             {
//                 int id;
//                 if((id = _Row[i]._entries[j])+1)
//                 {
//                     return _Row[i]._Block[id]->get_rows();
//                 }
//                 else
//                 {
//                     return -1;
//                 }
//             }
//
//             int get_cols(const int i, const int j)
//             {
//                 int id;
//                 if((id = _Row[i]._entries[j])+1)
//                 {
//                     return _Row[i]._Block[id]->get_cols();
//                 }
//                 else
//                 {
//                     return -1;
//                 }
//             }
//
//             double* get_mat(const int i, const int j, int& n, int &m)
//             {
//                 int id;
//                 if((id = _Row[i]._entries[j])+1)
//                 {
//                     n = _Row[i]._Block[id]->get_rows();
//                     m = _Row[i]._Block[id]->get_cols();
//                     return _Row[i]._Block[id]->get_mat();
//                 }
//                 else
//                 {
//                     return (double*)NULL;
//                 }
//             }
//
//             double get_val(const int i, const int j, const int id1, const int id2)
//             {
//                 int id;
//                 if((id = _Row[i]._entries[j])+1)
//                 {
//                     int m = _Row[i]._Block[id]->get_cols();
//                     return _Row[i]._Block[id]->get_val(id1,id2);
//                 }
//                 else
//                 {
//                     return (double) NULL;
//                 }
//             }
//
//             void rescale(const int i, const int j, double *scale)
//             {
//                 int id;
//                 if((id = _Row[i]._entries[j])+1)
//                 {
//                     _Row[i]._Block[id]->rescale(scale);
//                 }
//             }
//
//             void neg()
//             {
//                 int i,j;
//
//                 for(i=0; i < _max_rowblk; ++i)
//                 {
//                     for(j=0;j < _Row[i]._Block.size(); ++j)
//                     {
//                         _Row[i]._Block[j]->neg();
//                     }
//                 }
//             }
//
//             int cnt_blks()
//             {
//                 int i,cnt;
//                 for(cnt = i =0; i < _max_rowblk; ++i)
//                 {
//                     cnt += _Row[i]._Block.size();
//                 }
//
//                 return cnt;
//             }
//
//             double AmAt();
//
//
//         private:
//             int _max_colblk;
//             int _max_rowblk;
//             int *_offset;  ///< offset of dof entries in row
//
//             BlockVec *_Row; // Block entry in Row
//     };
//
// } // end of blockmat namespace
//
//
// #endif // NEKTAR_LIB_UTILITIES_BLOCK_MAT_H

/**
 * $Log: BlockMat.h,v $
 * Revision 1.2  2006/05/14 21:33:25  bnelson
 * Updates to fix compile errors and format for the coding convention.
 *
 * Revision 1.1  2006/05/04 18:57:41  kirby
 * *** empty log message ***
 *
 * Revision 1.3  2006/02/26 21:13:45  bnelson
 * Fixed a variety of compiler errors caused by updates to the coding standard.
 *
 * Revision 1.2  2006/01/31 13:51:12  bnelson
 * Updated for new configure.
 *
 **/

