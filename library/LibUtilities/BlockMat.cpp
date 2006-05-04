#include <LibUtilities/BlockMat.h>

namespace blockmat {

  void BlockVec::GenBlock(const int id, const int rows, const int cols,
              const double *mat){
    int en;
    BlockSubMat *B;

    if((en = _entries[id])+1){// entry already exists so add matrix
      _Block[en]->AddMat(rows,cols,mat);
    }
    else{ // add new entry
      B = new BlockSubMat(id,rows,cols,mat);

      _entries[id] = _Block.size();

      _Block.push_back(B);
    }
  }

  void BlockVec::GenBlock(const int id, const int rows, const int cols,
              const int id1, const int id2, const double val){
    int en;
    BlockSubMat *B;

    if((en = _entries[id])+1){// entry already exists so add matrix
    if((_Block[en]->get_rows() != rows)||(_Block[en]->get_cols() != cols))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockVec::GenBlock","add value to Block which "
          "does not have the same dimensions as calling argument");
      _Block[en]->AddVal(id1,id2,val);
    }
    else{ // add new entry
      B = new BlockSubMat(id,rows,cols);
      B->AddVal(id1,id2,val);

      _entries[id] = _Block.size();

      _Block.push_back(B);
    }
  }

  // (this) = a + b
  BlockVec& BlockVec::add(const BlockVec& a, const BlockVec& b){
    int i,id,en1,en,cnt=0;
    BlockSubMat *B;

    if(this == &b){ // need to to b-matrix first

      // check over 'a' row entries first
      for(i = 0; i < a._Block.size(); ++i){
    if((en = _entries[id = a._Block[i]->get_id()])+1){
      *(_Block[en]) += *(a._Block[i]);
    }
    else{ // need to declare memory and copy  in a matrix
      B = new BlockSubMat(*a._Block[i]);
      en = _entries[id] = _Block.size();
      _Block.push_back(B);
    }
      }
    }
    else{

      // check over 'a' row entries first
      for(i = 0; i < a._Block.size(); ++i){
    if((en = _entries[id = a._Block[i]->get_id()])+1)
      *(_Block[en]) = *(a._Block[i]);
    else{ // need to declare memory and copy  in a matrix
      B = new BlockSubMat(*a._Block[i]);
      en = _entries[id] = _Block.size();
      _Block.push_back(B);
    }

    //check to see if b entry exists and if so add
    if((en1 = b._entries[id])+1)
      *(_Block[en]) += *(b._Block[en1]);
      }

      cnt = a._Block.size();

      // check over 'b' vector  entries for any independent entries not in 'a'
      for(i = 0; i < b._Block.size(); ++i){
    if((a._entries[id = b._Block[i]->get_id()]) == -1){
      if((en = _entries[id])+1)
        *(_Block[en]) = *(b._Block[i]);
      else{
        B = new BlockSubMat(*b._Block[i]);
        _entries[id] = _Block.size();
        _Block.push_back(B);
      }
      ++cnt;
    }
      }


     //If original row is larger than new contributions then delete extra terms
      if(_Block.size() != cnt){
    for(i = 0; i < _Block.size(); ++i){
      id = _Block[i]->get_id();
      if((a._entries[id] == -1)&&(b._entries[id] == -1)){
        // remove element
        B = _Block[i] ;
        _Block.erase(_Block.begin()+i);
        _entries[i] = -1;
        delete B;
      }
    }
    if(_Block.size() != cnt)
      ErrorUtil::Error(ErrorUtil::efatal,"BlockVec::add","incorrect final size");
      }
    }

    return *this;
  }

  // (this) = a - b
  BlockVec& BlockVec::sub(const BlockVec& a, const BlockVec& b){
    int i,id,en1,en,cnt;
    BlockSubMat *B;

    if(this == &b){ // need to to b-matrix first

      // Negate b entries
      for(i = 0; i < b._Block.size(); ++i)
    _Block[i]->neg();


      // add 'a'  entries
      for(i = 0; i < a._Block.size(); ++i){
    if((en = _entries[id = a._Block[i]->get_id()])+1)
      *(_Block[en]) += *(a._Block[i]);
    else{ // need to declare memory and copy  in a matrix
      B = new BlockSubMat(*a._Block[i]);
      en = _entries[id] = _Block.size();
      _Block.push_back(B);
    }
      }

    }
    else{
      // check over 'a' row entries first
      for(i = 0; i < a._Block.size(); ++i){
    if((en = _entries[id = a._Block[i]->get_id()])+1)
      *(_Block[en]) = *(a._Block[i]);
    else{ // need to declare memory and copy  in a matrix
      B = new BlockSubMat(*a._Block[i]);
      en = _entries[id] = _Block.size();
      _Block.push_back(B);
    }
    //check to see if b entry exists and if so add
    if((en1 = b._entries[id])+1)
      *(_Block[en]) -= *(b._Block[en1]);
      }

      cnt = a._Block.size();

      // check over 'b' row entries for any independent entries
      for(i = 0; i < b._Block.size(); ++i){
    if((a._entries[id = b._Block[i]->get_id()]) == -1){
      if((en = _entries[id])+1){
        *(_Block[en]) = *(b._Block[i]);
        _Block[en]->neg();
      }
    else{
      B = new BlockSubMat(*b._Block[i]);
      B->neg();
      _entries[id] = _Block.size();
      _Block.push_back(B);
    }
      ++cnt;
    }
      }

      // If original row entries are  larger than new contributions
      // then delete extra terms
      if(_Block.size() != cnt){
    for(i = 0; i < _Block.size(); ++i){
      id = _Block[i]->get_id();
      if((a._entries[id] == -1)&&(b._entries[id] == -1)){
        // remove element
        B = _Block[i] ;
        _Block.erase(_Block.begin()+i);
        _entries[i] = -1;
        delete B;
      }
    }
    if(_Block.size() != cnt)
      ErrorUtil::Error(ErrorUtil::efatal,"BlockVec::sub","incorrect final size");
      }
    }

    return *this;
  }


  // (this) += a*v where a is a submatrix
  BlockVec& BlockVec::axpy(const double alpha, const BlockVec& A){
    int en,i,id;

    // check over 'A' row entries
    for(i = 0; i < A._Block.size(); ++i){
      if((en = _entries[id = A._Block[i]->get_id()])+1)
    _Block[en]->axpy(alpha,*(A._Block[i]));
      else{ // need to declare memory and copy  in a matrix
    BlockSubMat *B;
    B = new BlockSubMat(*A._Block[i]);
    B->scal(alpha);
    en = _entries[id] = _Block.size();
    _Block.push_back(B);
      }
    }

    return *this;
  }



  // y += a*v where a is a submatrix
  double* BlockVec::geMxv(MatStorage form, const double& alpha, const double *v,
            const double& beta, double *y, const int *offset){
    int i,id;

    switch(form){
    case RowMajor:

      // Multiply a by v and add to  y
      for(i = 0; i < _Block.size(); ++i){
    id = _Block[i]->get_id();
    _Block[i]->geMxv(form,alpha,v+offset[id],beta,y);
      }
      break;
    case ColMajor:

      // Multiply a by v and add to  y
      for(i = 0; i < _Block.size(); ++i){
    id = _Block[i]->get_id();
    _Block[i]->geMxv(form,alpha,v,beta,y+offset[id]);
      }
      break;
    default:
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxv",
              "matrix format not specified for matrix A");
      break;
    }

    return y;
  }


  // y += a*v where a is a submatrix
  double* BlockVec::Mxvpy(const int *offset, const double *v, double *y){
    int i,id;

    // Multiply a by v and add to  y
    for(i = 0; i < _Block.size(); ++i){
      id = _Block[i]->get_id();
      _Block[i]->Mxvpy(v+offset[id],y);
    }

    return y;
  }

  // y += a^t*v ywhere a is a submatrix
  double* BlockVec::Mtxvpy(const int *offset, const double *v, double *y){
    int i,id;

    // Multiply a by v and add to  y
    for(i = 0; i < _Block.size(); ++i){
      id = _Block[i]->get_id();
      _Block[i]->Mtxvpy(v,y+offset[id]);
    }

    return y;
  }


  // (this) += a*b where a is a submatrix
  BlockVec& BlockVec::MxMpM(const BlockSubMat& a, const BlockVec& b){
    int i,id,en;

    // Multiply a by all b entries and put into (*this)
    for(i = 0; i < b._Block.size(); ++i){
      if(b._Block[i]->get_rows() != a.get_cols())
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockVec::MxMpM","cols and rows do not match");

      if((en = _entries[id = b._Block[i]->get_id()])+1)
    _Block[en]->MxMpM(a,*(b._Block[i]));
      else{ // need to declare memory and copy  in a matrix
    BlockSubMat *B;
    B = new BlockSubMat(id,a.get_rows(),b._Block[i]->get_cols());
    B->MxMpM(a,*(b._Block[i]));
    en = _entries[id] = _Block.size();
    _Block.push_back(B);
      }
    }

    return *this;
  }


  // (this) += a^t*b where a is a submatrix
  BlockVec& BlockVec::MtxMpM(const BlockSubMat& a, const BlockVec& b){
    int i,id,en;

    // Multiply a by all b entries and put into (*this)
    for(i = 0; i < b._Block.size(); ++i){
      if(b._Block[i]->get_rows() != a.get_rows())
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockVec::axbpy","cols and rows do not match");

      if((en = _entries[id = b._Block[i]->get_id()])+1)
    _Block[en]->MtxMpM(a,*(b._Block[i]));
      else{ // need to declare memory and copy  in a matrix
    BlockSubMat *B;
    B = new BlockSubMat(id,a.get_cols(),b._Block[i]->get_cols());
    B->MtxMpM(a,*(b._Block[i]));
    en = _entries[id] = _Block.size();
    _Block.push_back(B);
      }
    }

    return *this;
  }


  //----------------------------------------------------------------

  // (this) = a + b
  BlockMat& BlockMat::add(const BlockMat& a, const BlockMat& b){

    if((a._max_rowblk != _max_rowblk)||(b._max_rowblk != _max_rowblk)||
       (a._max_colblk != _max_colblk)||(b._max_colblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::sub","Matrix rows/cols not compatible");

    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].add(a._Row[i],b._Row[i]);

    return *this;
  }

  // (this) = a - b
  BlockMat& BlockMat::sub(const BlockMat& a, const BlockMat& b){

    if((a._max_rowblk != _max_rowblk)||(b._max_rowblk != _max_rowblk)||
       (a._max_colblk != _max_colblk)||(b._max_colblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::sub","Matrix rows/cols not compatible");

    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].sub(a._Row[i],b._Row[i]);

    return *this;
  }


  // (this) = (this) + alpha*A
  BlockMat& BlockMat::axpy(const double alpha, const BlockMat& A){

    if((A._max_rowblk != _max_rowblk)||(A._max_colblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::axMpM","Matrix rows/cols not compatible");

    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].axpy(alpha,A._Row[i]);

    return *this;
  }


  // y = alpha A*v + beta y
  double* BlockMat::geMxv(MatStorage form, const double& alpha, const double *v,
            const double& beta, double *y){
    int i,n;

    if(!_offset) setup_offset();

    switch(form){
    case RowMajor:
      // multiply every row by v
      for(n=i=0; i < _max_rowblk; ++i){
    _Row[i].geMxv(form,alpha,v,beta,y+n,_offset);
    n += _Row[i]._Block[0]->get_rows();
      }
      break;
    case ColMajor:
      // multiply every row by v
      for(n=i=0; i < _max_rowblk; ++i){
    _Row[i].geMxv(form,alpha,v+n,beta,y,_offset);
    n += _Row[i]._Block[0]->get_rows();
      }
      break;
    default:
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxv",
            "matrix format not specified for matrix A");
      break;
    }

    return y;
  }

  // y  = A * v + y
  double* BlockMat::Mxvpy(double* v, double* y){
    int n,i;

    if(!_offset) setup_offset();

    // multiply every row by v
    for(n=i=0; i < _max_rowblk; ++i){
      _Row[i].Mxvpy(_offset,v,y+n);
      n += _Row[i]._Block[0]->get_rows();
    }

    return y;
  }

  // y  = A^t * v + y
  double* BlockMat::Mtxvpy(double* v, double* y){
    int n,i;

    if(!_offset) setup_offset();

    // multiply every row by v
    for(n=i=0; i < _max_rowblk; ++i){
      _Row[i].Mtxvpy(_offset,v+n,y);
      n += _Row[i]._Block[0]->get_rows();
    }
    return y;
  }

  // (this) = alpha*A*B + beta (this)
  BlockVec& BlockVec::geMxM(MatStorage formA, MatStorage formB,
            const double& alpha, const BlockSubMat& A,
            const BlockVec& B, const double& beta){
    int i,id,en;

    // Multiply A by all b entries and put into (*this)
    for(i = 0; i < B._Block.size(); ++i){
      if((en = _entries[id = B._Block[i]->get_id()])+1)
    _Block[en]->geMxM(formA,formB,alpha,A,*(B._Block[i]),beta);
      else{ // need to declare memory and copy  in a matrix
    BlockSubMat *Btmp;
    Btmp = new BlockSubMat(id,A.get_rows(),B._Block[i]->get_cols());
    Btmp->geMxM(formA,formB,alpha,A,*(B._Block[i]),beta);
    en = _entries[id] = _Block.size();
    _Block.push_back(Btmp);
      }
    }
    return *this;
  }

  // (this) =  A *(this) where A is block diagonal
  BlockMat& BlockMat::diagMxy(const BlockMat& A){
    int i,j,rows,cols;
    BlockSubMat *Btmp = (BlockSubMat *)NULL;

    rows = cols = -1;

    // multiply every entry in (this) by matrix diagonal of A
    for(i=0; i < _max_rowblk; ++i){
      // check rows are the same
      if(A._Row[i]._Block[0]->get_rows() !=
     _Row[i]._Block[0]->get_rows())
    ErrorUtil::Error(ErrorUtil::efatal,"Blockmat::diagMxy","Rows are not the same");

      for(j=0;j < _Row[i]._Block.size(); ++j){
    // set up temporary matrix
    if((rows != _Row[i]._Block[j]->get_rows())
       ||(cols != _Row[i]._Block[j]->get_cols())){
      if(Btmp) delete Btmp;
      rows = _Row[i]._Block[j]->get_rows();
      cols = _Row[i]._Block[j]->get_cols();
      Btmp = new BlockSubMat(j,rows,cols);
    }

    // do multiplication
    Btmp->geMxM(RowMajor,RowMajor,1.0,
            A._Row[i]._Block[A._Row[i]._entries[i]][0],
            _Row[i]._Block[j][0],0.0);

    // copy into original storage;
    *(_Row[i]._Block[j]) = *Btmp;
      }
    }
    if(Btmp) delete Btmp;

    return *this;
  }


  // (this) = alpha* A * B + beta (this)
  BlockMat& BlockMat::geMxM(MatStorage formA, MatStorage formB,const double& alpha,
            const BlockMat& A, const BlockMat& B, const double& beta){
    int id,ida,idb,i,j,k;
    BlockSubMat *Btmp;

    // if beta == 0 then initalise (*this) matrix by  deleting row blocks
    if(!beta)
      Reset_Rows();

    switch(formB){
    case RowMajor:
      switch(formA){
      case RowMajor:
    // multiply every entry in row of B by A matrix and put in C[i]
    for(i=0; i < A._max_rowblk; ++i)
      for(j=0;j < A._Row[i]._Block.size(); ++j){
        id = A._Row[i]._Block[j]->get_id();
        _Row[i].geMxM(formA,formB,alpha,A._Row[i]._Block[j][0],
              B._Row[id],beta);
      }
    break;
      case ColMajor:
    // multiply every entry in row of B by A matrix and put in C[id]
    for(i=0; i < A._max_rowblk; ++i)
      for(j=0;j < A._Row[i]._Block.size(); ++j){
        id = A._Row[i]._Block[j]->get_id();
        _Row[id].geMxM(formA,formB,alpha,A._Row[i]._Block[j][0],
               B._Row[i],beta);
      }
    break;
      default:
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxM",
            "matrix format not specified for matrix A");
    break;
      }
      break;
    case ColMajor:
      switch(formA){
      case RowMajor:
    // Can not use BlockVec since Block row ordering is not consistent
    for(i=0; i < A._max_rowblk; ++i)
      for(j=0;j < B._max_rowblk; ++j)

        for(k = 0; k < A._max_colblk; ++k){
          if(((ida=A._Row[i]._entries[k])+1)
         &&((idb=B._Row[j]._entries[k])+1))
        if((id = _Row[i]._entries[j])+1) // add local contribution
          _Row[i]._Block[id]->geMxM(formA,formB,alpha,
          A._Row[i]._Block[ida][0],B._Row[j]._Block[idb][0],beta);
        else{ // declare and put in matrix
          Btmp = new BlockSubMat(j,A._Row[i]._Block[ida]->get_rows(),
                    B._Row[j]._Block[idb]->get_rows());
          Btmp->geMxM(formA,formB,alpha,A._Row[i]._Block[ida][0],
               B._Row[j]._Block[idb][0],beta);

          _Row[i]._entries[j] = _Row[i]._Block.size();
          _Row[i]._Block.push_back(Btmp);
        }
        }
    break;
      case ColMajor:
    // Can not use BlockVec since Block row ordering is not consistent
    for(i=0; i < A._max_colblk; ++i)
      for(j=0;j < B._max_rowblk; ++j)

        for(k = 0; k < A._max_rowblk; ++k){
          if(((ida=A._Row[k]._entries[i])+1)
         &&((idb=B._Row[j]._entries[k])+1))
        if((id = _Row[i]._entries[j])+1) // add local contribution
          _Row[i]._Block[id]->geMxM(formA,formB,alpha,
          A._Row[k]._Block[ida][0],B._Row[j]._Block[idb][0],beta);
        else{ // declare and put in matrix
          Btmp = new BlockSubMat(j,A._Row[k]._Block[ida]->get_cols(),
                    B._Row[j]._Block[idb]->get_rows());
          Btmp->geMxM(formA,formB,alpha,A._Row[k]._Block[ida][0],
               B._Row[j]._Block[idb][0],beta);

          _Row[i]._entries[j] = _Row[i]._Block.size();
          _Row[i]._Block.push_back(Btmp);
        }
        }
    break;
      default:
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxM",
            "matrix format not specified for matrix A");
    break;
      }
      break;
    default:
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::geMxM",
              "matrix format not specified for matrix B");
      break;
    }

    if(!beta)
      setup_offset();

    return *this;
  }


  // (this) = a * b
  BlockMat& BlockMat::MxM(const BlockMat& a, const  BlockMat& b){
    int id,i,j;

    if((a._max_colblk != b._max_rowblk)||(a._max_rowblk != _max_rowblk)||
       (b._max_colblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::mult",
              "Matrix rows/cols are not compatible");

    // initalise (*this) matrix - delete row blocks
    Reset_Rows();

    // multiply every entry in row in each entry
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < a._Row[i]._Block.size(); ++j){
    id = a._Row[i]._Block[j]->get_id();
    _Row[i].MxMpM(a._Row[i]._Block[j][0],b._Row[id]);
      }

    setup_offset();
    return *this;
  }

  // (this) = a^T * b
  BlockMat& BlockMat::MtxM(const BlockMat& a, const BlockMat& b){
    int id,i,j;

    if((a._max_rowblk != b._max_rowblk)||(a._max_colblk != _max_rowblk)||
       (b._max_colblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::mult",
              "Matrix rows/cols are not compatible");

    // initalise (*this) matrix - delete row blocks
    Reset_Rows();

    // multiply every entry in row in each entry
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < a._Row[i]._Block.size(); ++j){
    id = a._Row[i]._Block[j]->get_id();
    _Row[id].MtxMpM(a._Row[i]._Block[j][0],b._Row[i]);
      }

    setup_offset();
    return *this;
  }


  // (this) = a * b^t
  BlockMat& BlockMat::MxMt(const BlockMat& a, const BlockMat& b){
    int id,ida,idb,i,j,k;

    if((a._max_colblk != b._max_colblk)||(a._max_rowblk != _max_rowblk)||
       (b._max_rowblk != _max_colblk))
      ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::mult",
              "Matrix rows/cols are not compatible");


    // initalise (*this) matrix - delete row blocks
    Reset_Rows();


    // Can not use BlockVec since Block row ordering is not consistent
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < b._max_rowblk; ++j)

    for(k = 0; k < a._max_colblk; ++k){
      if(((ida=a._Row[i]._entries[k])+1)&&((idb=b._Row[j]._entries[k])+1))
        if((id = _Row[i]._entries[j])+1) // add local contribution
          _Row[i]._Block[id]->MxMtpM(a._Row[i]._Block[ida][0],
                       b._Row[j]._Block[idb][0]);
        else{ // declare and put in matrix
          BlockSubMat *B;
          B = new BlockSubMat(j,a._Row[i]._Block[ida]->get_rows(),
                b._Row[j]._Block[idb]->get_rows());
          B->MxMtpM(a._Row[i]._Block[ida][0],b._Row[j]._Block[idb][0]);
          _Row[i]._entries[j] = _Row[i]._Block.size();
          _Row[i]._Block.push_back(B);
        }
    }

    setup_offset();
    return *this;
  }

  // invert diagonal matrices
  BlockMat& BlockMat::invert_diag(){
    int i,id;

    for(i = 0; i < _max_rowblk; ++i){
      if(_Row[i]._Block.size() != 1)
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::invert","matrix not diagonal");

      if((id = _Row[i]._entries[i])+1)
    _Row[i]._Block[id]->invert();
      else
    ErrorUtil::Error(ErrorUtil::ewarning,"BlockMat::invert","No diagonal component");
    }
    return *this;
  }


  double BlockMat::AmAt(){
    int i,j,k,l,nr,nc;
    double  sum = 0;
    double *mat1,*mat2;

    for(i = 0; i < _max_rowblk; ++i)
      for(j = 0; j < _max_colblk; ++j)
    if(mat1 = get_mat(i,j,nr,nc))
      if(mat2 = get_mat(j,i,nr,nc))
        for(k = 0; k < nr; ++k)
          for(l = 0; l < nc; ++l)
        sum += mat1[k*nc+l] - mat2[l*nc+k];
    return sum;
  }
}

