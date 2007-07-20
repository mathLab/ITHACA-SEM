//#ifndef H_NEKLINSYS
//#define H_NEKLINSYS
//
//#include <LibUtilities/NekLinAlg.hpp>
//#include <LibUtilities/Lapack.h>
//#include <LibUtilities/NekMatrix.hpp>
//#include <LibUtilities/NekVector.hpp>
//
//#include <boost/shared_ptr.hpp> 
//
//namespace LibUtilities
//{
//    namespace NekLinAlg
//    {
//
//    class NekLinSys
//    {
//    private:
//        boost::shared_ptr<NekMatrix<double> > _matrix;
//        
//        bool    _factored;      ///< Flag to identify if matrix is factored
//        int    *_ipiv;          ///< Pivoting array   
//        int     _lda;           ///< leading diagonal of matrix
//        int     _bwidth;        ///< Bandwdith for positive banded matrix \n
//        ///< Upper sug diagonals plus one for general
//        ///< banded matrix
//        int     _ldiag;         ///< Low sub diagonals for general banded matrix
//        double* _packed_matrix; ///< Inverse/Factorisation of Matrix pointer 
//        
//        
//        
//        // --------------------------------------------------------------------      
//        // Private Methods
//        // --------------------------------------------------------------------
//        
//        inline void Initiate()
//        {
//        _factored = false;
//        _ipiv = (int*)NULL;
//        _packed_matrix = (double *)NULL;
//        _bwidth = 0; _ldiag = 0;
//        }
//        
//        // declare the memory of _packed_matrix depending upon its definition
//        void SetMemPackedMatrix()
//        {
//        
//        if(_packed_matrix) 
//        {
//            return;
//        }
//        
//        ASSERTL0((_lda > 0), "NekLinSys::SetMemPackedMatrix","_lda not defined ");
//        
//        switch(_matrix->matform()){
//        case Symmetric:     case Symmetric_Positive:
//            _packed_matrix = new double [_lda*(_lda+1)/2];
//            Vmath::zero(_lda*(_lda+1)/2,_packed_matrix,1);
//            break;
//        case Symmetric_Positive_Banded:
//            ASSERTL0((_bwidth > 0 ), "NekLinSys::SetMemPackedMatrix",
//                 "_bwidth  not set");
//            _packed_matrix = new double [_lda*_bwidth];
//            Vmath::zero(_lda*_bwidth,_packed_matrix,1);
//            break;
//        case General_Banded:
//            ASSERTL0((_bwidth > 0 )&&(_ldiag > 0),
//                 "NekLinSys::SetMemPackedMatrix","_bwidth or _ldiag is  not set");
//            _packed_matrix = new double [_lda*2*(_ldiag+_bwidth)];
//            Vmath::zero(_lda*2*(_ldiag+_bwidth),_packed_matrix,1);
//            break;
//        case General_Full:
//            _packed_matrix = new double [_lda*_lda];
//            Vmath::zero(_lda*_lda,_packed_matrix,1);
//            break;
//        } 
//        }
//      
//        /** \brief factorise matrix depending upon definition stored in _mat_form.
//        
//        Options for mat_from are: 
//        
//        - _mat_form = Symmetric-Positive implies Cholesky factorization using
//        Lapack::dpptrf.
//        
//        - _mat_form = Symmetric implies Factorisation using Bunch-Kaufman
//        pivoting using Lapack::dsptrf.
//        
//        - _mat_form = Symmetric-Positive-Banded implies lower diagonal banded
//        cholesky factorisation using Lapack::dpbtrf.
//        
//        - _mat_form = General-Banded implies factoring using Lapack::dgbtrf.
//        
//        - _mat_form = General-Full implies factoring using Lapack::dgetrf.
//        
//        */
//        void Factor()
//        {
//        int info;
//        
//        if(_factored)
//        {
//            return;
//        }
//    
//        ASSERTL1(_lda, "NekLinSys::Factor","_lda is not set");
//        
//        if(!_packed_matrix)
//        {
//            FillPackedMatrix();
//        }
//        
//        switch(_matrix->matform())
//        {
//        case Symmetric:
//            _ipiv = new int[_lda];
//            Lapack::dsptrf('L',_lda,_packed_matrix,_ipiv,info);
//            ASSERTL0(info==0, "NekLinSys::Factor","matrix did not factor");
//            break;
//        case Symmetric_Positive:
//            Lapack::dpptrf('L', _lda, _packed_matrix, info);
//            ASSERTL0(info==0, "NekLinSys::Factor","matrix did not factor");
//            break;
//        case Symmetric_Positive_Banded:
//            Lapack::dpbtrf('L',_lda,_bwidth-1,_packed_matrix,_bwidth,info);
//            ASSERTL0(info==0, "NekLinSys::Factor","matrix did not factor");
//            break;
//        case General_Banded:
//            _ipiv = new int[_lda];
//            Lapack::dgbtrf(_lda,_lda,_ldiag,_bwidth-1,_packed_matrix,
//                   2*_ldiag+_bwidth,_ipiv,info);
//            ASSERTL0(info==0, "NekLinSys::Factor","matrix did not factor");
//            break;
//        case General_Full:
//            _ipiv = new int[_lda];
//            Lapack::Dgetrf(_lda,_lda,_packed_matrix,_lda,_ipiv,info);     
//            ASSERTL0(info==0, "NekLinSys::Factor","matrix did not factor");
//            break;
//        }
//        
//        _factored = true;    
//        }
//      
//      
//        /** \brief Fill _packed_matrix from _matrix 
//        
//        Note: Assume input matrix is stored in RowMajor 'C' type format. The
//        packed matrix however is stored in Column Major Fortran type format
//        so that the appropriate Lapack routines can be called. This is only
//        important for the General Matrix forms.
//        
//        */
//        
//        void FillPackedMatrix()
//        {
//        int i,j,cnt;
//        
//        ASSERTL1(_lda, "NekLinSys::FillPackedMatrix",
//             "lda is not set");
//        
//        // check to see if memory is declared and if not setup
//        if(!_packed_matrix)
//        {
//            SetMemPackedMatrix();
//        }
//    
//        switch(_matrix->matform())
//        {
//        case Symmetric_Positive: case Symmetric:
//            // store matrix in symmetric form 
//            cnt = 0;
//            for(i = 0; i < _lda; ++i)
//            {
//            for(j = i; j < _lda; ++j)
//            {
//                _packed_matrix[cnt++] = (*_matrix)(i,j);
//            }
//            }
//            break;
//        case Symmetric_Positive_Banded:
//            // store matrix in symmetric banded 
//            cnt = 0;
//            for(i = 0; i < _lda; ++i)
//            {
//            for(j  = i; j < i+_bwidth; ++j)
//            {
//                _packed_matrix[cnt++] = (*_matrix)(i,j);
//            }
//            }
//            break;
//        case General_Banded:
//            for(i = 0; i < _lda; ++i)
//            {
//            for(j  = max(i-_ldiag,0); j < i+_bwidth; ++j)
//            {
//                _packed_matrix[j*(2*_ldiag+_bwidth)+_ldiag+_bwidth-1+(i-j)]
//                = (*_matrix)(i,j);
//            }
//            }
//            break;
//        case General_Full:
//            // pack full matrix transposing matrix 
//            for(i = 0; i < _lda; ++i)
//            {
//            for(j = 0; j < _lda; ++j)
//            {
//                _packed_matrix[j+i*_lda] = (*_matrix)(i,j);
//            }
//            }
//            
//            break;
//        }
//        }
//        
//        
//    public:
//      
//        NekLinSys(int rows): _matrix(new NekMatrix<double>(rows,rows)){
//        Initiate();
//        _lda = rows;  
//        }
//        
//        NekLinSys(int rows, int columns): _matrix(new NekMatrix<double>(rows,columns)){
//        Initiate();
//        _lda = rows;      
//        }
//        
//        NekLinSys(MatrixForm matform, int rows): _matrix(new NekMatrix<double>(matform,rows,rows)){
//        Initiate();
//        _lda = rows;  
//      }
//      
//      
//      NekLinSys(boost::shared_ptr<NekMatrix<double> > mat): _matrix(mat){
//        Initiate();    
//    _lda = _matrix->rows();  
//      }
//      
//      
//      
//      // --------------------------------------------------------------------      
//      // Access (get and set) of private data
//      // --------------------------------------------------------------------
//      
//      /// return private data for _matform
//      MatrixForm matform() const { return _matrix->matform(); }
//      
//      /// return private data for _lda
//      int get_lda     () const { return _lda;}    
//      
//      /// return private data for lda
//      int get_bwidth    () const { return _bwidth;}
//      
//      /// set private data: bwidth
//      void set_bwidth (int val) {_bwidth = val;}
//      
//      /// return private data for _ldiag
//      int get_ldiag    () const { return _ldiag;}
//      
//      /// set private data: bwidth
//      void set_ldiag (int val) {_ldiag = val;}
//      
//      /// dump the matrix to file out
//      void dump_matrix           (FILE *out);
//      
//      /// show the matrix structure in file out
//      void show_matrix_structure (FILE *out);
//
//      inline int size() const { return _matrix->size(); }
//      inline int rows() const { return _matrix->rows(); }
//      inline int cols() const { return _matrix->cols(); }
//      
//      inline double operator()(int i, int j) const { 
//    return (*_matrix)(i,j); 
//      }
//      
//      inline double& operator()(int i, int j){ 
//    return (*_matrix)(i,j); 
//      }
//      
//      inline SharedArray<double> getdata() const {return _matrix->getdata();}
//      inline double * getdataptr() const {return _matrix->getdataptr();}
//      inline bool unique() const {return _matrix->unique();}
//      
//      
//      /**  \brief  evaulate \f$u \leftarrow _packed_matrix^{-1} u\f$ for \a nrhs solves
//       Options for mat_from are: 
//       
//       - _mat_form = Symmetric-Positive implies Cholesky back solve using
//       Lapack::dpptrs.
//       
//       - _mat_form = Symmetric implies Back solve from  using Lapack::dsptrs.
//       
//       - _mat_form = Symmetric-Positive-Banded implies lower diagonal banded
//       cholesky backsole using  Lapack::dpbtrs.
//       
//       - _mat_form = General-Banded implies LU back solve using Lapack::dgbtrs.
//       
//       - _mat_form = General-Full implies LU back solve using Lapack::dgetrs.
//       
//      */
//      
//      
//      void Solve(double *u, int nrhs){
//    int info;
//    
//    if(!_factored) Factor();
//    
//    switch(_matrix->matform()){
//    case Symmetric:
//      Lapack::dsptrs('L',_lda,nrhs,_packed_matrix,_ipiv,u,_lda,info);
//      ASSERTL0(info==0, "NekLinSys::Solve","matrix did not solve");
//      break;
//    case Symmetric_Positive:
//      Lapack::dpptrs('L', _lda,nrhs,_packed_matrix,u,_lda,info);
//      ASSERTL0(info==0, "NekLinSys::Solve","matrix did not solve");
//      break;
//    case Symmetric_Positive_Banded:
//      Lapack::dpbtrs('L',_lda,_bwidth-1,nrhs,_packed_matrix,_bwidth,u,_lda,info);
//      ASSERTL0(info==0, "NekLinSys::Solve","matrix did not solve");
//      break;
//    case General_Banded:
//      Lapack::dgbtrs('N',_lda,_ldiag,_bwidth-1,nrhs,_packed_matrix,
//             2*_ldiag+_bwidth,_ipiv,u,_lda,info);
//      ASSERTL0(info==0, "NekLinSys::Solve","matrix did not solve");
//      break;
//    case General_Full:
//      Lapack::Dgetrs('N',_lda,nrhs,_packed_matrix,_lda,_ipiv,u,_lda,info);     
//      ASSERTL0(info==0, "NekLinSys::Solve","matrix did not solve");
//      break;
//    }
//      }
//      
//      
//      double L2ConditionNo(){
//    double *er = new double [_lda];
//    double *ei = new double [_lda];
//    double max,min;
//    
//    EigenValues(er,ei,(double *)NULL);
//    
//    Vmath::vmul (_lda,er,1,er,1,er,1);
//    Vmath::vmul (_lda,ei,1,ei,1,ei,1);
//    Vmath::vadd (_lda,er,1,ei,1,er,1);
//    
//    max = sqrt(er[Vmath::imax(_lda,er,1)]);
//    min = sqrt(er[Vmath::imin(_lda,er,1)]);
//    
//    if(min < 1e-12){ // if min < 1e-12 find second smallest ev
//      fprintf(stderr,"Min ev < 1e-12 using second ev\n");
//      er[Vmath::imin(_lda,er,1)] += max;
//      min = sqrt(er[Vmath::imin(_lda,er,1)]);
//    }
//    
//    delete[] er;
//    delete[] ei;
//    return max/min;
//      }
//      
//      double MaxEigenValue(){
//    double *er = new double[_lda];
//    double *ei = new double[_lda];
//    double max;
//    
//    EigenValues(er,ei,(double *)NULL);
//    
//    Vmath::vmul(_lda,er,1,er,1,er,1);
//    Vmath::vmul(_lda,ei,1,ei,1,ei,1);
//    Vmath::vadd(_lda,er,1,ei,1,er,1);
//    
//    max = sqrt(er[Vmath::imax(_lda,er,1)]);
//    delete[] er;
//    delete[] ei;
//    return max;
//      }
//      
//      /** 
//      Return the null space of the Packed Matrix 
//      where the null space is defined by the number of eigenvalues
//      with magnitude smaller than tol.
//      */
//      
//      int NullSpaceDim(const double tol){
//    double *er = new double [_lda];
//    double *ei = new double [_lda];
//    int i,ndim; 
//    
//    EigenValues(er,ei,(double *)NULL);
//    
//    Vmath::vmul (_lda,er,1,er,1,er,1);
//    Vmath::vmul (_lda,ei,1,ei,1,ei,1);
//    Vmath::vadd (_lda,er,1,ei,1,er,1);
//    Vmath::vsqrt(_lda,er,1,er,1);
//    
//    ndim = 0;
//    for(i = 0; i < _lda; ++i)
//      if(er[i] < tol)
//        ++ndim;
//    
//    delete[] er;
//    delete[] ei;
//    
//    return ndim;
//      }
//      
//      
//      /** 
//      Write the eigenvalues of matrix _packed_matrix to file \a file
//      */
//      void EigenValues(const char file[]){
//    int i;
//    double *er = new double [_lda];
//    double *ei = new double [_lda];
//    FILE *fp;
//    
//    fp = fopen(file,"w");
//    
//    EigenValues(er,ei,(double *)NULL);
//    
//    fprintf(fp,"# Real Imag Magnitude\n");
//    for(i = 0; i < _lda; ++i)
//      fprintf(fp,"%lg %lg %lg \n",er[i],ei[i],sqrt(er[i]*er[i]+ei[i]*ei[i]));
//    
//    fclose(fp);
//    delete[] er;
//    delete[] ei;
//      }
//      
//      /** 
//      Determine the eigen specturm of the Packed matrix by  returning its
//      read values in \a er and its imaginary values in \a ei. 
//      
//      if evecs is given then the eigenvectors are also returned. 
//      This option is currently only setup for a General_Full matrix.
//      
//      */
//      
//      void EigenValues(double *er, double *ei, double *evecs){
//    double dum;
//    int    info;
//    
//    switch(_matrix->matform()){
//    case Symmetric_Positive: case Symmetric:
//      {
//        double *work  = new double [3*_lda];
//        Vmath::zero(_lda,ei,1);
//        Lapack::dspev('N','L',_lda,_packed_matrix,er,&dum,1,work,info);
//        ASSERTL0(info==0, "NekLinSys::EigenValues","info is not zero");
//        delete[] work;
//        break;
//      }
//    case Symmetric_Positive_Banded:
//      {
//        double *work  = new double [3*_lda];
//        Vmath::zero(_lda,ei,1);
//        Lapack::dsbev('N','L',_lda,_bwidth-1,_packed_matrix,_bwidth,er,
//              &dum,1,work,info);
//        ASSERTL0(info==0, "NekLinSys::EigenValues","info is not zero");
//        delete[] work;
//        break;
//      }
//    case General_Banded:
//      ErrorUtil::Error(fatal,"GMatrix::Spectrum","Eigenvalue evaluation "
//               "for genaral baneded matrix needs coding");
//      break;
//    case General_Full:
//      {
//        double *work  = new double [4*_lda];
//        if(evecs)
//          Lapack::dgeev('N','V',_lda,_packed_matrix,_lda,er,ei,&dum,1,
//                evecs,_lda,work,4*_lda,info);
//        else{
//          double dum1;
//          Lapack::dgeev('N','N',_lda,_packed_matrix,_lda,er,ei,&dum,1,
//                &dum1,1,work,4*_lda,info);
//        }      
//        ASSERTL0(info==0, "NekLinSys::EigenValues","info is not zero");
//        delete[] work;
//        break;
//      }
//    }
//      }      
//    };
//  } // end of namespace
//} // end of namespace
//
//#endif
