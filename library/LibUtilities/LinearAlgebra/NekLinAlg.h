#ifndef H_NEKLINALG
#define H_NEKLINALG

#include <xpr1.h> //from MET library
#include <xpr2.h> //from MET library

#include <LibUtilities/Vmath.h>
#include <LibUtilities/ErrorUtil.h>

#include <boost/shared_array.hpp> 
#include <boost/static_assert.hpp>
#include <boost/call_traits.hpp>


namespace LibUtilities{
  namespace NekLinAlg{
    
    /// NekPoint Class
    template<typename DataType, int dim>
    class NekPoint{
    protected:
      DataType _data[dim];
      
    public:
      
      typedef typename boost::call_traits<DataType>::param_type ParamDataType;
      typedef typename boost::call_traits<DataType>::value_type ValueDataType;
      typedef typename boost::call_traits<DataType>::const_reference ConstReferenceDataType;
      typedef typename boost::call_traits<DataType>::reference ReferenceDataType;

      /// Default constructor
      NekPoint(){
	for(int i=0;i<dim;++i)
	  _data[i] = 0;
      }
      
      /// Constructor
      NekPoint(const DataType &a){
	Vmath::fill(dim,a,_data,1);
      }

      /// Copy Constructor
      NekPoint(const NekPoint<DataType,dim>& rhs){
	for(int i = 0; i<dim; ++i) 
	  _data[i] = rhs._data[i];
      }

      /// Destructor
      virtual ~NekPoint(){
      }

      /// Dimension       
      inline int dimension() const { return dim; }

      inline DataType * getdataptr() const {return &(_data[0]);}
          
      inline ReferenceDataType operator()(int i) { 
	ASSERTL1( (i>=0) && (i<dim), "NekPoint::operator", "Invalid access to _data via parenthesis operator");
	return _data[i];
      }
   
      inline ConstReferenceDataType operator()(int i) const  { 
	ASSERTL1( (i>=0) && (i<dim), "NekPoint::operator", "Invalid access to _data via parenthesis operator");
	return _data[i];
      }

      /// Assignment Operator
      NekPoint& operator  = (const NekPoint<DataType,dim>& rhs){
	for(int i = 0; i<dim; ++i) 
	  _data[i] = rhs._data[i];

	return *this;
      }

      ConstReferenceDataType x() const {
	BOOST_STATIC_ASSERT(dim >= 1);
	return (*this)(0);
      }

      ConstReferenceDataType y() const {
	BOOST_STATIC_ASSERT(dim >= 2);
	return (*this)(1);
      }

      ConstReferenceDataType z() const {
	BOOST_STATIC_ASSERT(dim >= 3);
	return (*this)(2);
      }

      ReferenceDataType x() {
	BOOST_STATIC_ASSERT(dim >= 1);
	return (*this)(0);
      }

      ReferenceDataType y() {
	BOOST_STATIC_ASSERT(dim >= 2);
	return (*this)(1);
      }

      ReferenceDataType z() {
	BOOST_STATIC_ASSERT(dim >= 3);
	return (*this)(2);
      }

      /// Arithmetic Routines

      // Unitary operators
      NekPoint operator-() const{
	NekPoint<DataType,dim> temp(*this);
	for(int i=0;i<dim;++i)
	  temp(i) = -temp(i);
	return temp;
      }

      // Binary Operators
      NekPoint& operator += (const NekPoint<DataType,dim>& rhs){
	for(int i=0;i<dim;++i)
	  _data[i] += rhs._data[i];
	return *this;
      }

      NekPoint& operator -= (const NekPoint<DataType,dim>& rhs){
	for(int i=0;i<dim;++i)
	  _data[i] -= rhs._data[i];
	return *this;
      }

      NekPoint& operator *= (const DataType& rhs){
	for(int i=0;i<dim;++i)
	  _data[i] *= rhs;
	return *this;
      }


      NekPoint operator+(const NekPoint<DataType,dim>& P0){
	NekPoint<DataType,dim> temp(P0);
	for(int i=0;i<dim;++i)
	  temp += *this;
	return temp;
      }

      
      NekPoint operator-(const NekPoint<DataType,dim>& P0){
	NekPoint<DataType,dim> temp(P0);
	for(int i=0;i<dim;++i)
	  temp -= *this;
	return temp;
      }
   
    };


    /////////////////////////////////////////////////////////////////////////////////////////    

    enum VectorForm{
      UninitialisedVector=0,        ///< Uninitialized vector  
      ZeroVector,                   ///< Vector containing all zeros
      ConstantVector,               ///< Vector containing a single value
      GeneralVector                 ///< General entries
    };


    /// NekVector Class
    template<typename DataType>
    class NekVector : public Dim1<DataType,NekVector<DataType> >{
    protected:
      int _size;                             ///< length of the vector
      VectorForm _vecform;                   ///< type of the vector
      boost::shared_array<DataType> _data;   ///< smart pointer to data
      
    public:

      // //////////////////////////
      /// Default Constructor
      //  Comments: 
      //
      // /////////////////////////

      NekVector() : Dim1<DataType,NekVector<DataType> >(), _size(0), _vecform(UninitialisedVector){	
	
      }

      // //////////////////////////
      /// Constructor B
      //  Comments: 
      //
      // /////////////////////////

      NekVector(int size, VectorForm vecform) : Dim1<DataType,NekVector<DataType> >(),
						_size(size), _vecform(vecform){
	ASSERTL0( ((_size>0)&&(_vecform!=UninitialisedVector)) || ((_size==0)&&(_vecform==UninitialisedVector)), 
		  "NekVector::NekVector", "Constructor B Error");
	
	switch(_vecform){
	case ConstantVector:
	  _data.reset(new DataType[1]);
	  break;
	case GeneralVector:
	  _data.reset(new DataType[_size]);
	  break;
	}
      }

      // //////////////////////////
      /// Constructor C
      //  Comments: 
      //
      // /////////////////////////

      NekVector(int size) : Dim1<DataType,NekVector<DataType> >(),
			    _size(size), _data(new DataType[_size]), _vecform(GeneralVector){
	ASSERTL0( (_size>0), "NekVector::NekVector", "Constructor C Error");
      }
      

      // //////////////////////////
      /// Constructor D
      //  Comments: 
      //
      // /////////////////////////

      
      NekVector(int size, boost::shared_array<DataType> ptr) : Dim1<DataType,NekVector<DataType> >(),
							       _size(size), _data(ptr), _vecform(GeneralVector){
	ASSERTL0( (_size>0), "NekVector::NekVector", "Constructor D Error");
     }


      // //////////////////////////
      /// Copy Constructor
      //  Comments: 
      //
      // /////////////////////////
      
      NekVector(const NekVector<DataType>& rhs) : Dim1<DataType,NekVector<DataType> >(), 
						  _size(rhs._size), _vecform(rhs._vecform){ 
	ASSERTL0( (_size>0), "NekVector::NekVector", "Copy Constructor Error");

	switch(_vecform){
	case ConstantVector:
	  _data.reset(new DataType[1]);
	  _data[0] = _rhs._data[0];
	  break;
	case GeneralVector:
	  _data.reset(new DataType[_size]);
	  assignFrom(rhs);
	  break;
	}
      }

      // //////////////////////////
      /// Destructor
      //  Comments: 
      //
      // /////////////////////////

      
      virtual ~NekVector(){}


      // //////////////////////////
      /// Initialise A
      //  Comments: 
      //
      // /////////////////////////
      
      void Initialise(int size, VectorForm vecform){
	ASSERTL0(_vecform==UninitialisedVector, "NekVector::Initialize A", "Can only initialize an uninitialized vector");
	ASSERTL0( (size>0) && (vecform!=UninitialisedVector), "NekVector::Initialize A", "Can only initialize to a positive size");

	_vecform = vecform;
	_size = size;

	switch(_vecform){
	case ConstantVector:
	  _data.reset(new DataType[1]);
	  break;
	case GeneralVector:
	  _data.reset(new DataType[_size]);
	  break;
	}
      }

      // //////////////////////////
      /// Initialise B
      //  Comments: 
      //
      // /////////////////////////

      void Initialise(int size, boost::shared_array<DataType> ptr){
	ASSERTL0(_vecform==UninitialisedVector, "NekVector::Initialize B", "Can only initialize an uninitialized vector");
	ASSERTL0( (size>0) && (vecform!=UninitialisedVector), "NekVector::Initialize B", "Can only initialize to a positive size");

	_vecform = GeneralVector;
	_size = size;
	_data.reset(ptr); 
      }


      // //////////////////////////
      /// Operator= A
      //  Comments: 
      //
      // /////////////////////////

      template <typename X> NekVector<DataType>& operator=(const Xpr1<DataType,X>& rhs) {
	ASSERTL0(_vecform==GeneralVector, "NekVector::operator= A", "Can only assign to a GeneralVector");

	return assignFrom(rhs);
      }

      // //////////////////////////
      /// Operator= B
      //  Comments: 
      //
      // /////////////////////////
      
      template <typename V> NekVector<DataType>& operator=(const Dim1<DataType,V>& rhs) {
	ASSERTL0(_vecform==GeneralVector, "NekVector::operator= B", "Can only assign to a GeneralVector");
	ASSERTL0(_size==rhs.size(), "NekVector::operator= B", "Sizes must match");
	
	return assignFrom(rhs);
      }

      // //////////////////////////
      /// Operator= C
      //  Comments: 
      //
      // /////////////////////////
      
      NekVector<DataType>& operator=(const NekVector<DataType>& rhs) { 
	ASSERTL0(_vecform==GeneralVector, "NekVector::operator= C", "Can only assign to a GeneralVector");
	ASSERTL0(_size==rhs.size(), "NekVector::operator= B", "Sizes must match");
	
	return assignFrom(rhs); 
      }

      // //////////////////////////
      /// Operator= D
      //  Comments: 
      //
      // /////////////////////////
      
      NekVector<DataType>& operator=(DataType rhs) { 
	ASSERTL0((_vecform==GeneralVector)||(_vecform==ConstantVector), "NekVector::operator= D", 
		 "Can only assign to a ConstantVector or a GeneralVector");
	switch(_vecform){
	case ConstantVector:
	  _data[0] = rhs;
	  return *this;
	  break;
	case GeneralVector:
	  return assignFrom(rhs);
	  break;
	} 
      }
      

      // //////////////////////////
      /// Operator() A
      //  Comments: 
      //
      // /////////////////////////

      inline DataType operator()(int i) const { 
	ASSERTL0(_vecform!=UninitialisedVector, "NekVector::operator() A", "Invalid query of uninitialised vector");	
	ASSERTL1( (i>=0) && (i<_size), "NekVector::operator", "Invalid access to _data via parenthesis operator");
	
	switch(_vecform){
	case ZeroVector:
	  return 0;
	  break;
	case ConstantVector:
	  return _data[0];
	  break;
	case GeneralVector:
	  return _data[i]; 
	  break;
	}
      }

      // //////////////////////////
      /// Operator() B
      //  Comments: 
      //
      // /////////////////////////

      inline DataType& operator()(int i){ 
	ASSERTL0((_vecform==GeneralVector)||(_vecform==ConstantVector), "NekVector::operator= D", 
		 "Can only assign to a ConstantVector or a GeneralVector");	
	ASSERTL1( (i>=0) && (i<_size), "NekVector::operator", "Invalid access to _data via parenthesis operator");

	switch(_vecform){
	case ConstantVector:
	  ASSERTL0(i==0, "NekVector::operator= D", "Can only assign to a ConstantVector zeroeth element");	
	  return _data[0];
	  break;
	case GeneralVector:
	  return _data[i]; 
	  break;
	}
      }
     
      /// Necessary information methods

      inline boost::shared_array<DataType> getdata() const {return _data;}
      inline DataType * getdataptr() const {return _data.get();}
      inline bool unique() const {return _data.unique();};
      inline int size() const { return _size; }
      inline VectorForm vecform() const { return _vecform;}

      /// Extra Operators
    }; 


    /////////////////////////////////////////////////////////////////////////////////////////    

    enum MatrixForm{
      UninitializedMatrix=0,        ///< Uninitialized system                
      RectangleMatrix,              ///< Rectangular System
      ZeroMatrix,                   ///< Matrix containing all zeros
      DiagMatrix,                   ///< Square Matrix with only diagonal entries
      Symmetric,                    ///< Square Symmetric matrix
      Symmetric_Positive,           ///< Symmetric positive definite matrix
      Symmetric_Positive_Banded,    ///< Symmetric positive definite banded matrix
      General_Full,                 ///< Square General full 
      General_Banded,               ///< General Banded System
    };

    /// NekMatrix Class
    template <typename DataType>
    class NekMatrix : public Dim2<DataType,NekMatrix<DataType> > {

    private:
      int _rowsize, _columnsize, _totalsize;
      MatrixForm _matform;
      boost::shared_array<DataType> _data;
      boost::shared_array<DataType*> _mat;

    public:

      // //////////////////////////
      /// Default Constructor
      //  Comments: 
      //
      // /////////////////////////

      NekMatrix(): Dim2<DataType,NekMatrix<DataType> >(),
		   _rowsize(0), _columnsize(0), _totalsize(0), 
		   _matform(UninitializedMatrix){
     
      }


      // //////////////////////////
      /// Constructor A
      //  Comments: 
      //
      // /////////////////////////

      NekMatrix(int rows) : Dim2<DataType,NekMatrix<DataType> >(),
			    _rowsize(rows), _columnsize(rows), _totalsize(_rowsize*_columnsize), 
			    _data(new DataType[_totalsize]), _mat(new DataType*[_rowsize]), _matform(General_Full) {
	for (int i=0; i<_rowsize; ++i) 
	  _mat[i] = _data.get() + i*_columnsize;
      }


      // //////////////////////////
      /// Constructor B
      //  Comments: 
      //
      // /////////////////////////

      NekMatrix(MatrixForm matform, int rows) : Dim2<DataType,NekMatrix<DataType> >(),
						_rowsize(rows), _columnsize(rows), _totalsize(_rowsize*_columnsize), 
						_matform(matform){
	int i;

	switch(_matform){
	case RectangleMatrix: case General_Full:
	  _data.reset(new DataType[_totalsize]);
	  _mat.reset(new DataType*[_rowsize]);
	  for (i=0; i<_rowsize; ++i) 
	    _mat[i] = _data.get() + i*_columnsize;
	  break;
	  
	case ZeroMatrix:
	  // No memory allocated
	  break;
	  
	case DiagMatrix:
	  _data.reset(new DataType[_rowsize]);
	  break;
	  
	case Symmetric: 
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  _data.reset(new DataType[((_rowsize+1)*(_columnsize))/2]);
	  break;
	  
	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");
	  
	}
      }
      
      // //////////////////////////
      /// Constructor C
      //  Comments: 
      //
      // /////////////////////////
      
      NekMatrix(int rows, int columns) : Dim2<DataType,NekMatrix<DataType> >(),
					 _rowsize(rows), _columnsize(columns), _totalsize(_rowsize*_columnsize), 
					 _data(new DataType[_totalsize]), _mat(new DataType*[_rowsize]), _matform(RectangleMatrix) {
	
	for (int i=0; i<_rowsize; ++i) 
	  _mat[i] = _data.get() + i*_columnsize;
      }
      
      // //////////////////////////
      /// Constructor D
      //  Comments: 
      //
      // /////////////////////////
      
      NekMatrix(MatrixForm matform, int rows, int columns) : 
	Dim2<DataType,NekMatrix<DataType> >(),
	_rowsize(rows), _columnsize(columns), _totalsize(_rowsize*_columnsize), 
	_data(new DataType[_totalsize]), _mat(new DataType*[_rowsize]) {

	int i;

	switch(_matform){
	case RectangleMatrix: case General_Full:
	  _data.reset(new DataType[_totalsize]);
	  _mat.reset(new DataType*[_rowsize]);
	  for (i=0; i<_rowsize; ++i) 
	    _mat[i] = _data.get() + i*_columnsize;
	  break;
	  
	case ZeroMatrix:
	  // No memory allocated
	  break;
	 
	case DiagMatrix:
	  ASSERTL0(_rowsize != _columnsize, 
		   "NekMatrix::NekMatrix", "DiagMatrix must be a square matrix");
	  
	  _data.reset(new DataType[_rowsize]);
	  break;

	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  ASSERTL0(_rowsize != _columnsize, 
		   "NekMatrix::NekMatrix", "Symmetric Matrix must be a square matrix");
	  
	  _data.reset(new DataType[((_rowsize+1)*(_columnsize))/2]);
	  break;

	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");
	}
      }
      

      // //////////////////////////
      /// Constructor E
      //  Comments: 
      //
      // /////////////////////////
      
      NekMatrix(int rows, int columns, 
		boost::shared_array<DataType> ptr) : Dim2<DataType,NekMatrix<DataType> >(),
						     _rowsize(rows), _columnsize(columns), 
						     _totalsize(_rowsize*_columnsize), 
						     _data(ptr), _mat(new DataType*[_rowsize]) {
	
	if(rows == columns)
	  _matform = General_Full;
	else
	  _matform = RectangleMatrix;

	for (int i=0; i<_rowsize; ++i) 
	  _mat[i] = _data.get() + i*_columnsize;
      }
      

      // //////////////////////////
      /// Copy Constructor
      //  Comments: 
      //
      // /////////////////////////
      
      NekMatrix(const NekMatrix& rhs) : Dim2<DataType,NekMatrix<DataType> >(),
					_rowsize(rhs._rowsize), _columnsize(rhs._columnsize),
					_matform(rhs._matform) {
	
	int i,j;

	switch(_matform){
	case Rectangle: case General_Full:
	  _data.reset(new DataType[_totalsize]);
	  _mat.reset(new DataType*[_rowsize]);
	  for (i=0; i<_rowsize; ++i) 
	    _mat[i] = _data.get() + i*_columnsize;

	  assignFrom(rhs);	
	  break;
	  
	case ZeroMatrix:
	  // No memory allocated
	  break;
	  
	case DiagMatrix:
	  _data.reset(new DataType[_rowsize]);
	  for(i=0;i<_rowsize;i++)
	    _data[i] = rhs._data[i];
	  break;
  
	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  _data.reset(new DataType[((_rowsize+1)*(_columnsize))/2]);

	  for(i=0;i<_rowsize;i++)
	    for(j=i;j<_columnsize;j++)
	      (*this)(i,j) = rhs(i,j);

	  break;

	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");

	}
	
      }


      // //////////////////////////
      /// Destructor
      //  Comments: 
      //
      // /////////////////////////

      virtual ~NekMatrix(){}


      // //////////////////////////
      /// Initialise A
      //  Comments: 
      //
      // /////////////////////////

      void Initialise(MatrixForm matform, int rowsize){
	int i;

	ASSERTL0(_matform==UninitializedMatrix, "NekMatrix::Initialize", "Can only initialize an uninitialized matrix");
	
	_matform = matform;
	_rowsize = rowsize;
	_columnsize = rowsize;
	_totalsize = _rowsize*_columnsize;
	
	switch(_matform){
	case Rectangle: case General_Full:
	  _data.reset(new DataType[_totalsize]);
	  _mat.reset(new DataType*[_rowsize]);
	  for (i=0; i<_rowsize; ++i) 
	    _mat[i] = _data.get() + i*_columnsize;
	  break;
	  
	case ZeroMatrix:
	  // No memory allocated
	  break;
	  
	case DiagMatrix:
	  _data.reset(new DataType[_rowsize]);
	  break;
  
	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  _data.reset(new DataType[((_rowsize+1)*(_columnsize))/2]);
	  break;

	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");

	}
      }


      // //////////////////////////
      /// Initialise B
      //  Comments: 
      //
      // /////////////////////////

      void Initialise(const NekMatrix& rhs){
     
	ASSERTL0(_matform==UninitializedMatrix, "NekMatrix::Initialize", "Can only initialize an uninitialized matrix");
	
	_matform = rhs._matform;
	_rowsize = rhs._rowsize;
	_columnsize = rhs._columnsize;
	_totalsize = _rowsize*_columnsize;

	int i,j;

	switch(_matform){
	case Rectangle: case General_Full:
	  _data.reset(new DataType[_totalsize]);
	  _mat.reset(new DataType*[_rowsize]);
	  for (i=0; i<_rowsize; ++i) 
	    _mat[i] = _data.get() + i*_columnsize;

	  assignFrom(rhs);	
	  break;
	  
	case ZeroMatrix:
	  // No memory allocated
	  break;
	  
	case DiagMatrix:
	  _data.reset(new DataType[_rowsize]);
	  for(i=0;i<_rowsize;i++)
	    _data[i] = rhs._data[i];
	  break;
  
	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  _data.reset(new DataType[((_rowsize+1)*(_columnsize))/2]);

	  for(i=0;i<_rowsize;i++)
	    for(j=i;j<_columnsize;j++)
	      (*this)(i,j) = rhs(i,j);

	  break;

	default:
	  ASSERTL0(0, "NekMatrix::Initialise B", "Default in Switch Reached");
	}	
      }
      
      // //////////////////////////
      /// Operator= A
      //  Comments: 
      //
      // /////////////////////////

      template <typename X> NekMatrix<DataType>& operator=(const Xpr2<DataType,X>& rhs) {
	ASSERTL0(_vecform==General, "NekMatrix::operator= A", "Can only assign to a GeneralMatrix");
	return assignFrom(rhs);
      }
      
      // //////////////////////////
      /// Operator= B
      //  Comments: 
      //
      // /////////////////////////

      template <typename M> NekMatrix<DataType>& operator=(const Dim2<DataType,M>& rhs) {
	ASSERTL0(_vecform==General, "NekMatrix::operator= B", "Can only assign to a GeneralVector");
	ASSERTL0((_rowsize==rhs.rowsize())&&(_columnsize==rhs.columnsize()), "NekMatrix::operator= B", "Sizes must match");
	return assignFrom(rhs);
      }


      // //////////////////////////
      /// Operator() A
      //  Comments: 
      //
      // /////////////////////////

      inline DataType operator()(int i, int j) const { 
	ASSERTL1( (i>=0) && (i<_rowsize) && (j>=0) && (j<_columnsize), 
		  "NekMatrix::operator", "Invalid access to _data via parenthesis operator");

	switch(_matform){
	case Rectangle: case General_Full:
	  return _mat[i][j]; 
	  break;
	  
	case ZeroMatrix:
	  //return (DataType)Null;
	  break;
	  
	case DiagMatrix:
	  ASSERTL0(i!=j, "NekMatrix::NekMatrix", "Off Diagonal of a Diagonal Matrix Requested");
	  return _data[i];
	  break;
	  
	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  int itemp;
	  
	  if(j>=i)
	    itemp =  (_rowsize*(_rowsize+1) - ((_rowsize-i)*(_rowsize-i+1)))/2 + (j-i);
	  else
	    itemp =  (_rowsize*(_rowsize+1) - ((_rowsize-j)*(_rowsize-j+1)))/2 + (i-j);
	  
	  return _data[itemp];
	  break;
	  
	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");
  
	}
      }

      // //////////////////////////
      /// Operator() B
      //  Comments: 
      //
      // /////////////////////////
      
      inline DataType& operator()(int i, int j){ 
	ASSERTL1( (i>=0) && (i<_rowsize) && (j>=0) && (j<_columnsize), 
		  "NekMatrix::operator", "Invalid access to _data via parenthesis operator");

	int k1,k2;
	
	switch(_matform){
	case RectangleMatrix: case General_Full:
	  return _mat[i][j]; 
	  break;
	  
	case ZeroMatrix:
	  _data.reset(new DataType[_totalsize]);
	  for(k1=0;k1<_totalsize;++k1)
	    _data[k1] = (DataType)0;

	  _mat.reset(new DataType*[_rowsize]);
	  for (k1=0; k1<_rowsize; ++k1) 
	    _mat[k1] = _data.get() + k1*_columnsize;
	  
	  return _mat[i][j];
	  break;
	  
	case DiagMatrix:
	  if(i==j)
	    return _data[i];
	  else{
	    DataType * tempdata = new DataType[_rowsize];
	    for(k1=0;k1<_rowsize;++k1)
	      tempdata[k1] = _data[k1];

	    _data.reset(new DataType[_totalsize]);
	    for(k1=0;k1<_totalsize;++k1)
	      _data[k1] = (DataType)0;

	    _mat.reset(new DataType*[_rowsize]);
	    for (k1=0; k1<_rowsize; ++k1){
	      _mat[k1] = _data.get() + k1*_columnsize;
	      _mat[k1][k1] = tempdata[k1];
	    }
	      
	    delete[] tempdata;
	    return _mat[i][j];
	  }
	  break;
  
	case Symmetric:
	case Symmetric_Positive:
	case Symmetric_Positive_Banded:

	  if(j>=i)
	    return _data[ (_rowsize*(_rowsize+1) - ((_rowsize-i)*(_rowsize-i+1)))/2 + (j-i)];
	  else{
	    DataType * tempdata = new DataType[_rowsize*(_rowsize+1)/2];
	    for(k1=0;k1<_rowsize*(_rowsize+1)/2;++k1)
	      tempdata[k1] = _data[k1];

	    _data.reset(new DataType[_totalsize]);
	    for(k1=0;k1<_totalsize;++k1)
	      _data[k1] = (DataType)0;

	    _mat.reset(new DataType*[_rowsize]);
	    for (k1=0; k1<_rowsize; ++k1){
	      _mat[k1] = _data.get() + k1*_columnsize;
	      for(k2=0;k2<_columnsize;++k2){
		int itemp;
		if(k2>=k2)
		  itemp =  ((_rowsize*(_rowsize+1)) - ((_rowsize-k1)*(_rowsize-k1+1)))/2 + (k2-k1);
		else
		  itemp =  ((_rowsize*(_rowsize+1)) - ((_rowsize-k2)*(_rowsize-k2+1)))/2 + (k1-k2);
	  
		_mat[k1][k2] = tempdata[itemp];
	      }
	    }
	    
	    delete[] tempdata;
	    return _mat[i][j];
	  }
	  break;


	default:
	  ASSERTL0(0, "NekMatrix::NekMatrix", "Default in Switch Reached");
	}
      }


      /// Necessary Functions for Matrix Operations
      inline int size() const { return _totalsize; }
      inline int rows() const { return _rowsize; }
      inline int cols() const { return _columnsize; }
      inline MatrixForm matform() const { return _matform; }
      inline boost::shared_array<DataType> getdata() const {return _data;}
      inline DataType * getdataptr() const {return _data.get();}
      inline bool unique() const {return _data.unique();};
      
      /// Extra Operators
      


    };


    /////////////////////////////////////////////////////////////////////////////////////////    

    /// NekBlockVector Class
    template<typename DataType>
    class NekBlockVector {
    private:
      int _size;
      boost::shared_array<NekVector<DataType> > _data;
      
    public:
      
      NekBlockVector(int size): _size(size), _data(new NekVector<DataType>[_size]){	

      }

      
      
      /// Extra Operators
    }; 



    /////////////////////////////////////////////////////////////////////////////////////////    

    /// NekBlockMatrix Class
    template <typename DataType>
    class NekBlockMatrix {
    private:
      int _rowsize, _columnsize, _totalsize;
      boost::shared_array<NekMatrix<DataType> > _data;
      boost::shared_array<NekMatrix<DataType>* > _mat;

    public:

      NekBlockMatrix(int rows): _rowsize(rows), _columnsize(rows), _totalsize(_rowsize*_columnsize), 
				_data(new NekMatrix<DataType>[_totalsize]), _mat(new NekMatrix<DataType>*[_rowsize]){
	
	for (int i=0; i<_rowsize; ++i) 
	  _mat[i] = _data.get() + i*_columnsize;
				  
      }


      NekBlockMatrix(int rows, int columns): _rowsize(rows), _columnsize(colums), _totalsize(_rowsize*_columnsize), 
					     _data(new NekMatrix<DataType>[_totalsize]), _mat(new NekMatrix<DataType>*[_rowsize]){
	
	for (int i=0; i<_rowsize; ++i) 
	  _mat[i] = _data.get() + i*_columnsize;
				  
      }


      /// Necessary Functions for Matrix Operations
      inline int size() const { return _totalsize; }
      inline int rows() const { return _rowsize; }
      inline int cols() const { return _columnsize; }
      
      inline NekMatrix<DataType> operator()(int i, int j) const { 
	ASSERTL1( (i>=0) && (i<_rowsize) && (j>=0) && (j<_columnsize), 
		  "NekBlockMatrix::operator", "Invalid access to _data via parenthesis operator");
	
	return _mat[i][j]; 
      }
      
      
      inline NekMatrix<DataType>& operator()(int i, int j){ 
	ASSERTL1( (i>=0) && (i<_rowsize) && (j>=0) && (j<_columnsize), 
		  "NekBlockMatrix::operator", "Invalid access to _data via parenthesis operator");
	
	return _mat[i][j]; 
      }

      inline boost::shared_array<DataType> getdata() const {return _data;}
      inline DataType * getdataptr() const {return _data.get();}
      inline bool unique() const {return _data.unique();};
    };

    /////////////////////////////////////////////////////////////////////////////////////////    

    
    template<typename DataType, int dim>
    NekPoint<DataType,dim> operator*(const DataType &s, const NekPoint<DataType,dim>& P){
      NekPoint<DataType,dim> temp(P);
      for(int i=0;i<dim;++i)
	temp(i) *= s;
      return temp;
    }

    template<typename DataType, int dim>
    NekPoint<DataType,dim> operator*(const NekPoint<DataType,dim>& P, const DataType &s){
      NekPoint<DataType,dim> temp(P);
      for(int i=0;i<dim;++i)
	temp(i) *= s;
      return temp;
    }
    
    template<typename DataType, int dim>
    NekPoint<DataType,dim> operator+(const NekPoint<DataType,dim>& P, const NekVector<DataType>& V){
      ASSERTL0(P.dimension() == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
	
      NekPoint<DataType,dim> temp(P);
      for(int i=0;i<P.dimension();++i)
      	temp(i) = temp(i) + V(i); 
      return temp;
    }
    
    template<typename DataType, int dim>
    NekPoint<DataType,dim> operator+(const NekVector<DataType>& V, const NekPoint<DataType,dim>& P){
      ASSERTL0(dim == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
      
      NekPoint<DataType,dim> temp(P);
      for(int i=0;i<dim;++i)
	temp(i) += V(i);
      return temp;	
    }
    
    template<typename DataType, int dim>

    NekPoint<DataType,dim> operator-(const NekPoint<DataType,dim>& P, const NekVector<DataType>& V){
      ASSERTL0(dim == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
      
      NekPoint<DataType,dim> temp(P);
      for(int i=0;i<dim;++i)
	temp(i) -= V(i);
      return temp;      
    }
  
  }//end namespace
}//end namespace

#endif






