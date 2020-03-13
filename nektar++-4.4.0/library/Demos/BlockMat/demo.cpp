#include "BlockMat.h"

/*
  g++  -I../../../include -g -o demo demo.cpp -L../../ -lBlockMat -lblas -lg2c -llapack
*/

using namespace blockmat;

main(){
  int i;
  BlockMat *A,*B,*C;
  double *mat;

  /* declare matrix
      | 1   2  -   - | 
  A = | 3   4  -   - |   using submatrix mat = | 1  2 |
      | -   -  1   2 |                         | 3  4 |
      | -   -  3   4 | 
  */

  mat = new double [4];
  mat[0] = 1; mat[1] = 2; mat[2] = 3; mat[3] = 4;
  
  cout << "A: " << endl;
  A = new BlockMat(2,2);
  A->GenBlock(0,0,2,2,mat);
  A->GenBlock(0,1,2,2,mat);
  A->GenBlock(1,1,2,2,mat);
  A->PrintBlocks();  

  cout << endl << "B: " << endl;
  B = new BlockMat(2,2);
  B->GenBlock(0,0,2,2,mat);
  B->GenBlock(1,1,2,2,mat);
  B->PrintBlocks();  

  cout << endl << "C=A+B: " << endl;
  C = new BlockMat(2,2);
  C->add(*A,*B);
  C->PrintBlocks();


  cout << endl << "C=A-B: " << endl;
  C->sub(*A,*B);
  C->PrintBlocks();

  cout << endl << "C=A*B: " << endl;
  C->MxM(*A,*B);
  C->PrintBlocks();


  double *y = new(double)[4];
  double *v = new(double)[4];
  Vmath::fill(4,1.0,v,1);
  Vmath::zero(4,y,1);

  cout << endl << "y = A*v: "<< endl;
  C->Mxvpy(v,y);
  for(i = 0; i < 4; ++i)
    cout << y[i] << " ";
  cout << endl;

  Vmath::zero(4,y,1);
  cout << endl << "y = A^T*v: "<< endl;
  //C->Mtxvpy(v,y);
  C->geMxv(ColMajor,1,v,1,y);
  for(i = 0; i < 4; ++i)
    cout << y[i] << " ";
  cout << endl;

  delete A;
  delete B;
  delete C;
  delete[] y;
  delete[] v;
  delete[] mat;

  return  0;
}
