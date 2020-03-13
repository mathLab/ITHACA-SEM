#include "BlockMat.h"

/*
 g++ -I../../../include -g -o SCdemo SCdemo.cpp -L../../ -lBlockMat -lblas -llapack -lg2c
*/

using namespace blockmat;

/* declare matrix
   
      | 1   2  -   - | 
  A = | 3   4  -   - |                      using submatrix mat = | 1  2 |
      | -   -  1   2 |                                            | 3  4 |
      | -   -  3   4 | 

      | 1   2   3   -   -   -   -   -   - | 
  B = | 4   5   6   -   -   -   -   -   - | using submatrix mat = | 1  2  3 |
      | 1   2   3   1   2   3   1   2   3 |                       | 4  5  6 | 
      | 4   5   6   4   5   6   4   5   6 | 

      | 1   2   3   -   -   -   -   -   - | 
      | 4   5   6   -   -   -   -   -   - | 
  C = | 7   0   0   -   -   -   -   -   - | using submatrix mat = | 1  2  3 |
      | -   -   -   1   2   3   -   -   - |                       | 4  5  6 | 
      | -   -   -   4   5   6   -   -   - |                       | 7  0  0 | 
      | -   -   -   7   0   0   -   -   - | 
      | -   -   -   -   -   -   1   2   3 | 
      | -   -   -   -   -   -   4   5   6 | 
      | -   -   -   -   -   -   7   0   0 | 

      | 1   2   1   2 | 
  D = | 3   4   3   4 |                     using submatrix mat = | 1  2 |
      | 5   6   5   6 |                                           | 3  4 | 
      | -   -   1   2 |                                           | 5  6 | 
      | -   -   3   4 | 
      | -   -   5   6 | 
      | -   -   1   2 |      
      | -   -   3   4 | 
      | -   -   5   6 | 

      Calculate DC = A - B*C*D
  */

main(){
  BlockMat *A,*B,*C,*D,*T,*SC;
  double *mat;

  mat = new double [9];
  mat[0] = 1; mat[1] = 2; mat[2] = 3; 
  mat[3] = 4; mat[4] = 5; mat[5] = 6;
  mat[6] = 7; mat[7] = 0; mat[8] = 0; 
  
  cout << "A: " << endl;
  A = new BlockMat(2,2);
  A->GenBlock(0,0,2,2,mat);
  A->GenBlock(1,1,2,2,mat);
  A->PrintBlocks();  

  cout << endl << "B: " << endl;
  B = new BlockMat(2,3);
  B->GenBlock(0,0,2,3,mat);
  B->GenBlock(1,0,2,3,mat);
  B->GenBlock(1,1,2,3,mat);
  B->GenBlock(1,2,2,3,mat);
  B->PrintBlocks();  

  cout << endl << "C: " << endl;
  C = new BlockMat(3,3);
  C->GenBlock(0,0,3,3,mat);
  C->GenBlock(1,1,3,3,mat);
  C->GenBlock(2,2,3,3,mat);
  C->PrintBlocks();
  
  cout << endl << "C^{-1}: " << endl;
  C->invert_diag();
  C->PrintBlocks();

  cout << endl << "D: " << endl;
  D = new BlockMat(3,2);
  D->GenBlock(0,0,3,2,mat);
  D->GenBlock(0,1,3,2,mat);
  D->GenBlock(1,1,3,2,mat);
  D->GenBlock(2,1,3,2,mat);
  D->PrintBlocks();

  cout << endl << "SC=A-B*C*D: " << endl;
  SC = new BlockMat(2,2);
  T  = new BlockMat(3,2);
  // T->geMxM(RowMajor,RowMajor,1,*C,*D,0);  
  T->MxM(*C,*D);  
  SC->sub(*A,SC->MxM(*B,*T));
  SC->PrintBlocks();

  cout << endl << "SC=A-D^T*C*D: " << endl;
  T->MxM(*C,*D); 
  SC->sub(*A,SC->MtxM(*D,*T));
  SC->PrintBlocks();

  cout << endl << "SC=A-B*C*B^T: " << endl;
  T->MxMt(*C,*B); 
  //T->geMxM(RowMajor,ColMajor,1,*C,*B,0); 
  SC->sub(*A,SC->MxM(*B,*T));
  SC->PrintBlocks();

  delete A;
  delete B;
  delete C;
  delete D;
  delete T;
  delete SC;
  delete[] mat;

  return 0;
}
