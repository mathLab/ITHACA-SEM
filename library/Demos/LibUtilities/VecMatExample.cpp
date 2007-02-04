#include <cstdio>
#include <cstdlib>

const int size = 10;

int main(int argc, char *argv[])
{
    return 0;
    
//   int i,j;
//   double dx = 2.0*M_PI/(double)(size-1);
// 
//   // Array of doubles
//   boost::shared_array<double> dblptr1(new double[size]), dblptr2(new double[size]);
//   
//   // NekPoint statically allocated
//   NekPoint<double,size> pt1,pt2;
//   
//   // NekPoint dynamically allocated
//   NekPoint<double,size> * pt3 = new NekPoint<double,size>,
//     *pt4 = new NekPoint<double,size>;
//   
//   // NekVector statically allocated
//   NekVector<double> vec1(size),vec2(size);
//   
//   // NekVector dynamically allocated
//   NekVector<double> * vec3 = new NekVector<double>(size),
//     *vec4 = new NekVector<double>(size);
//   
//   // NekMatrix staticallyy allocated
//   NekMatrix<double> mat1(size),mat2(size,size);
//   
//   // NekMatrix dynamically allocated
//   NekMatrix<double> *mat3 = new NekMatrix<double>(size),
//     *mat4 = new NekMatrix<double>(size,size);
// 
//   // Arrays to populate variables
//   for(i=0;i<size;i++){
//       pt1(1) = 1.0;
//       pt2(i) = 1.0 + i*dx; 
//       (*pt3)(i) = sin(i*dx);  
//       (*pt4)(i) = cos(i*dx);  
//       vec1(i) = 1.0;
//       vec2(i) = 1.0 + i*dx; 
//       (*vec3)(i) = sin(i*dx);  
//       (*vec4)(i) = cos(i*dx);  
//     
//     for(j=0;j<size;j++){
//       mat1(i,j) = mat2(i,j) = (*mat3)(i,j) = (*mat4)(i,j) = pow(i*dx,j);
//     }
// 
//   }
//   
// 
//   // Add two points together
//   (*pt4) = pt1 + pt2;
// 
//   // Subtract two points
//   (*pt4) = pt1 - pt2;
//   
//   // Scale a point by a number
//   (*pt3) = 2.0*(*pt4);
//   (*pt4) = (*pt3)*2.0;
// 
//   // Add a vector to a point
//   (*pt3) = pt1 + vec1;
//   (*pt4) = vec2 + pt2;
//   
//   // Subtract a vector from a point
//   pt1 = pt2 - vec2;
//  
//   // Add two vectors together
//   (*vec3) = vec1 + vec2;
// 
//   // Subtract two points
//   (*vec4) = vec1 - vec2;
//   
//   // Scale a point by a number
//   (*vec3) = 2.0*(*vec4);
//   (*vec4) = (*vec3)*2.0;
// 
//   // add mutiple vectors together
//   vec1 = vec1 + vec2 + (*vec3) + (*vec4);
//   
//   // add matrices
//   mat1 = mat2 + (*mat3);
// 
//   // subtract matrices
//   mat2 = mat1 - (*mat4);
  
}
  
  
