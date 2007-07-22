#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;


#include <LibUtilities/Foundations/NodalUtil.h>

#include "LibUtilities/Foundations/Foundations.hpp"
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace Nektar;
using namespace boost;
using namespace Nektar::LibUtilities;

int main(int argc, char *argv[]){

   // Argument check: Display a help message if the count is wrong
    if(argc != 3){
        cerr << "Usage: NodalTriFeketeDemo Points2D-Type nPtsPerSide" << endl;

        cerr << "Where type is an integer value which dictates the basis as:\n";
        for(int i=0; i<SIZE_PointsType; ++i){
            cerr << setw(30) << PointsTypeMap[i] << " =" << i << endl;
        }

        cerr << "Note type = 1 ~ 14 are for one dimensional basis, and 15 and 16 are for two dimensional basis " << endl;
        
        cerr << "\n Example: NodalTriFeketeDemo 16 3" << endl;
        
        cerr << "\n\t Tests 2D NodalTriFekete on 3 points per side" << endl;

        return 1; // Aborts main() function
    }

    // Read in the type for the points from the caller
    PointsType pointsType = (PointsType) atoi(argv[1]);
    if(pointsType == eNoPointsType){
        cerr << "pointsType = " << pointsType << endl;
        cerr << "PointsTypeMap[" <<pointsType<< "]=" << PointsTypeMap[pointsType] << endl;
        ErrorUtil::Error(ErrorUtil::efatal,__FILE__, __LINE__,
                         "No Points Type requested",0);
    }

    // Read in the number of points per side at which the function will be evaluated
    int nPtsPerSide = atoi(argv[2]);

     // Show the set up to the user
    cout << "Points Type:            " << PointsTypeMap[pointsType] << endl;
    cout << "Number of points per side:       " << nPtsPerSide << endl;

    // Display the example test function to the user
    if( pointsType == eNodalTriFekete){
        cout << "Uniform grid on a Triangle points" << endl;
    }


    int degree = nPtsPerSide-1;

    boost::shared_ptr<Points<NekDouble> > points = PointsManager()[PointsKey(nPtsPerSide, pointsType)];

    NodalTriFekete *ntf = dynamic_cast<NodalTriFekete*>(points.get());
    ConstArray<OneD, NekDouble> ax, ay;
    ntf->GetPoints(ax,ay);
    NekVector<NekDouble> x = toVector(ax), y = toVector(ay);

    
    // /////////////////////////////////////////////
    // Test Interpolation
    // 
    // Make a uniform grid on the triangle
    int nGridPtsPerSide = 5, nGridPts = nGridPtsPerSide * (nGridPtsPerSide + 1) / 2;
    Array<OneD, NekDouble> axi(nGridPts), ayi(nGridPts);
    for( int i = 0, n = 0; i < nGridPtsPerSide; ++i ) {
        for( int j = 0; j < nGridPtsPerSide - i; ++j, ++n ) {
            axi[n] = -1.0 + 2.0*j / (nGridPtsPerSide - 1);
            ayi[n] = -1.0 + 2.0*i / (nGridPtsPerSide - 1);
        }
    }
    
    NekVector<NekDouble> xi = toVector(axi), yi = toVector(ayi);
    boost::shared_ptr<NekMatrix<NekDouble> > Iptr = ntf->GetI(axi, ayi);
    const NekMatrix<NekDouble> & I = *Iptr;

    NekMatrix<NekDouble> Vnn = getMonomialVandermonde(x, y);
    NekMatrix<NekDouble> Vmn = getMonomialVandermonde(xi, yi, degree);
    NekMatrix<NekDouble> VmnTilde = I * Vnn;
    NekMatrix<NekDouble> E(xi.GetRows(), x.GetRows());
    NekMatrix<NekDouble> relativeError(getSize(xi), getSize(x));
    long double epsilon = 1e-15;

    for( int i = 0; i < int(xi.GetRows()); ++i ) {
        for( int j = 0; j < int(x.GetRows()); ++j ) {
            E(i,j) = LibUtilities::round(1e16 * fabs(VmnTilde(i,j) - Vmn(i,j)))/1e16;

                // Compute relative error
            relativeError(i,j) = (Vmn(i,j) - VmnTilde(i,j))/Vmn(i,j);
            if( fabs(Vmn(i,j)) < numeric_limits<double>::epsilon() ) {
                relativeError(i,j) = Vmn(i,j) - VmnTilde(i,j);
         }
      }
   }
   
    cout << "------------------------------- NodalTriFekete Interpolation Test -------------------------------" << endl;
    cout << "\n Result of NumericInterpolation = \n" << VmnTilde << endl;
    cout << "\n Result of I matrix = \n" << I << endl;
    cout << "\n epsilon = \n" << epsilon << endl;
    cout << "\n relativeError : Interpolation = \n" << relativeError << endl;
    cout << "\n Error : abs(NumericInterpolation - exact) = \n" << E << endl;
    cout << "---------------------------------- End of Interpolation Test ------------------------------------" << endl;


    // /////////////////////////////////////////////
    // Test X Derivative 
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dptr = points->GetD();
    const NekMatrix<NekDouble> & D = *Dptr;

    NekMatrix<NekDouble> Vx = getXDerivativeOfMonomialVandermonde(x, y);
    NekMatrix<NekDouble> NumericXDerivative = D * Vnn;
    NekMatrix<NekDouble> Error(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDx(x.GetRows(), y.GetRows());
    long double epsilonDx = 1e-14;
    for(int i=0; i< int(x.GetRows()); ++i ){
        for(int j=0; j< int(y.GetRows()); ++j){

                Error(i,j) = LibUtilities::round(1e15 * fabs(NumericXDerivative(i,j) - Vx(i, j)))/1e15;

            // Compute relative error
            relativeErrorDx(i,j) = (Vx(i,j) - NumericXDerivative(i,j))/Vx(i,j);
            if( fabs(Vx(i,j)) < numeric_limits<double>::epsilon() ) {
                relativeErrorDx(i,j) = Vx(i,j) - NumericXDerivative(i,j);
            }
        }
    }

    cout << "------------------ NodalTriFekete X Derivative Floating Point Error Precision ---------------------------" << endl;
    cout << "\n Result of NumericXDerivative = \n" << NumericXDerivative << endl;
    cout << "\n Result of D matrix = \n" << D << endl;
    cout << "epsilon = \n" << epsilonDx <<endl;
    cout << "\n relativeError : X Derivative = \n" << relativeErrorDx << endl;
    cout << "\n Error : abs(exact - NumericXDerivative) = \n" << Error << endl;
    cout << "\n --------------- End of Testing X Derivative Matrix                          --------------------------" << endl;


     // /////////////////////////////////////////////
    // Test Integral
    // 
    const ConstArray<OneD,NekDouble> &W = points->GetW();
    NekVector<NekDouble> Vwi = getIntegralOfMonomialVandermonde(degree);
    NekVector<NekDouble> NumericIntegral = getTranspose(Vnn) * toVector(W);
    NekVector<NekDouble> ErrorIntegral = NumericIntegral - Vwi;
    NekVector<NekDouble> relativeErrorIntegral(getSize(x));
    long double epsilonIntegral = 1e-16;
    
    for(int i=0; i<int(x.GetRows()); ++i){
        relativeErrorIntegral(i) = (Vwi(i) -  NumericIntegral(i))/Vwi(i);
        if(fabs(Vwi(i)) < epsilonIntegral){
            relativeErrorIntegral(i) = (Vwi(i) -  NumericIntegral(i));
        }
    }

    cout << "------------------------------- NodalTriFekete Integral Test -------------------------------" << endl;
    cout << "\n Result of NumericIntegral = \n" << NumericIntegral << endl;
    cout << "epsilon = \n" << epsilonIntegral << endl;
    cout << "\n relativeError : Integral = \n" << relativeErrorIntegral << endl;
    cout << "\n ErrorIntegral : (NumericIntegral - exact) = \n" << ErrorIntegral << endl;
    cout << "\n W = \n" << toVector(W) << endl;
    cout << "------------------------------- End of Integral Test ---------------------------------------" << endl;
 
}
