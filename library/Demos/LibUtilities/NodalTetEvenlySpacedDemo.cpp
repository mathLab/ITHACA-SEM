
#include <iostream>
#include <iomanip>
#include <iosfwd>

using namespace std;


#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/NodalTetEvenlySpaced.h>

#include "LibUtilities/Foundations/Foundations.hpp"
#include <LibUtilities/Foundations/Points.h>


using namespace Nektar;
using namespace boost;
using namespace Nektar::LibUtilities;

int main(int argc, char *argv[]){ 

   // Argument check: Display a help message if the count is wrong
    if(argc != 3)
    {
        cerr << "Usage: NodalTetEvenlySpacedDemo Points3D-Type nPtsPerSide" << endl;

        cerr << "Where type is an integer value which dictates the basis as:\n";
        for(int i=0; i<SIZE_PointsType; ++i)
        {
            cerr << setw(30) << kPointsTypeStr[i] << " =" << i << endl;
        }

        cerr << "Note type = 1 ~ 14 are for one dimensional basis, and 15, 16, and 17 are for two dimensional basis " << endl;
        
        cerr << "\n Example: NodalTetEvenlySpacedDemo 18 4" << endl;
        
        cerr << "\n\t Tests 3D NodalTetEvenlySpaced on 4 points per side" << endl;

        return 1; // Aborts main() function
    }

    // Read in the type for the points from the caller
    PointsType pointsType = (PointsType) atoi(argv[1]);
    if(pointsType == eNoPointsType)
    {
        cerr << "pointsType = " << pointsType << endl;
        cerr << "kPointsTypeStr[" <<pointsType<< "]=" << kPointsTypeStr[pointsType] << endl;
        ErrorUtil::Error(ErrorUtil::efatal,__FILE__, __LINE__,
                         "No Points Type requested",0);
    }

    // Read in the number of points per side at which the function will be evaluated
    int nPtsPerSide = atoi(argv[2]);

     // Show the set up to the user
    cout << "Points Type:            " << kPointsTypeStr[pointsType] << endl;
    cout << "Number of points per side:       " << nPtsPerSide << endl;

    // Display the example test function to the user
    if( pointsType == eNodalTetEvenlySpaced)
    {
        cout << "Uniform grid on a Tetrahedron points" << endl;
    }


    int degree = nPtsPerSide-1;

    boost::shared_ptr<Points<NekDouble> > points = PointsManager()[PointsKey(nPtsPerSide, pointsType)];

    NodalTetEvenlySpaced *nTete = dynamic_cast<NodalTetEvenlySpaced*>(points.get());
    Array<OneD, const NekDouble> ax, ay, az;
    nTete->GetPoints(ax,ay, az);
    NekVector<NekDouble> x = ToVector(ax), y = ToVector(ay), z = ToVector(az);

    
    // /////////////////////////////////////////////
    // Test Interpolation
    // 
    
    // Make a uniform grid on the Tetrahedron
    int nGridPtsPerSide =6;
    int nGridPts = nGridPtsPerSide * (nGridPtsPerSide + 1) * (nGridPtsPerSide + 2) / 6;
            
            Array<OneD, NekDouble> axi(nGridPts), ayi(nGridPts), azi(nGridPts);
            for( int i = 0, n = 0; i < nGridPtsPerSide; ++i )
            {
                for( int j = 0; j < nGridPtsPerSide - i; ++j )
                {
                    for(int k = 0; k < nGridPtsPerSide - i - j; ++k, ++n )
                    {
                        axi[n] = -1.0 + 2.0*j / (nGridPtsPerSide - 1);
                        ayi[n] = -1.0 + 2.0*i / (nGridPtsPerSide - 1);
                        azi[n] = -1.0 + 2.0*k / (nGridPtsPerSide - 1);
                    }
                }
            }
            // Done with constructing the uniform grid of points on the Tetrahedron
    
    NekVector<NekDouble> xi = ToVector(axi), yi = ToVector(ayi), zi = ToVector(azi);
    boost::shared_ptr<NekMatrix<NekDouble> > Iptr = nTete->GetI(axi, ayi, azi);
    const NekMatrix<NekDouble> & interpMat = *Iptr;

    NekMatrix<NekDouble> matVnn = GetMonomialVandermonde(x, y, z);
    NekMatrix<NekDouble> matVmn = GetMonomialVandermonde(xi, yi, zi, degree);
    NekMatrix<NekDouble> matVmnTilde = interpMat * matVnn;
    NekMatrix<NekDouble> err(xi.GetRows(), x.GetRows());
    NekMatrix<NekDouble> relativeError(GetSize(xi), GetSize(x));
    long double epsilon = 1e-15;

    for( int i = 0; i < int(xi.GetRows()); ++i )
    {
        for( int j = 0; j < int(x.GetRows()); ++j )
        {
            err(i,j) = LibUtilities::MakeRound(1e16 * fabs(matVmnTilde(i,j) - matVmn(i,j)))/1e16;

                // Compute relative error
            relativeError(i,j) = (matVmn(i,j) - matVmnTilde(i,j))/matVmn(i,j);
            if( fabs(matVmn(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeError(i,j) = matVmn(i,j) - matVmnTilde(i,j);
         }
      }
   }

    cout << "\n\n\n************************* NodalTetEvenlySpaced Interpolation Demo *******************************" << endl;
    cout << "\n Result of NumericInterpolation = \n" << MatrixToString(matVmnTilde) << endl;
    cout << "\n Result of Exact Interpolation = \n" << MatrixToString(matVmn) << endl;
    cout << "\n Result of I matrix = \n" << MatrixToString(interpMat) << endl;
    cout << "\n epsilon = \n" << epsilon << endl;
    cout << "\n relativeError : Interpolation = \n" << MatrixToString(relativeError) << endl;
    cout << "\n Error : abs(NumericInterpolation - exact) = \n" << MatrixToString(err) << endl;
    cout << "---------------------------------- End of Interpolation Demo ------------------------------------" << endl;


    // /////////////////////////////////////////////
    //  X Derivative
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dxptr = points->GetD(xDir);
    const NekMatrix<NekDouble> & xderivativeMat = *Dxptr;

    NekMatrix<NekDouble> matVx = GetTetXDerivativeOfMonomialVandermonde(x, y, z);
    NekMatrix<NekDouble> numericXDerivative = xderivativeMat * matVnn;
    NekMatrix<NekDouble> error(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDx(x.GetRows(), y.GetRows());
    long double epsilonDx = 1e-14;

    for(int i=0; i< int(x.GetRows()); ++i )
    {
        for(int j=0; j< int(y.GetRows()); ++j)
        {
                error(i,j) = LibUtilities::MakeRound(1e15 * fabs(numericXDerivative(i,j) - matVx(i, j)))/1e15;

            // Compute relative error
            relativeErrorDx(i,j) = (matVx(i,j) - numericXDerivative(i,j))/matVx(i,j);
            if( fabs(matVx(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeErrorDx(i,j) = matVx(i,j) - numericXDerivative(i,j);
            }
        }
    }

    cout << "\n\n\n**************** NodalTetEvenlySpaced X Derivative Floating Point Error Precision ***************" << endl;
    cout << "\n Result of NumericXDerivative = \n" << MatrixToString(numericXDerivative) << endl;
    cout << "\n Result of exact XDerivative = \n" << MatrixToString(matVx) << endl;
    cout << "\n Result of D matrix = \n" << MatrixToString(xderivativeMat) << endl;
    cout << "epsilon = \n" << epsilonDx <<endl;
    cout << "\n relativeError : X Derivative = \n" << MatrixToString(relativeErrorDx) << endl;
    cout << "\n Error : abs(exact - NumericXDerivative) = \n" << MatrixToString(error) << endl;
    cout << "\n --------------- End of X Derivative Matrix Demo        --------------------------" << endl;


    // /////////////////////////////////////////////
    //  Y Derivative
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dyptr = points->GetD(yDir);
    const NekMatrix<NekDouble> & yderivativeMat = *Dyptr;

    NekMatrix<NekDouble> matVy = GetTetYDerivativeOfMonomialVandermonde(x, y, z);
    NekMatrix<NekDouble> numericYDerivative = yderivativeMat * matVnn;
    NekMatrix<NekDouble> errorDy(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDy(x.GetRows(), y.GetRows());
    long double epsilonDy = 1e-14;
    
    for(int i=0; i< int(x.GetRows()); ++i )
    {
        for(int j=0; j< int(y.GetRows()); ++j)
        {

                error(i,j) = LibUtilities::MakeRound(1e15 * fabs(numericXDerivative(i,j) - matVx(i, j)))/1e15;

            // Compute relative error
            relativeErrorDy(i,j) = (matVy(i,j) - numericYDerivative(i,j))/matVy(i,j);
            if( fabs(matVy(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeErrorDy(i,j) = matVy(i,j) - numericYDerivative(i,j);
            }
        }
    }

    cout << "\n\n\n*************** NodalTetEvenlySpaced Y Derivative Floating Point Error Precision **************" << endl;
    cout << "\n Result of NumericYDerivative = \n" << MatrixToString(numericYDerivative) << endl;
    cout << "\n Result of exact Y Derivative = \n" << MatrixToString(matVy) << endl;
    cout << "\n Result of D matrix = \n" << MatrixToString(yderivativeMat) << endl;
    cout << "epsilon = \n" << epsilonDy <<endl;
    cout << "\n relativeError : Y Derivative = \n" << MatrixToString(relativeErrorDy) << endl;
    cout << "\n Error : abs(exact - NumericYDerivative) = \n" << MatrixToString(error) << endl;
    cout << "\n --------------- End of Y Derivative Matrix          --------------------------" << endl;


     // /////////////////////////////////////////////
    //  Z Derivative
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dzptr = points->GetD(zDir);
    const NekMatrix<NekDouble> & zderivativeMat = *Dzptr;

    NekMatrix<NekDouble> matVz = GetTetZDerivativeOfMonomialVandermonde(x, y, z);
    NekMatrix<NekDouble> numericZDerivative = zderivativeMat * matVnn;
    NekMatrix<NekDouble> errorDz(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDz(x.GetRows(), y.GetRows());
    long double epsilonDz = 1e-14;
    for(int i=0; i< int(x.GetRows()); ++i )
    {
        for(int j=0; j< int(y.GetRows()); ++j)
        {
                error(i,j) = LibUtilities::MakeRound(1e15 * fabs(numericZDerivative(i,j) - matVz(i, j)))/1e15;

            // Compute relative error
            relativeErrorDz(i,j) = (matVz(i,j) - numericZDerivative(i,j))/matVz(i,j);
            if( fabs(matVz(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeErrorDz(i,j) = matVz(i,j) - numericZDerivative(i,j);
            }
        }
    }

    cout << "\n\n\n*************** NodalTetEvenlySpaced Z Derivative Floating Point Error Precision **************" << endl;
    cout << "\n Result of NumericZDerivative = \n" << MatrixToString(numericZDerivative) << endl;
    cout << "\n Result of exact Z Derivative = \n" << MatrixToString(matVz)<< endl;
    cout << "\n Result of D matrix = \n" << MatrixToString(zderivativeMat) << endl;
    cout << "epsilon = \n" << epsilonDz <<endl;
    cout << "\n relativeError : Z Derivative = \n" << MatrixToString(relativeErrorDz) << endl;
    cout << "\n Error : abs(exact - NumericZDerivative) = \n" << MatrixToString(error) << endl;
    cout << "\n --------------- End of Z Derivative Matrix      --------------------------" << endl;

    

//      // /////////////////////////////////////////////
//     // Test Integral
//     // 
//     const Array<OneD, const NekDouble> &weight = points->GetW();
//     NekVector<NekDouble> integMVandermonde = GetIntegralOfMonomialVandermonde(degree);
//     NekVector<NekDouble> numericIntegral = GetTranspose(matVnn) * ToVector(weight);
//     NekVector<NekDouble> errorIntegral = numericIntegral - integMVandermonde;
//     NekVector<NekDouble> relativeErrorIntegral(GetSize(x));
//     long double epsilonIntegral = 1e-16;
//     
//     for(int i=0; i<int(x.GetRows()); ++i){
//         relativeErrorIntegral(i) = (integMVandermonde(i) -  numericIntegral(i))/integMVandermonde(i);
//         if(fabs(integMVandermonde(i)) < epsilonIntegral){
//             relativeErrorIntegral(i) = (integMVandermonde(i) -  numericIntegral(i));
//         }
//     }
// 
//     cout << "------------------------------- NodalTetElec Integral Test -------------------------------" << endl;
//     cout << "\n Result of NumericIntegral = \n" << numericIntegral << endl;
//     cout << "epsilon = \n" << epsilonIntegral << endl;
//     cout << "\n relativeError : Integral = \n" << relativeErrorIntegral << endl;
//     cout << "\n ErrorIntegral : (NumericIntegral - exact) = \n" << errorIntegral << endl;
//     cout << "\n W = \n" << ToVector(weight) << endl;
//     cout << "------------------------------- End of Integral Test ---------------------------------------" << endl;
 
}
