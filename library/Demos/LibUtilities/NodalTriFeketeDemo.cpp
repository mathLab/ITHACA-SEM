
#include <iostream>
#include <iomanip>

using namespace std;


#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>

using namespace boost;
using namespace Nektar;
using namespace Nektar::LibUtilities;

int main(int argc, char *argv[])
{

   // Argument check: Display a help message if the count is wrong
    if(argc != 3)
    {
        cerr << "Usage: NodalTriFeketeDemo Points2D-Type nPtsPerSide" << endl;

        cerr << "Where type is an integer value which dictates the basis as:" << endl;
        for(int i=0; i<SIZE_PointsType; ++i)
        {
            cerr << setw(30) << kPointsTypeStr[i] << " =" << i << endl;
            
        }

        cerr << "Note type = 1 ~ 14 are for one dimensional basis, and 15 and 16 are for two dimensional basis " << endl;
        
        cerr << "\n Example: NodalTriFeketeDemo 16 3" << endl;
        
        cerr << "\n\t Tests 2D NodalTriFekete on 3 points per side" << endl;

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
    if( pointsType == eNodalTriFekete)
    {
        cout << "Uniform grid on a Triangle points" << endl;
    }
    cout << "Fekete!!!" << endl;

    int degree = nPtsPerSide-1;

    boost::shared_ptr<Points<NekDouble> > points = PointsManager()[PointsKey(nPtsPerSide, pointsType)];

    NodalTriFekete *ntf = dynamic_cast<NodalTriFekete*>(points.get());
    Array<OneD, const NekDouble> ax, ay;
    ntf->GetPoints(ax,ay);
    NekVector<NekDouble> x = ToVector(ax), y = ToVector(ay);

    
    // /////////////////////////////////////////////
    // Interpolation
    // 
    // Make a uniform grid on the triangle
    int nGridPtsPerSide = 5, nGridPts = nGridPtsPerSide * (nGridPtsPerSide + 1) / 2;
    Array<OneD, NekDouble> axi(nGridPts), ayi(nGridPts);
    for( int i = 0, n = 0; i < nGridPtsPerSide; ++i )
    {
        for( int j = 0; j < nGridPtsPerSide - i; ++j, ++n )
        {
            axi[n] = -1.0 + 2.0*j / (nGridPtsPerSide - 1);
            ayi[n] = -1.0 + 2.0*i / (nGridPtsPerSide - 1);
        }
    }
    
    NekVector<NekDouble> xi = ToVector(axi), yi = ToVector(ayi);
    boost::shared_ptr<NekMatrix<NekDouble> > Iptr = ntf->GetI(axi, ayi);
    const NekMatrix<NekDouble> & interpMat = *Iptr;

    NekMatrix<NekDouble> matVnn = GetMonomialVandermonde(x, y);
    NekMatrix<NekDouble> matVmn = GetMonomialVandermonde(xi, yi, degree);
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
 
    cout << "\n\n\n***********************   NodalTriFekete Interpolation Demo ***********************" << endl;
    cout << "\n Result of NumericInterpolation = \n" << MatrixToString(matVmnTilde) << endl;
    cout << "\n Result of I matrix = \n" << MatrixToString(interpMat) << endl;
    cout << "\n epsilon = \n" << epsilon << endl;
    cout << "\n relativeError : Interpolation = \n" << MatrixToString(relativeError) << endl;
    cout << "\n Error : abs(NumericInterpolation - exact) = \n" << MatrixToString(err) << endl;
    cout << "---------------------------------- End of Interpolation Demo ------------------------------------" << endl;


    // /////////////////////////////////////////////
    // X Derivative
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dxptr = points->GetD(xDir);
    const NekMatrix<NekDouble> & xderivativeMat = *Dxptr;

    NekMatrix<NekDouble> matVx = GetXDerivativeOfMonomialVandermonde(x, y);
    NekMatrix<NekDouble> numericXDerivative = xderivativeMat * matVnn;
    NekMatrix<NekDouble> errorDx(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDx(x.GetRows(), y.GetRows());
    long double epsilonDx = 1e-14;
    for(int i=0; i< int(x.GetRows()); ++i )
    {
        for(int j=0; j< int(y.GetRows()); ++j)
        {

                errorDx(i,j) = LibUtilities::MakeRound(1e15 * fabs(numericXDerivative(i,j) - matVx(i, j)))/1e15;

            // Compute relative error
            relativeErrorDx(i,j) = (matVx(i,j) - numericXDerivative(i,j))/matVx(i,j);
            if( fabs(matVx(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeErrorDx(i,j) = matVx(i,j) - numericXDerivative(i,j);
            }
        }
    }

    cout << "\n\n\n************  NodalTriFekete X Derivative Floating Point Error Precision ***********" << endl;
    cout << "\n Result of NumericXDerivative = \n" << MatrixToString(numericXDerivative) << endl;
    cout << "\n Result of D matrix = \n" << MatrixToString(xderivativeMat) << endl;
    cout << "epsilon = \n" << epsilonDx <<endl;
    cout << "\n relativeError : X Derivative = \n" << MatrixToString(relativeErrorDx) << endl;
    cout << "\n Error : abs(exact - NumericXDerivative) = \n" << MatrixToString(errorDx) << endl;
    cout << "\n --------------- End of X Derivative Matrix             --------------------------" << endl;


  // /////////////////////////////////////////////
    // Y Derivative
    //    
    boost::shared_ptr<NekMatrix<NekDouble> > Dyptr = points->GetD(yDir);
    const NekMatrix<NekDouble> & yderivativeMat = *Dyptr;

    NekMatrix<NekDouble> matVy = GetYDerivativeOfMonomialVandermonde(x, y);
    NekMatrix<NekDouble> numericYDerivative = yderivativeMat * matVnn;
    NekMatrix<NekDouble> errorDy(x.GetRows(), y.GetRows());
    NekMatrix<NekDouble> relativeErrorDy(x.GetRows(), y.GetRows());
    long double epsilonDy = 1e-14;
    for(int i=0; i< int(x.GetRows()); ++i )
    {
        for(int j=0; j< int(y.GetRows()); ++j)
        {

                errorDy(i,j) = LibUtilities::MakeRound(1e15 * fabs(numericYDerivative(i,j) - matVy(i, j)))/1e15;

            // Compute relative error
            relativeErrorDy(i,j) = (matVy(i,j) - numericYDerivative(i,j))/matVy(i,j);
            if( fabs(matVy(i,j)) < numeric_limits<double>::epsilon() )
            {
                relativeErrorDy(i,j) = matVy(i,j) - numericYDerivative(i,j);
            }
        }
    }

    cout << "\n\n\n************  NodalTriFekete Y Derivative Floating Point Error Precision  ***********" << endl;
    cout << "\n Result of NumericYDerivative = \n" << MatrixToString(numericYDerivative) << endl;
    cout << "\n Result of D matrix = \n" << MatrixToString(yderivativeMat) << endl;
    cout << "epsilon = \n" << epsilonDy <<endl;
    cout << "\n relativeError : Y Derivative = \n" << MatrixToString(relativeErrorDy) << endl;
    cout << "\n Error : abs(exact - NumericYDerivative) = \n" << MatrixToString(errorDy) << endl;
    cout << "\n --------------- End of Y Derivative Matrix           --------------------------" << endl;


    // /////////////////////////////////////////////
    //  Integral
    // 
    const Array<OneD, const NekDouble> &weight = points->GetW();
    NekVector<NekDouble> integMVandermonde = GetIntegralOfMonomialVandermonde(degree);
    NekVector<NekDouble> numericIntegral = GetTranspose(matVnn) * ToVector(weight);
    NekVector<NekDouble> errorIntegral = numericIntegral - integMVandermonde;
    NekVector<NekDouble> relativeErrorIntegral(GetSize(x));
    long double epsilonIntegral = 1e-16;
    
    for(int i=0; i<int(x.GetRows()); ++i)
    {
        relativeErrorIntegral(i) = (integMVandermonde(i) -  numericIntegral(i))/integMVandermonde(i);
        if(fabs(integMVandermonde(i)) < epsilonIntegral)
        {
            relativeErrorIntegral(i) = (integMVandermonde(i) -  numericIntegral(i));
        }
    }

    cout << "\n\n\n *******************  NodalTriFekete Integral Demo  ******************************" << endl;    
    cout << "\n Result of NumericIntegral = \n" << VectorToString(numericIntegral) << endl;
    cout << "epsilon = \n" << epsilonIntegral << endl;
    cout << "\n relativeError : Integral = \n" << VectorToString(relativeErrorIntegral) << endl;
    cout << "\n ErrorIntegral : (NumericIntegral - exact) = \n" << VectorToString(errorIntegral) << endl;
    cout << "\n W = \n" << ToVector(weight) << endl;
    cout << "------------------------------- End of Integral Demo ---------------------------------------" << endl;
 
}
