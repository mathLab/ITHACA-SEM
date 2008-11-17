///////////////////////////////////////////////////////////////////////////////
//
// File StreamFunctionWaves.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Class for computing semi-analytical solutions of the stream 
//              function waves theory. A good tool for verification of wave 
//              propagation codes based on potential theory
//
///////////////////////////////////////////////////////////////////////////////

#include <BoussinesqSolver/StreamFunctionWaves.h>
#include <BoussinesqSolver/BoussinesqEquations.h>

namespace Nektar
{
  StreamFunctionWaves::StreamFunctionWaves():
    m_init(),
    m_f(),
    m_f1(),
    m_f2()
  {
  }

  
  void StreamFunctionWaves::SetUpParameters(NekDouble Time, NekDouble Percent, NekDouble WaveLength, NekDouble SWL, NekDouble g, int nCoeffs, int nSteps)
  {
    m_time       = Time;
    m_d          = SWL;
    m_wavelength = WaveLength;
    m_wavenumber = (2.0*M_PI)/WaveLength;
    m_kd         = m_wavenumber*m_d;
    m_coeffs     = nCoeffs;
    m_steps      = nSteps;
    m_g          = g;


    // if the non-dimensional wavenumber larger than 3
    // then outside the range of validity of weakly 
    // dispersive Boussinesq models

    cout << scientific << "Dimensionless depth: " << m_kd << endl;
    cout << "Wave Steepness: " << Percent << " %" << endl;


    if (m_kd > 3.0)
      {
	cout << "Dispersion outside the range of validity: kd = "<<m_kd << endl;
	exit(1);
      }
      
    m_waveheight = (Percent/100.0)*m_wavelength*(0.141063+0.0095721*(2.0*M_PI/m_kd)
					       +0.0077829*pow(2.0*M_PI/m_kd,2.0))/
      (1.0+0.0788340*(2.0*M_PI/m_kd)+0.0317567*pow(2.0*M_PI/m_kd,2.0)+0.0093407*pow(2.0*M_PI/m_kd,3.0));
    
    
    // initialize arrays
    m_coefficientsA = Array<OneD, NekDouble>(m_coeffs,0.0);
    m_coefficientsB = Array<OneD,NekDouble>(m_coeffs,0.0);
    m_init = Array<OneD, NekDouble>(2*m_coeffs+5,0.0);
    m_f    = Array<OneD, NekDouble>(2*m_coeffs+5,0.0);
    m_f1   = Array<OneD, NekDouble>(2*m_coeffs+5,0.0);
    m_f2   = Array<OneD, NekDouble>(2*m_coeffs+5,0.0);
  
    // set init
    NekDouble WavePeriod = (2.0*M_PI)/sqrt(m_g*m_wavenumber*tanh(m_kd));
    m_init[0] = m_wavelength/WavePeriod;
    m_init[1] = m_wavelength/WavePeriod;
    m_init[2] = (m_g/(2.0*m_wavenumber))*tanh(m_kd);
    m_init[3] = 0.0;
    
    m_hi = m_waveheight/m_steps;
    //std::cout << "Hi = " << m_hi << std::endl;
    for (int i = 0; i < m_coeffs+1; ++i)
      {
	m_init[i+4] = 0.5*m_hi*cos((i+1)*(m_wavelength/(2.0*m_coeffs)));
      }
       
    m_init[5+m_coeffs] = m_g*0.5*m_hi*WavePeriod/m_wavelength;
 
  }
  
  
  void StreamFunctionWaves::StreamFunctionSolve(NekDouble tol, int nMaxIterations)
  {
  
    for (int i = 1; i < m_steps+1; ++i)
      {

	if (i == 1)
	  {
	    cout << "SF Waves - step: "<< i<< endl;
	    NonLinearSolve(m_init,m_f1,nMaxIterations,tol);
	    NonLinearSolve(m_f1,m_f2,nMaxIterations,tol);
	  }
 	else
 	  {
	    cout << "SF Waves - step: "<< i<< endl;
	    m_hi = i*(m_waveheight/m_steps);
	    ExtrapolateSolution();
	    NonLinearSolve(m_init,m_f,nMaxIterations,tol);
	    UpdateSolution();
	  }
      }
    
  // Extract the data
  m_c    = m_f2[0];
  m_ubar = m_f2[1];
  //R    = m_f2[2];
  //Q    = m_f2[3];
    
    Array<OneD, NekDouble> eta(m_coeffs+1);

  for (int i = 0; i < m_coeffs+1; ++i)
    {
      eta[i] = m_f2[i+4];
      // cout << "eta["<<i<<"] = " << eta[i] << endl;

    }
 
  // cout <<"wlength" << m_wavelength <<" wnumb = " << m_wavenumber << endl;
 
  for (int i = 0; i < m_coeffs; ++i)
    {
      m_coefficientsA[i] = (2.0/m_coeffs)*(0.5*(eta[0]+eta[m_coeffs]*cos(m_wavenumber*(i+1)*0.5*m_wavelength)));
      m_coefficientsB[i] = m_f2[m_coeffs+5+i];
      
      for (int j = 1; j < m_coeffs; ++j)
	{
	  m_coefficientsA[i] = m_coefficientsA[i] + (2.0/m_coeffs)*eta[j]*cos((j)*m_wavenumber*(i+1)*m_wavelength/(2.0*m_coeffs));
	}
    }
 
 
  for (int j = 0; j < m_coeffs; ++j)
    {
      cout << "Coefficient for A["<<j<<"] = " << m_coefficientsA[j] << endl;
    }
  
  for (int j = 0; j < m_coeffs; ++j)
    {
      cout << "Coefficient for B["<<j<<"] = " << m_coefficientsB[j] << endl;
    }
  
  
}



  
  void StreamFunctionWaves::EvaluateFunction(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &f)
  {
    
    //----------------------------------------
    // extract values from arguments
    
    NekDouble WavePeriod = (2.0*M_PI)/sqrt(m_g*m_wavenumber*tanh(m_kd));
    NekDouble c    = x[0];
    NekDouble ubar = x[1];
    NekDouble R    = x[2];
    NekDouble Q    = x[3];
    

    Array<OneD, NekDouble> eta(m_coeffs+1);
    Array<OneD, NekDouble> B(m_coeffs);

    for (int i = 0; i < m_coeffs+1; ++i)
      {
	eta[i] = x[4+i];
      }

    for (int i = 0; i < m_coeffs; ++i)
      {
	B[i]   = x[5+m_coeffs+i];
      }
   
    //----------------------------------------

    //----------------------------------------
    // compute temporary vectors and arrays
    
    Array<OneD, NekDouble> tmpCoshkd(m_coeffs);


    for (int i = 0; i < m_coeffs; ++i)
      {
 	tmpCoshkd[i] = cosh((i+1)*m_kd);
      }
    
    
    Array<TwoD, NekDouble> tmp1(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmp2(m_coeffs,m_coeffs+1);
    for (int i = 0; i < m_coeffs; ++i)
      for (int j = 0; j < m_coeffs+1; ++j)
	{
	  tmp1[i][j] = (i+1)*m_wavenumber*(eta[j]+m_d);
	  tmp2[i][j] = (i+1)*j*(M_PI/m_coeffs);
	}
    
    
    Array<TwoD, NekDouble> tmpCosh(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmpSinh(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmpCos(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmpSin(m_coeffs,m_coeffs+1);
    for (int i = 0; i < m_coeffs; ++i)
      for (int j = 0; j < m_coeffs+1; ++j)
	{
	  tmpCosh[i][j] = cosh(tmp1[i][j]);
	  tmpSinh[i][j] = sinh(tmp1[i][j]);
	  tmpCos[i][j]  = cos(tmp2[i][j]);
	  tmpSin[i][j]  = sin(tmp2[i][j]);
	}

    
    Array<TwoD, NekDouble> tmpSinhCos(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmpCoshCos(m_coeffs,m_coeffs+1);
    Array<TwoD, NekDouble> tmpSinhSin(m_coeffs,m_coeffs+1);
    
    for (int i = 0; i < m_coeffs; ++i)
      for (int j = 0; j < m_coeffs+1; ++j)
	{
	  tmpSinhCos[i][j] = tmpSinh[i][j]*tmpCos[i][j];
	  tmpCoshCos[i][j] = tmpCosh[i][j]*tmpCos[i][j];
	  tmpSinhSin[i][j] = tmpSinh[i][j]*tmpSin[i][j];
	}
    
    //------------------------------------------------------


    //-----------------------------------------------
    // compute estimates
    
    Array<OneD, NekDouble> tmpB1(m_coeffs);
    Array<OneD, NekDouble> tmpB2(m_coeffs);

    for (int i = 0; i < m_coeffs; ++i)
      {
	tmpB1[i] = B[i]/tmpCoshkd[i];
	tmpB2[i] = tmpB1[i]/((i+1)*m_wavenumber);
      }
	

    Array<OneD, NekDouble> psi(m_coeffs+1,0.0);
    Array<OneD, NekDouble> u(m_coeffs+1,-ubar);
    Array<OneD, NekDouble> w(m_coeffs+1,0.0);

    Vmath::Smul(m_coeffs+1,-ubar,eta,1,psi,1);
      
    for (int i = 0; i < m_coeffs+1; ++i)
      for (int j = 0; j < m_coeffs; ++j)
	{
	  psi[i] += tmpB2[j]*tmpSinhCos[j][i];
	  u[i]   += tmpB1[j]*tmpCoshCos[j][i];
	  w[i]   += tmpB1[j]*tmpSinhSin[j][i];
	}
    
    
    // fill function

    f[0] = eta[0] - eta[m_coeffs] - m_hi;
    f[1] = eta[0] + eta[m_coeffs] + 2.0 * Vmath::Vsum(m_coeffs-1,eta+1,1);
    f[2] = ubar - c;

    for (int i = 0; i < m_coeffs+1; ++i)
      {
	f[3+i]          = psi[i] - Q;
	f[4+m_coeffs+i] = m_g * eta[i] + 0.5*(u[i]*u[i]+w[i]*w[i]) - R;
      } 
    
  }


  // using forward difference formulae
  void StreamFunctionWaves::EvaluateJacobian(Array<OneD, NekDouble> &x,
					     Array<OneD, NekDouble> &f,
					     Array<TwoD, NekDouble> &jac)
  {
    
    NekDouble tmp, h;
    
    Array<OneD, NekDouble> xh(2*m_coeffs+5);
    Vmath::Smul(2*m_coeffs+5,1.0,x,1,xh,1);
    
    Array<OneD, NekDouble> fh(2*m_coeffs+5,0.0);

    for (int j = 0; j < 2*m_coeffs+5; ++j)
      {

	h = 0.001 * xh[j];
	
	if (fabs(xh[j]) < 1.0e-4)
	  h = 1.0e-5;

	xh[j] += h;
 	
	EvaluateFunction(xh,fh);

	xh[j] -= h; 
	
	for (int i = 0; i < 2*m_coeffs+5; ++i)
	  {
	    jac[i][j] = (fh[i] - f[i])/h;
	    // cout << "h = "<<h<<"   jac["<<i<<"]["<<j<<"] = " << jac[i][j] << endl;
	  }
      }
  }


  void StreamFunctionWaves::NonLinearSolve(Array<OneD, NekDouble> &xIn, Array<OneD, NekDouble> & xOut, int nTrial, NekDouble tol)
  {

    int matrixRows = 2*m_coeffs+5;
    int matrixCols = 2*m_coeffs+5;
    
    Array<TwoD, NekDouble>jac(matrixRows,matrixCols,0.0);
    Array<OneD, NekDouble>f(matrixRows,0.0);
    
    double *p         = (double *)malloc(matrixRows*sizeof(double));
    double *lapackJac = (double *)malloc((matrixRows*matrixCols)*sizeof(double));
    int    *ipivot    = (int *)malloc(matrixRows*sizeof(int));
    
    NekDouble errorf, errorx;
    NekDouble tolf =tol;
    NekDouble tolx =tol;
        
    int info=0;

    Vmath::Vcopy(matrixRows,xIn,1,xOut,1);

    for (int i = 0; i < nTrial; ++i)
      {
 	EvaluateFunction(xOut,f);
 	EvaluateJacobian(xOut,f,jac);

	errorf = 0.0;
	
	for (int j = 0; j < matrixRows; ++j)
	  {
	    errorf += fabs(f[j]);
	  }
	
 	if (errorf <= tolf)
 	  return;
	
 	for (int j = 0; j < matrixRows; ++j)
 	  {
 	    p[j] = -f[j];
 	  }

	//------------------------------------
 	// lapack solve goes here
 	for (int j = 0; j < matrixRows; ++j)
 	  for (int k = 0; k < matrixCols; ++k)
 	    lapackJac[j*matrixRows+k] = jac[j][k];
	
	Lapack::Dgetrf(matrixRows, matrixCols, lapackJac, matrixRows, ipivot, info);
	
	if (info > 0)
	  {
	    cout<<"dgetrf failed:"<< info << endl;
	    exit(1);
	  }
	else if (info < 0)
	  {
	    cout<<"dgetrf failed:"<< info <<endl;
	    exit(1);
	  }
	
 
	int nrhs = 1; // Only 1 right hand side.
	// note: use transformed factorized matrix since row-major in C++
	Lapack::Dgetrs('T', matrixRows, nrhs, lapackJac, matrixRows, ipivot, p, matrixRows, info);
	

	if (info > 0)
	  {
	    cout<<"dgetrs failed"<< endl;
	    exit(1);
	  }
	else if (info < 0)
	  {
	    cout<<"dgetrs failed"<< endl;
	    exit(1);
	  }

	//------------------------------------

	errorx = 0.0;
	for (int j = 0; j < matrixRows; ++j)
	  {
	    errorx += fabs(p[j]);
	    xOut[j] += p[j];
 	  }
	if (errorx <= tolx)
	  return;
	    
      }
    
    free(p);
    free(lapackJac);
    free(ipivot);
    
  }
    

  void StreamFunctionWaves::ExtrapolateSolution(void)
  {
    for (int i = 0; i < 2*m_coeffs+5; ++i)
      {
	m_init[i] = 2.0 * m_f2[i] - m_f1[i];
      }
  }

  void StreamFunctionWaves::UpdateSolution(void)
  {
    for (int i = 0; i < 2*m_coeffs+5; ++i)
      {
	m_f1[i] = m_f2[i];
	m_f2[i] = m_f[i];
      }
  }


  void StreamFunctionWaves::SetConsVariables(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &h,
					     Array<OneD, NekDouble> &hu, Array<OneD, NekDouble> &hv)
  {
    
    int nTotQuadPoints = x.num_elements();
    NekDouble d = -m_d; // z = -d for integration
    
    // move to moving frame 
    Vmath::Sadd(nTotQuadPoints,-m_c*m_time,x,1,x,1);
    
    // zero all variables
    Vmath::Fill(nTotQuadPoints,0.0,h,1);
    Vmath::Fill(nTotQuadPoints,0.0,hu,1);
    Vmath::Fill(nTotQuadPoints,0.0,hv,1);

    
    // compute eta 
    for (int i = 0; i < nTotQuadPoints; ++i)
      for (int j = 0; j < m_coeffs; ++j)
	{
	  h[i] += cos(m_wavenumber*x[i]*(j+1)) * m_coefficientsA[j];
	}

   
    // get gll points and weights for the u levels 
    int nLevelPoints = 40;
    Array<OneD, NekDouble> quadZeros;
    Array<OneD, NekDouble> quadWeights;
    LibUtilities::PointsType quadPointsType = LibUtilities::eGaussGaussLegendre;
    const LibUtilities::PointsKey quadPointsKey(nLevelPoints, quadPointsType);
    quadZeros   = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    quadWeights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

    NekDouble u, z, jacobian;
    	
    for (int i = 0; i< nTotQuadPoints; ++i)
      {
	// Analytically evaluate the Jacobian (note h contains eta).
	jacobian = 0.5 * (h[i] - d);
			
	// do the numerical integration
	for(int j = 0; j < nLevelPoints; j++)
	  {
	    // Calculate the local coordinate of quadrature zeros
	    z = d * (1-quadZeros[j])/2.0 + h[i] * (1+quadZeros[j])/2.0;
	    //   cout <<"z["<<j<<"] = " << z << endl;
	    
	    // Calculate the values of u(z)
	    u = m_c - m_ubar;
	    for (int k = 0; k < m_coeffs; ++k)
	      {
		u += (k+1)*m_coefficientsB[k]*cos((k+1)*m_wavenumber*x[i])*
		  ((cosh((k+1)*m_wavenumber*(m_d+z)))/(cosh((k+1)*m_wavenumber*m_d)));
	      }
	    // integration 
	    hu[i] += quadWeights[j] * u * jacobian;
	  }
      }
    
    
    // compute the total water depth
    Vmath::Sadd(nTotQuadPoints,m_d,h,1,h,1);
    

  }
  



} //end of namespace
