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
// Description: Class for computing Stream Function Waves to be used for 
//              validation of waves codes based on potential theory
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_STREAMFUNCTIONWAVES_H
#define NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_STREAMFUNCTIONWAVES_H

#include <BoussinesqSolver/BoussinesqEquations.h>

namespace Nektar
{
    /**
     * \brief This class provides semi-analytical solution
     * for the stream function wave theory (Fenton ,1982). 
     **/
  
  class StreamFunctionWaves

  {
  public:
    /**
     * \brief constructor
     */ 
    StreamFunctionWaves(void);

    void SetUpParameters(NekDouble Time, NekDouble Percent, NekDouble WaveLenght, NekDouble SWL, NekDouble g, int nCoeffs, int nSteps);

    void StreamFunctionSolve(NekDouble tol, int nMaxIterations);
  
    void SetConsVariables(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &h,
			  Array<OneD, NekDouble> &hu, Array<OneD, NekDouble> &hv);
    
      inline void SetWaveLength(NekDouble WaveLength)
    {
      m_wavelength = WaveLength;
    }

    inline void Setkd(NekDouble kd)
    {
      m_kd = kd;
      cout << "kd = " << m_kd << endl;
    }
    
  protected:
    
    
  private:
    void ComputeFirstSolutionStep(void);

    void EvaluateFunction(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &f);
  
    void EvaluateJacobian(Array<OneD, NekDouble> &x, Array<OneD, NekDouble> &f, Array<TwoD, NekDouble> &jac);

    void NonLinearSolve(Array<OneD, NekDouble> &xIn, Array<OneD, NekDouble> &xOut, int nTrial, NekDouble tol);

    void ExtrapolateSolution(void);
      
    void UpdateSolution(void);

    NekDouble m_time;
    NekDouble m_wavelength;
    NekDouble m_wavenumber;
    NekDouble m_d;
    NekDouble m_kd;
    NekDouble m_waveheight;
    int       m_coeffs;
    int       m_steps;
    NekDouble m_g;
    NekDouble m_hi;
    NekDouble m_c;
    NekDouble m_ubar;
    Array<OneD, NekDouble> m_coefficientsA;
    Array<OneD, NekDouble> m_coefficientsB;
    Array<OneD, NekDouble> m_init;
    Array<OneD, NekDouble> m_f;
    Array<OneD, NekDouble> m_f1;
    Array<OneD, NekDouble> m_f2;
  };
  
  
} //end of namespace

#endif // NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_STREAMFUNCTIONWAVES_H
  
  
  
