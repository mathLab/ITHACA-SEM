///////////////////////////////////////////////////////////////////////////////
//
// File VortexWaveInteraction.h
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
// Description: Vortex Wave Interaction class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VORTEXWAVEINTERACTION_H
#define NEKTAR_SOLVERS_VORTEXWAVEINTERACTION_H

#include <Auxiliary/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{     

    
    enum VWIIterationType
    {
        eFixedAlphaWaveForcing,
        eFixedWaveForcing,
        eFixedWaveForcingWithSubIterationOnAlpha,
        eVWIIterationTypeSize
    };

    const std::string VWIIterationTypeMap[] = 
    {
        "FixedAlphaWaveForcing",
        "FixedWaveForcing",
        "FixedWaveForcingWithSubIterationOnAlpha"
    };

    class VortexWaveInteraction
    {
    public:
        VortexWaveInteraction(int argc, char *argv[]);

        ~VortexWaveInteraction(void);
        
        void ExecuteRoll(void);
        void ExecuteStreak(void);
        void ExecuteWave(void);

        void ExecuteLoop(bool CalcWaveForce = true);
        void SaveLoopDetails(string dir, int i);

        void CalcNonLinearWaveForce(void);
        void SaveFile(string fileend, string dir, int n);
        void MoveFile(string fileend, string dir, int n);
        void CopyFile(string file1end, string file2end);


        bool CheckEigIsStationary(void);
        bool CheckIfAtNeutralPoint(void);
        void UpdateAlpha(int n);

        void AppendEvlToFile(std::string file, int n);

        int GetIterStart()
        {
            return m_iterStart;
        }

        int GetIterEnd()
        {
            return m_iterEnd;
        }

        VWIIterationType GetVWIIterationType(void)
        {
            return m_VWIIterationType;
        }
        
        int GetNOuterIterations(void)
        {
            return m_nOuterIterations;
        }

        int GetMaxOuterIterations(void)
        {
            return m_maxOuterIterations;
        }

        Array<OneD, int> GetReflectionIndex(void);

    protected:

    private:
        int m_iterStart; // Start iterations of inner loop
        int m_iterEnd;   // End iterations of inner loop
        
        int m_nOuterIterations; 
        int m_maxOuterIterations; // Maximum number of outer iterations        

        NekDouble m_waveForceMag;

        Array<OneD, NekDouble> m_leading_real_evl;   /// < Leading real eigenvalue 
        Array<OneD, NekDouble> m_leading_imag_evl;   /// < Leading imaginary eigenvalue

        Array<OneD, NekDouble> m_alpha; 

        NekDouble m_alphaStep;
        NekDouble m_neutralPointTol; 
        NekDouble m_eigRelTol; 
        NekDouble m_vwiRelaxation; 

        VWIIterationType m_VWIIterationType;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_waveVelocities;
        MultiRegions::ExpListSharedPtr              m_wavePressure;
        
        Array<OneD, Array<OneD, NekDouble > >  m_vwiForcing; 

        string m_sessionName;
        LibUtilities::SessionReaderSharedPtr m_sessionVWI; 

        LibUtilities::SessionReaderSharedPtr m_sessionRoll; 
        EquationSystemSharedPtr m_solverRoll;

        LibUtilities::SessionReaderSharedPtr m_sessionStreak; 
        LibUtilities::SessionReaderSharedPtr m_sessionWave; 
    };
}

#endif
