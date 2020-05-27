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

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/ForcingProgrammatic.h>
#include <MultiRegions/ExpList.h>
#include <string>

using namespace Nektar::SolverUtils;

#if defined(_MSC_VER) && defined(MoveFile)
#undef MoveFile
#endif

namespace Nektar
{     

    
    enum VWIIterationType
    {
        eFixedAlpha,
        eFixedWaveForcing,
        eFixedAlphaWaveForcing,
        eFixedWaveForcingWithSubIterationOnAlpha,
        eVWIIterationTypeSize
    };

    const std::string VWIIterationTypeMap[] = 
    {
        "FixedAlpha",
        "FixedWaveForcing",
        "FixedAlphaWaveForcing",
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
        void SaveLoopDetails(std::string dir, int i);

        void CalcNonLinearWaveForce(void);
        void CalcL2ToLinfPressure  (void);
        
        void SaveFile(std::string fileend, std::string dir, int n);
        void MoveFile(std::string fileend, std::string dir, int n);
        void CopyFile(std::string file1end, std::string file2end);


        bool CheckEigIsStationary(bool reset = false);
        bool CheckIfAtNeutralPoint(void);
        void UpdateAlpha(int n);
        void UpdateWaveForceMag(int n);
        void UpdateDAlphaDWaveForceMag(NekDouble alphainit);
        


        void AppendEvlToFile(std::string file, int n);
        void AppendEvlToFile(std::string file, NekDouble WaveForceMag);

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

        NekDouble GetAlpha(void)
        {
            return m_alpha[0];
        }
            

        NekDouble GetAlphaStep(void)
        {
            return m_alphaStep;
        }

        NekDouble GetWaveForceMag(void)
        {
            return m_waveForceMag[0];
        }

        NekDouble GetWaveForceMagStep(void)
        {
            return m_waveForceMagStep;
        }

        NekDouble GetDAlphaDWaveForceMag(void)
        {
            return m_dAlphaDWaveForceMag; 
        }

        int GetMaxWaveForceMagIter(void)
        {
            return m_maxWaveForceMagIter;
        }
        
        
	NekDouble GetEigRelTol(void)
            
	{
            return m_eigRelTol;
	}
	
	int GetMinInnerIterations(void)
	{
            return  m_minInnerIterations;
	}

        NekDouble GetPrevAlpha(void)
        {
            return m_prevAlpha;
        }

        void SetAlpha(NekDouble alpha)
        {
            m_alpha[0] = alpha; 
        }


        void SetWaveForceMag(NekDouble mag)
        {
            m_waveForceMag[0] = mag; 
        }

	void SetEigRelTol(NekDouble tol)
	{
	    m_eigRelTol = tol;
	}

        void  SetAlphaStep(NekDouble step)
        {
	    m_alphaStep = step;
        }

	void  SetMinInnerIterations(int niter)
	{
	  m_minInnerIterations = niter;
	}

        void SetPrevAlpha(NekDouble alpha)
        {
            m_prevAlpha = alpha; 
        }

        bool IfIterInterface(void)
        {
            return m_iterinterface;
        }

        Array<OneD, int> GetReflectionIndex(void);

        void FileRelaxation( int reg);

    protected:

    private:
        int m_iterStart; // Start iterations of inner loop
        int m_iterEnd;   // End iterations of inner loop
        
        int m_nOuterIterations; 
        int m_maxOuterIterations; // Maximum number of outer iterations        
        int m_minInnerIterations; // Minimum number of iterations in inner loop - based on relaxation factor
        int m_maxWaveForceMagIter; 

        bool m_deltaFcnApprox;  // Activate delta function approximation around wave 
        bool m_useLinfPressureNorm; // Activate if use Pressure Linf Normalisation

        bool m_moveMeshToCriticalLayer; // move mesh to critical layer 

        NekDouble m_deltaFcnDecay;   // Delta function decay level 

        Array<OneD, NekDouble>  m_waveForceMag;
        NekDouble m_waveForceMagStep;
        
        NekDouble m_rollForceScale; 

        Array<OneD, NekDouble> m_leading_real_evl;   /// < Leading real eigenvalue 
        Array<OneD, NekDouble> m_leading_imag_evl;   /// < Leading imaginary eigenvalue

        Array<OneD, NekDouble> m_alpha; 

        NekDouble m_alphaStep;
        NekDouble m_neutralPointTol; 
        NekDouble m_eigRelTol; 
        NekDouble m_vwiRelaxation; 
        NekDouble m_dAlphaDWaveForceMag;
        NekDouble m_prevAlpha;

        bool m_iterinterface;

        VWIIterationType m_VWIIterationType;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_waveVelocities;
        MultiRegions::ExpListSharedPtr              m_wavePressure;
        
        Array<OneD, Array<OneD, NekDouble > >  m_vwiForcing; 
        SolverUtils::ForcingProgrammaticSharedPtr m_vwiForcingObj;

        Array<OneD, Array<OneD, NekDouble > >  m_bcsForcing;
        
        Array<OneD, MultiRegions::ExpListSharedPtr> m_streakField; 

        Array<OneD, MultiRegions::ExpListSharedPtr> m_rollField; 
        
        std::string m_sessionName;
        LibUtilities::SessionReaderSharedPtr m_sessionVWI; 

        LibUtilities::SessionReaderSharedPtr m_sessionRoll;
        SpatialDomains::MeshGraphSharedPtr m_graphRoll;
        EquationSystemSharedPtr m_solverRoll;

        LibUtilities::SessionReaderSharedPtr m_sessionStreak;
        SpatialDomains::MeshGraphSharedPtr m_graphStreak;
        LibUtilities::SessionReaderSharedPtr m_sessionWave;
        SpatialDomains::MeshGraphSharedPtr m_graphWave;
    };
}

#endif
