///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSystemOutput.cpp
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
// Description: Output routines for  Pulse Wave Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <PulseWaveSolver/Utilities/PulseWaveSystemOutput.h>

using namespace std;

namespace Nektar
{
    
    string PulseWaveSystemOutput::className = GetEquationSystemFactory().RegisterCreatorFunction("PulseWaveSystemOutput", PulseWaveSystemOutput::create, "Pulse Wave Propagation output.");

    /**
     *  @class PulseWaveSystemOutput
     *
     *  Initialises the arterial subdomains in m_vessels and sets up
     *  all domain-linking conditions (bifurcations, junctions,
     *  merging flows). Detects the network structure and assigns
     *  boundary conditons. Also provides the underlying timestepping
     *  framework for pulse wave solvers including the general
     *  timestepping routines.
     */
    
    /**
     *  Processes SolverInfo parameters from the session file and sets
     *  up timestepping-specific code.
     *
     *  @param   m_Session        Session object to read parameters from.
     */
    PulseWaveSystemOutput::PulseWaveSystemOutput(const LibUtilities::SessionReaderSharedPtr& m_session)
        : PulseWaveSystem(m_session)
    {
    }
    
    /**
     *  Destructor
     */
    PulseWaveSystemOutput::~PulseWaveSystemOutput()
    {
    }
    
    /** o
     *
     * Duplicates PulseWaveSystem InitObject but does not set 
     * 
     */
    void PulseWaveSystemOutput::v_InitObject()
    {       
        m_session->RegisterCmdLineArgument("SetToOneSpaceDimension","False","Redefine mesh to be aligned to x-axis");
        
        PulseWaveSystem::v_InitObject();
    }
}
