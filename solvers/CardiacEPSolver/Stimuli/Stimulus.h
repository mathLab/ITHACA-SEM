///////////////////////////////////////////////////////////////////////////////
//
// File Stimulus.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Stimulus base class header.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUS
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUS

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <CardiacEPSolver/Stimuli/Protocol.h>
#include <SolverUtils/Core/Misc.h>

namespace Nektar
{
    // Forward declaration
    class Stimulus;

    /// A shared pointer to an EquationSystem object
    typedef std::shared_ptr<Stimulus> StimulusSharedPtr;

    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory< std::string, Stimulus,
                const LibUtilities::SessionReaderSharedPtr&,
                const MultiRegions::ExpListSharedPtr&,
                const TiXmlElement*> StimulusFactory;

    StimulusFactory& GetStimulusFactory();


    /// Stimulus base class.
    class Stimulus
    {
    public:
        virtual ~Stimulus() {}

        /// Initialise the stimulus storage and set initial conditions
        void Initialise();

        /// Updates RHS of outarray by adding a stimulus to it
        void Update(Array<OneD, Array<OneD, NekDouble> >&outarray,
                    const NekDouble time)
        {
            v_Update(outarray, time);
        }

        /// Print a summary of the outarray
        void GenerateSummary(SolverUtils::SummaryList& s)
        {
            v_GenerateSummary(s);
        }

        static std::vector<StimulusSharedPtr> LoadStimuli(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const MultiRegions::ExpListSharedPtr& pField);

    protected:
        /// Session
        LibUtilities::SessionReaderSharedPtr m_session;
        /// Transmembrane potential field from PDE system
        MultiRegions::ExpListSharedPtr m_field;
        /// Number of physical points.
        int m_nq;
        /// Stimulus protocol to apply
        ProtocolSharedPtr m_Protocol;

        Stimulus(const LibUtilities::SessionReaderSharedPtr& pSession,
                  const MultiRegions::ExpListSharedPtr& pField,
                  const TiXmlElement* pXml);

        virtual void v_Update(Array<OneD, Array<OneD, NekDouble> >&outarray,
                              const NekDouble time) = 0;

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s) = 0;

    };
}

#endif
