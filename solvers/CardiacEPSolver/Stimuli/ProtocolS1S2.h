///////////////////////////////////////////////////////////////////////////////
//
// File ProtocolS1S2.h
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
// Description: Protocol S1S2 stimulus.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOLS1S2
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOLS1S2
#include <CardiacEPSolver/Stimuli/Protocol.h>

namespace Nektar
{
    // Forward declaration
    class ProtocolS1S2;

    /// Protocol base class.
    class ProtocolS1S2 : public Protocol
    {
    public:
        /// Creates an instance of this class
        static ProtocolSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const TiXmlElement* pXml)
        {
            return MemoryManager<ProtocolS1S2>
                    ::AllocateSharedPtr(pSession, pXml);
        }

        /// Name of class
        static std::string className;

        friend class MemoryManager<ProtocolS1S2>;

        virtual ~ProtocolS1S2() {}

        /// Initialise the protocol storage and set initial conditions
        void Initialise();

    protected:
        NekDouble m_start;
        NekDouble m_dur;
        NekDouble m_s1cyclelength;
        NekDouble m_num_s1;
        NekDouble m_s2cyclelength;
        NekDouble m_s2start;

        virtual NekDouble v_GetAmplitude(const NekDouble time);

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        virtual void v_SetInitialConditions();

    private:
        ProtocolS1S2(const LibUtilities::SessionReaderSharedPtr& pSession,
                     const TiXmlElement* pXml);
    };

}

#endif
