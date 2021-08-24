///////////////////////////////////////////////////////////////////////////////
//
// File ProtocolSingle.h
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
// Description: Protocol Single stimulus header.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOLSINGLE
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOLSINGLE

#include <CardiacEPSolver/Stimuli/Protocol.h>

namespace Nektar
{
    // Forward declaration
    class ProtocolSingle;

    /// Protocol base class.
    class ProtocolSingle : public Protocol
    {
    public:
        /// Creates an instance of this class
        static ProtocolSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const TiXmlElement* pXml)
        {
            return MemoryManager<ProtocolSingle>
                    ::AllocateSharedPtr(pSession, pXml);
        }

        /// Name of class
        static std::string className;

        friend class MemoryManager<ProtocolSingle>;

        virtual ~ProtocolSingle() {}

        /// Initialise the protocol storage and set initial conditions
        void Initialise();

    protected:
        NekDouble m_start;
        NekDouble m_dur;

        virtual NekDouble v_GetAmplitude(const NekDouble time);

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        virtual void v_SetInitialConditions();

    private:
        ProtocolSingle(const LibUtilities::SessionReaderSharedPtr& pSession,
                       const TiXmlElement* pXml);
    };

}

#endif
