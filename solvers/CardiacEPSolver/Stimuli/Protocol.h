///////////////////////////////////////////////////////////////////////////////
//
// File Protocol.h
//
// For more information, please see: http://www.nektar.info
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
// Description: Stimulus protocol base class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOL
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_PROTOCOL

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/Core/Misc.h>

namespace Nektar
{
    // Forward declaration
    class Protocol;

    /// A shared pointer to an EquationSystem object
    typedef std::shared_ptr<Protocol> ProtocolSharedPtr;

    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory< std::string, Protocol,
                const LibUtilities::SessionReaderSharedPtr&,
                const TiXmlElement*> ProtocolFactory;

    ProtocolFactory& GetProtocolFactory();


    /// Protocol base class.
    class Protocol
    {
    public:
        virtual ~Protocol() {}

        /// Initialise the protocol storage and set initial conditions
        void Initialise();

        /// Returns amplitude of stimulus (1 or 0) at given time
        NekDouble GetAmplitude(const NekDouble time)
        {
            return v_GetAmplitude(time);
        }

        /// Print a summary of the cell model
        void GenerateSummary(SolverUtils::SummaryList& s)
        {
            v_GenerateSummary(s);
        }

    protected:
        /// Session
        LibUtilities::SessionReaderSharedPtr m_session;

        Protocol(const LibUtilities::SessionReaderSharedPtr& pSession,
                 const TiXmlElement* pXml);

        virtual NekDouble v_GetAmplitude(const NekDouble time) = 0;

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s) = 0;

    };

}

#endif /*PROTOCOL_H_ */
