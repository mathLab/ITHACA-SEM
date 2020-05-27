///////////////////////////////////////////////////////////////////////////////
//
// File Protocol.cpp
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
// Description: Protocol base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <CardiacEPSolver/Stimuli/Protocol.h>

namespace Nektar
{
    ProtocolFactory& GetProtocolFactory()
    {
        static ProtocolFactory instance;
        return instance;
    }

    /**
     * @class Protocol
     *
     * The Stimuli class and derived classes implement a range of stimuli.
     * The stimulus contains input stimuli that can be applied throughout the
     * domain, on specified regions determined by the derived classes of
     * Stimulus, at specified frequencies determined by the derived classes of
     * Protocol.
     */

    /**
     * Protocol base class constructor.
     */
    Protocol::Protocol(const LibUtilities::SessionReaderSharedPtr& pSession,
                       const TiXmlElement* pXml)
        : m_session(pSession)
    {
    }


    /**
     * Initialise the protocol. Allocate workspace and variable storage.
     */
    void Protocol::Initialise()
    {
    }

}
