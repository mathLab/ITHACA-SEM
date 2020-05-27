///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingStabilityCoupledLNS.h
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
// Description: Forcing to copy velocity field in steady coupled LNS stability 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGSTABILITYCOUPLEDLNS
#define NEKTAR_SOLVERUTILS_FORCINGSTABILITYCOUPLEDLNS

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>

namespace Nektar
{

class ForcingStabilityCoupledLNS : public SolverUtils::Forcing
{
    public:

        friend class MemoryManager<ForcingStabilityCoupledLNS>;

        /// Creates an instance of this class
        static SolverUtils::ForcingSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce)
        {
            SolverUtils::ForcingSharedPtr p =
                MemoryManager<ForcingStabilityCoupledLNS>::
                AllocateSharedPtr(pSession, pEquation);
            p->InitObject(pFields, pNumForcingFields, pForce);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:
        virtual void v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int&                         pNumForcingFields,
            const TiXmlElement*                         pForce);

        virtual void v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
            const Array<OneD, Array<OneD, NekDouble> >& inarray,
                  Array<OneD, Array<OneD, NekDouble> >& outarray,
            const NekDouble&                            time);

    private:
        ForcingStabilityCoupledLNS(
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                const std::weak_ptr<SolverUtils::EquationSystem> &pEquation);
};

}

#endif
