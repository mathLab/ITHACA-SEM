///////////////////////////////////////////////////////////////////////////////
//
// File Coupling.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Coupling class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COUPLING
#define NEKTAR_COUPLING

#include <FieldUtils/Interpolator.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace SolverUtils
{

class Coupling;

SOLVER_UTILS_EXPORT typedef std::shared_ptr<Coupling> CouplingSharedPtr;

/// Declaration of the Coupling factory
SOLVER_UTILS_EXPORT typedef LibUtilities::
    NekFactory<std::string, Coupling, MultiRegions::ExpListSharedPtr>
        CouplingFactory;

/// Declaration of the Coupling factory singleton
SOLVER_UTILS_EXPORT CouplingFactory &GetCouplingFactory();

class Coupling
{

public:
    typedef std::map<std::string, std::string> CouplingConfigMap;

    SOLVER_UTILS_EXPORT virtual ~Coupling()
    {
    }

    SOLVER_UTILS_EXPORT inline void Init()
    {
        v_Init();
    }

    SOLVER_UTILS_EXPORT inline const std::map<std::string, std::string> GetConfig()
    {
        return m_config;
    }

    SOLVER_UTILS_EXPORT inline std::vector<std::string> GetSendFieldNames()
    {
        return m_sendFieldNames;
    }

    SOLVER_UTILS_EXPORT inline std::vector<std::string> GetRecvFieldNames()
    {
        return m_recvFieldNames;
    }

    SOLVER_UTILS_EXPORT inline void Finalize()
    {
        v_Finalize();
    }

    SOLVER_UTILS_EXPORT inline void Send(
        const int step,
        const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames)
    {
        v_Send(step, time, field, varNames);
    }

    SOLVER_UTILS_EXPORT inline void Receive(
        const int step,
        const NekDouble time,
        Array<OneD, Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames)
    {
        v_Receive(step, time, field, varNames);
    }

protected:
    std::string m_couplingName;

    CouplingConfigMap m_config;

    MultiRegions::ExpListSharedPtr m_evalField;

    int m_nSendVars;
    std::vector<std::string> m_sendFieldNames;
    int m_sendSteps;

    int m_nRecvVars;
    std::vector<std::string> m_recvFieldNames;
    int m_recvSteps;

    SOLVER_UTILS_EXPORT Coupling(MultiRegions::ExpListSharedPtr field);

    SOLVER_UTILS_EXPORT virtual void v_Init();

    SOLVER_UTILS_EXPORT virtual void v_Send(
        const int step,
        const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames) = 0;

    SOLVER_UTILS_EXPORT virtual void v_Receive(
        const int step,
        const NekDouble time,
        Array<OneD, Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames) = 0;

    SOLVER_UTILS_EXPORT virtual void v_Finalize()
    {
    }

    SOLVER_UTILS_EXPORT std::vector<int> GenerateVariableMapping(
        std::vector<std::string> &vars, std::vector<std::string> &transVars);
};
}
}

#endif
