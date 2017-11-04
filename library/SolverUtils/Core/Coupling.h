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
// Description: Coupling class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COUPLING
#define NEKTAR_COUPLING

#include <FieldUtils/Interpolator.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace SolverUtils
{

class Coupling;

typedef std::shared_ptr<Coupling> CouplingSharedPointer;

/// Declaration of the Coupling factory
typedef LibUtilities::NekFactory<std::string,
                                 Coupling,
                                 MultiRegions::ExpListSharedPtr>
    CouplingFactory;

/// Declaration of the Coupling factory singleton
CouplingFactory &GetCouplingFactory();

class Coupling
{

public:
    typedef std::map<std::string, std::string> CouplingConfigMap;

    virtual ~Coupling(){};

    inline void Init()
    {
        v_Init();
    };

    inline const std::map<std::string, std::string> GetConfig()
    {
        return m_config;
    }

    inline std::vector<std::string> GetSendFieldNames()
    {
        return m_sendFieldNames;
    }

    inline std::vector<std::string> GetRecvFieldNames()
    {
        return m_recvFieldNames;
    }

    inline void Finalize()
    {
        v_Finalize();
    };

    inline void Send(const int step,
                     const NekDouble time,
                     const Array<OneD, const Array<OneD, NekDouble> > &field,
                     vector<string> &varNames)
    {
        v_Send(step, time, field, varNames);
    };

    inline void Receive(const int step,
                        const NekDouble time,
                        Array<OneD, Array<OneD, NekDouble> > &field,
                        vector<string> &varNames)
    {
        v_Receive(step, time, field, varNames);
    };

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

    Coupling(MultiRegions::ExpListSharedPtr field);

    virtual void v_Init();

    virtual void v_Send(const int step,
                        const NekDouble time,
                        const Array<OneD, const Array<OneD, NekDouble> > &field,
                        vector<string> &varNames) = 0;

    virtual void v_Receive(
        const int step,
        const NekDouble time,
        Array<OneD, Array<OneD, NekDouble> > &field,
        vector<string> &varNames) = 0;

    virtual void v_Finalize()
    {
    };

    std::vector<int> GenerateVariableMapping(std::vector<std::string> &vars, std::vector<std::string> &transVars);

};
}
}

#endif
