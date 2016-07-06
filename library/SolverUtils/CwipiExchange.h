///////////////////////////////////////////////////////////////////////////////
//
// File CwipiExchange.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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
// Description: CWIPI Exchange class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_CWIPIEXCHANGE
#define NEKTAR_CWIPIEXCHANGE

#include <SolverUtils/Interpolator.h>
#include <SolverUtils/EquationSystem.h>

#include <cwipi.h>

namespace Nektar
{
namespace SolverUtils
{

typedef std::map<std::string, std::string> CouplingConfigMap;

class CwipiCoupling
{

public:
    CwipiCoupling(){};

    CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                  string name,
                  int outputFreq,
                  double geomTol);

    ~CwipiCoupling();

    const std::map<std::string, std::string> GetConfig()
    {
        return m_config;
    }

    std::vector<std::string> GetSendFieldNames()
    {
        return m_sendFieldNames;
    }

    std::vector<std::string> GetRecvFieldNames()
    {
        return m_recvFieldNames;
    }

    inline void FinalizeCoupling()
    {
        v_FinalizeCoupling();
    }

    void SendFields(const int step,
                    const NekDouble time,
                    const NekDouble timestep,
                    const Array<OneD, const Array<OneD, NekDouble> > &field);

    void ReceiveFields(const int step,
                       const NekDouble time,
                       const NekDouble timestep,
                       Array<OneD, Array<OneD, NekDouble> > &field);

    void SendCallback(Array<OneD, Array<OneD, NekDouble> > &interpField,
                      Array<OneD, Array<OneD, NekDouble> > distCoords);

    void PrintProgressbar(const int position, const int goal) const;

protected:
    string m_couplingName;

    CouplingConfigMap m_config;
    NekDouble m_filtWidth;

    Array<OneD, Array<OneD, NekDouble> > m_sendField;

    MultiRegions::ExpListSharedPtr m_evalField;
    MultiRegions::ExpListSharedPtr m_recvField;

    Array<OneD, Array<OneD, NekDouble> > m_oldFields;
    Array<OneD, Array<OneD, NekDouble> > m_newFields;

    int m_nSendVars;
    std::vector<std::string> m_sendFieldNames;
    int m_sendSteps;

    int m_nRecvVars;
    std::vector<std::string> m_recvFieldNames;
    int m_recvSteps;

    NekDouble m_lastUpdate;

    int m_spacedim;
    int m_nPoints;
    double *m_points;
    double *m_coords;
    int *m_connecIdx;
    int *m_connec;
    double *m_rValsInterl;
    double *m_sValsInterl;

    map<int, int> m_vertMap;

    SolverUtils::InterpolatorSharedPtr m_sendInterpolator;

    virtual void v_FinalizeCoupling(void);

    const NekDouble GetSendField(const int i, const int j) const
    {
        return m_sendField[i][j];
    }

private:
    void ReadConfig(LibUtilities::SessionReaderSharedPtr session);

    void SetupReceive();

    void SetupSend();

    void EvaluateFields(Array<OneD, Array<OneD, NekDouble> > interpField,
                        Array<OneD, Array<OneD, NekDouble> > distCoords);

    void SetupSendInterpolation();

    void AnnounceMesh();

    void DumpRawFields(const NekDouble time,
                       Array<OneD, Array<OneD, NekDouble> > rVals);

    template <typename T>
    void AddElementsToMesh(T geom,
                           int &coordsPos,
                           int &connecPos,
                           int &conidxPos);

    void FetchFields(const int step,
                     const NekDouble time,
                     Array<OneD, Array<OneD, NekDouble> > &field);
};

typedef boost::shared_ptr<CwipiCoupling> CwipiCouplingSharedPointer;

typedef boost::function<void(Array<OneD, Array<OneD, NekDouble> > interpField,
                             Array<OneD, Array<OneD, NekDouble> > distCoords)>
    SendCallbackType;
std::map<std::string, SendCallbackType> SenderCouplings;
}
}

#endif
