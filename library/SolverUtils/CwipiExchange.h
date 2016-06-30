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

#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
namespace SolverUtils
{

class CwipiCoupling
{

public:
    CwipiCoupling(){};

    CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                  string name,
                  int outputFreq,
                  double geomTol);

    ~CwipiCoupling();

    MultiRegions::ExpListSharedPtr GetEvalField()
    {
        return m_evalField;
    };

    Array<OneD, MultiRegions::ExpListSharedPtr> GetRecvFields()
    {
        return m_recvFields;
    };

    string GetName()
    {
        return m_couplingName;
    };

    const std::map<std::string, std::string> GetConfig()
    {
        return m_config;
    }

    inline void FinalizeCoupling()
    {
        v_FinalizeCoupling();
    }

    void SendFields(const int step,
                    const NekDouble time,
                    Array<OneD, Array<OneD, NekDouble> > &field)
    {
        v_SendFields(step, time, field);
    }

    void ReceiveFields(const int step,
                       const NekDouble time,
                       Array<OneD, Array<OneD, NekDouble> > &field)
    {
        v_ReceiveFields(step, time, field);
    }

    void PrintProgressbar(const int position, const int goal) const;

protected:
    string m_couplingName;

    std::map<std::string, std::string> m_config;
    NekDouble m_filtWidth;

    MultiRegions::ExpListSharedPtr m_evalField;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_recvFields;

    int m_nSendVars;
    int m_nRecvVars;

    int m_spacedim;
    int m_nPoints;
    double *m_points;
    double *m_coords;
    int *m_connecIdx;
    int *m_connec;
    double *m_rValsInterl;
    double *m_sValsInterl;

    map<int, int> m_vertMap;

    virtual void v_FinalizeCoupling(void);

    virtual void v_SendFields(const int step,
                              const NekDouble time,
                              Array<OneD, Array<OneD, NekDouble> > &field);

    virtual void v_ReceiveFields(const int step,
                                 const NekDouble time,
                                 Array<OneD, Array<OneD, NekDouble> > &field);

private:
    void ReadConfig(LibUtilities::SessionReaderSharedPtr session);

    void SetupRecvFields();

    void AnnounceMesh();

    void AnnounceRecvPoints();

    template <typename T>
    void AddElementsToMesh(T geom,
                           int &coordsPos,
                           int &connecPos,
                           int &conidxPos);
};

typedef boost::shared_ptr<CwipiCoupling> CwipiCouplingSharedPointer;
}
}

#endif
