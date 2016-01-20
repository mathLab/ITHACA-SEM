///////////////////////////////////////////////////////////////////////////////
//
// File Exchange.h
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
// Description: openPALM Exchange class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_EXCHANGE
#define NEKTAR_SOLVERUTILS_EXCHANGE

#include <string>

#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>

using namespace std;

namespace Nektar
{

namespace SolverUtils
{

class Coupling
{

public:

    Coupling()
    {
    };

    Coupling(MultiRegions::ExpListSharedPtr field, string name) :
        m_field(field),
        m_points(NULL),
        m_name(name)
    {
    };

    ~Coupling()
    {
    };

    MultiRegions::ExpListSharedPtr GetField()
    {
        return m_field;
    };

    string GetName()
    {
        return m_name;
    };

    int GetNPoints()
    {
        return m_nPoints;
    };

    void GetPoints(double *points)
    {
        points = m_points;
    }

    inline void FinalizeCoupling()
    {
        v_FinalizeCoupling();
    }

protected:

    MultiRegions::ExpListSharedPtr m_field;
    int m_nPoints;
    double *m_points;

    string m_name;

    virtual void v_FinalizeCoupling(void)
    {
        ASSERTL0(false, "This method is not valid for the base class");
    };
};

typedef boost::shared_ptr<Coupling> CouplingSharedPointer;


class Exchange
{
public:

    Exchange()
    {
    };

    Exchange(CouplingSharedPointer coupling, string name) :
        m_coupling(coupling),
        m_name(name)
    {
    };

    ~Exchange()
    {
    };

    virtual void SendFields(const int step, const NekDouble time,
                                Array<OneD, Array<OneD, NekDouble> > &field)
    {
        v_SendFields(step, time, field);
    }

    virtual void ReceiveFields(const int step, const NekDouble time,
                                Array<OneD, Array<OneD, NekDouble> > &field)
    {
        v_ReceiveFields(step, time, field);
    }

protected:

    CouplingSharedPointer m_coupling;

    string m_name;

    virtual void v_SendFields(const int step, const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > &field)
    {
        ASSERTL0(false, "not valid for this class");
    };

    virtual void v_ReceiveFields(const int step, const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > &field)
    {
        ASSERTL0(false, "not valid for this class");
    };

};


typedef boost::shared_ptr<Exchange> ExchangeSharedPtr;

}
}
#endif
