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
#include <SolverUtils/Exchange.h>

namespace Nektar
{
namespace SolverUtils
{

class CwipiCoupling : public SolverUtils::Coupling
{

public:

    CwipiCoupling()
    {
    };

    CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                     string name, string distAppname, int outputFreq, double geomTol, int oversamp = 0);

    ~CwipiCoupling();

    void PrintProgressbar(const int position, const int goal) const;

protected:

    string m_distAppname;
    string m_outputFormat;
    string m_outputFormatOption;
    int m_outputFreq;
    double m_geomTol;

    double *m_coords;
    int *m_connecIdx;
    int *m_connec;

    map< int, int > m_vertMap;

    virtual void v_FinalizeCoupling(void);

private:

    template <typename T>
    void AddElementsToMesh(T geom, int& coordsPos, int& connecPos, int& conidxPos);

    Array<OneD, Array<OneD, float> > m_weights;
    Array<OneD, Array<OneD, unsigned int> > m_neighbourInds;
};



typedef boost::shared_ptr<CwipiCoupling> CwipiCouplingSharedPointer;


class CwipiExchange : public SolverUtils::Exchange
{
public:

    CwipiExchange()
    {
    };

    CwipiExchange(SolverUtils::CouplingSharedPointer coupling, string name,
                     int nEVars);

    ~CwipiExchange();


protected:

    int m_nEVars;

    string m_sendFieldName;
    string m_recvFieldName;

    double *m_rValsInterl;


    virtual void v_SendFields(const int step, const NekDouble time,
                              Array<OneD, Array<OneD, NekDouble> > &field);

    virtual void v_ReceiveFields(const int step, const NekDouble time,
                                 Array<OneD, Array<OneD, NekDouble> > &field);


};


typedef boost::shared_ptr<CwipiExchange> CwipiExchangeSharedPtr;

}
}
#endif
