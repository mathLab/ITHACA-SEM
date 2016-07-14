///////////////////////////////////////////////////////////////////////////////
//
// File APE_coupled.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
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
// Description: Coupled APE1/APE4 (Acoustic Perturbation Equations)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_COUPLED_H
#define NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_COUPLED_H

#include <APESolver/EquationSystems/APE.h>

#include <SolverUtils/CwipiExchange.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class APE_coupled : public APE
{
public:
    friend class MemoryManager<APE_coupled>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        EquationSystemSharedPtr p =
            MemoryManager<APE_coupled>::AllocateSharedPtr(pSession);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Destructor
    virtual ~APE_coupled();

protected:
    SolverUtils::CwipiCouplingSharedPointer m_coupling;

    ForcingSharedPtr                        m_extForcing;

    /// Initialises UnsteadySystem class members.
    APE_coupled(const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual void v_InitObject();

    virtual bool v_PreIntegrate(int step);

    virtual void v_Output(void);

private:
    void receiveFields(int step);
};
}

#endif
