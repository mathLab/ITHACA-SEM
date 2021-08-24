///////////////////////////////////////////////////////////////////////////////
//
// File LinearisedAdvection.h
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
// Description: TBA
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_NAVIERSTOKESADVECTION_H
#define NEKTAR_SOLVERS_NAVIERSTOKESADVECTION_H

#include <SolverUtils/Advection/Advection.h>


namespace Nektar
{     

class NavierStokesAdvection: public SolverUtils::Advection

{
public:
    friend class MemoryManager<NavierStokesAdvection>;

    /// Creates an instance of this class
    static SolverUtils::AdvectionSharedPtr create(std::string) {
        return MemoryManager<NavierStokesAdvection>::AllocateSharedPtr();
    }

    /// Name of class
    static std::string className;
    static std::string className2;

    void SetSpecHPDealiasing(bool value)
    {
        m_specHP_dealiasing = value;
    }

protected:

    NavierStokesAdvection();

    virtual ~NavierStokesAdvection();

    virtual void v_InitObject(
              LibUtilities::SessionReaderSharedPtr         pSession,
              Array<OneD, MultiRegions::ExpListSharedPtr>  pFields);

    virtual void v_Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayOfArray,
        const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayOfArray);

        static std::string navierStokesAdvectionTypeLookupIds[];

private:
    bool m_specHP_dealiasing;
    bool m_homogen_dealiasing;
    bool m_SingleMode;
    bool m_HalfMode;
};
    
    
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
