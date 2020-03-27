///////////////////////////////////////////////////////////////////////////////
//
// File: SubSteppingExtrapolateWeakPressure.h
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
// Description: Abstract base class for SubSteppingExtrapolateWeakPressure.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SUBSTEPPINGEXTRAPOLATEWEAKPRESSURE_H
#define NEKTAR_SOLVERS_SUBSTEPPINGEXTRAPOLATEWEAKPRESSURE_H

#include <IncNavierStokesSolver/EquationSystems/SubSteppingExtrapolate.h>

namespace Nektar
{
    //-------------
    // Substepping
    //-------------
    
    class SubSteppingExtrapolateWeakPressure;
    
    typedef std::shared_ptr<SubSteppingExtrapolateWeakPressure> SubSteppingExtrapolateWeakPressureSharedPtr;
    
    class SubSteppingExtrapolateWeakPressure : public SubSteppingExtrapolate
    {
    public:
        /// Creates an instance of this class
        static ExtrapolateSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            MultiRegions::ExpListSharedPtr              &pPressure,
            const Array<OneD, int>                      &pVel,
            const SolverUtils::AdvectionSharedPtr       &advObject)
        {
            ExtrapolateSharedPtr
                p = MemoryManager<SubSteppingExtrapolateWeakPressure>
                ::AllocateSharedPtr(pSession,pFields,pPressure,pVel,advObject);
            return p;
        }

        /// Name of class
        static std::string className;

        SubSteppingExtrapolateWeakPressure(
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            MultiRegions::ExpListSharedPtr              pPressure,
            const Array<OneD, int> pVel,
            const SolverUtils::AdvectionSharedPtr advObject);

        virtual ~SubSteppingExtrapolateWeakPressure();
        
    protected:
        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            NekDouble Aii_Dt,
            NekDouble kinvis);

        virtual void v_AddNormVelOnOBC(const int nbcoeffs, const int nreg,
                                       Array<OneD, Array<OneD, NekDouble> > &u);
    };
}
#endif

