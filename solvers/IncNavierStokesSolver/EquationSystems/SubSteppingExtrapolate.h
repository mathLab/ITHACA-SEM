///////////////////////////////////////////////////////////////////////////////
//
// File: SubSteppingExtrapolate.h
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
// Description: Abstract base class for SubSteppingExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SUBSTEPPINGEXTRAPOLATE_H
#define NEKTAR_SOLVERS_SUBSTEPPINGEXTRAPOLATE_H

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/AdvectionSystem.h>

namespace Nektar
{
    //--------
    // Substepping
    // --------
    
    class SubSteppingExtrapolate;
    
    typedef std::shared_ptr<SubSteppingExtrapolate> SubSteppingExtrapolateSharedPtr;
    
    class SubSteppingExtrapolate : public Extrapolate
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
            ExtrapolateSharedPtr p = MemoryManager<SubSteppingExtrapolate>
                ::AllocateSharedPtr(pSession,pFields,pPressure,pVel,advObject);
            return p;
        }

        /// Name of class
        static std::string className;

        SubSteppingExtrapolate(
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            MultiRegions::ExpListSharedPtr              pPressure,
            const Array<OneD, int> pVel,
            const SolverUtils::AdvectionSharedPtr advObject);

        virtual ~SubSteppingExtrapolate();
        
    protected:
        virtual void v_EvaluatePressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &fields,
                                           const Array<OneD, const Array<OneD, NekDouble> >  &N,
                                           NekDouble kinvis);

        virtual void v_SubSteppingTimeIntegration(
            const LibUtilities::TimeIntegrationSchemeSharedPtr & IntegrationScheme );
 
        virtual void v_SubStepSaveFields(
            int nstep);

        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            NekDouble Aii_Dt,
            NekDouble kinvis);

        virtual void v_SubStepAdvance(
            const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr &integrationSoln, 
            int nstep, 
            NekDouble time);

        virtual void v_MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection);
        
        virtual std::string v_GetSubStepName(void);

        void SubStepAdvection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
            Array<OneD, Array<OneD,NekDouble> > &outarray,
            const NekDouble time);

        void SubStepProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
            Array<OneD, Array<OneD, NekDouble> > &outarray, 
            const NekDouble time);

        void SubStepExtrapolateField(
            NekDouble toff, 
            Array< OneD, Array<OneD, NekDouble> > &ExtVel);

        void AddAdvectionPenaltyFlux(
            const Array<OneD, const Array<OneD, NekDouble> > &velfield,
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,Array<OneD, 
            Array<OneD, NekDouble> > &outarray);

        NekDouble GetSubstepTimeStep();

        LibUtilities::TimeIntegrationSchemeSharedPtr m_subStepIntegrationScheme;
        LibUtilities::TimeIntegrationSchemeOperators m_subStepIntegrationOps;

        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;

        NekDouble m_cflSafetyFactor;
        int m_infosteps;
        int m_minsubsteps;
        int m_maxsubsteps;
    };
}

#endif

