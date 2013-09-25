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
// Description: Abstract base class for SubSteppingExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SUPSTEPPINGEXTRAPOLATE_H
#define NEKTAR_SOLVERS_SUBSTEPPINGEXTRAPOLATE_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>
#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>

namespace Nektar
{
    //--------
    // Substepping
    // --------
    
    class SubSteppingExtrapolate;
    
    typedef boost::shared_ptr<SubSteppingExtrapolate> SubSteppingExtrapolateSharedPtr;
    
    class SubSteppingExtrapolate : public Extrapolate
    {
    public:

        /// Creates an instance of this class
        static ExtrapolateSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            Array<OneD, int> pVel,
            AdvectionTermSharedPtr advObject)
        {
            ExtrapolateSharedPtr p = MemoryManager<SubSteppingExtrapolate>::AllocateSharedPtr(pSession,pFields,pVel,advObject);
            return p;
        }

        /// Name of class
        static std::string className;

        SubSteppingExtrapolate(
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            Array<OneD, int> pVel,
            AdvectionTermSharedPtr advObject);

        virtual ~SubSteppingExtrapolate();
        
    protected:

        virtual void v_SubSteppingTimeIntegration(
            int intMethod);
 
        virtual void v_SubStepSaveFields(
            const int nstep);

        virtual void v_SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        const NekDouble Aii_Dt,
        NekDouble kinvis);

        virtual void v_SubStepAdvance(
            const int nstep, 
            NekDouble time);

        virtual void v_MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, Array<OneD, NekDouble> &Q, 
            Array<OneD, NekDouble> &Advection);
        
        void AddDuDt(
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble Aii_Dt);
        void AddDuDt2D(
            const Array<OneD, const Array<OneD, NekDouble> >  &N, 
            NekDouble Aii_Dt);
        void AddDuDt3D(
            const Array<OneD, const Array<OneD, NekDouble> >  &N, 
            NekDouble Aii_Dt);

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

        Array<OneD,NekDouble> GetMaxStdVelocity(
            const Array<OneD, Array<OneD,NekDouble> > inarray);

        LibUtilities::TimeIntegrationWrapperSharedPtr m_subStepIntegrationScheme;
        LibUtilities::TimeIntegrationSchemeOperators m_subStepIntegrationOps;
        LibUtilities::TimeIntegrationSolutionSharedPtr  m_integrationSoln;

        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;

        int nConvectiveFields;

        /// Get the number of coefficients
        int ncoeffs;

        Array<OneD, int>  vel_loc;

        NekDouble cflSafetyFactor;


        //int m_intSteps;

        /// Number of time steps between outputting status information.
        //int m_infosteps;

        /// CFL safety factor (comprise between 0 to 1).
        //NekDouble m_cflSafetyFactor;

        /// Time step size
        //NekDouble m_timestep;

        //bool m_HalfMode;

        //bool m_SingleMode;
        
    };
}

#endif

