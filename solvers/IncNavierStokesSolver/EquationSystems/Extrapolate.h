///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.h
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_EXTRAPOLATE_H
#define NEKTAR_SOLVERS_EXTRAPOLATE_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>

namespace Nektar
{
    // Forward declaration
    class Extrapolate;

    typedef boost::shared_ptr<Extrapolate> ExtrapolateSharedPtr;
    
    typedef LibUtilities::NekFactory< std::string, Extrapolate,
        const LibUtilities::SessionReaderSharedPtr& ,
        Array<OneD, MultiRegions::ExpListSharedPtr>& ,
        Array<OneD, int> > ExtrapolateFactory;

    ExtrapolateFactory& GetExtrapolateFactory();


    class Extrapolate
    {
    public:
        Extrapolate(        
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            Array<OneD, int> pVel);

        virtual ~Extrapolate();

        inline void SubSteppingTimeIntegration(
            int intMethod,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields);
        inline void SubStepSaveFields(const int nstep);
        inline void SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            const NekDouble Aii_DT);
        inline void SubStepAdvance(const int nstep, NekDouble m_time);
        inline void AddDuDt(
            const Array<OneD, const Array<OneD, NekDouble> >  &N, 
            NekDouble Aii_Dt);

    protected:
        virtual void v_SubSteppingTimeIntegration(
            int intMethod,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields) = 0;
        virtual void v_SubStepSaveFields(const int nstep)=0;
        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT)=0;
        virtual void v_SubStepAdvance(const int nstep, NekDouble m_time)=0;
        virtual void v_AddDuDt(
            const Array<OneD, const Array<OneD, NekDouble> >  &N, 
            NekDouble Aii_Dt)=0;
        
        LibUtilities::SessionReaderSharedPtr        m_session;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

        Array<OneD, int> m_velocity;

        LibUtilities::CommSharedPtr m_comm;

        AdvectionTermSharedPtr m_advObject;


        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;

    private:
        static std::string def;
    };

    /**
     *
     */
    inline void Extrapolate::SubSteppingTimeIntegration(
        int intMethod,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
    {
        v_SubSteppingTimeIntegration(intMethod,pFields);
    }

    /**
     *
     */
    inline void Extrapolate::SubStepSaveFields(const int nstep)
    {
        v_SubStepSaveFields(nstep);
    }

    /**
     *
     */
    inline void Extrapolate::SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT)
    {
        v_SubStepSetPressureBCs(inarray,Aii_DT);
    }

    /**
     *
     */
    inline void Extrapolate::AddDuDt(
        const Array<OneD, const Array<OneD, NekDouble> >  &N, 
        NekDouble Aii_Dt)
    {
        v_AddDuDt(N,Aii_Dt);
    }

    /**
     *
     */
    inline void Extrapolate::SubStepAdvance(const int nstep, NekDouble m_time)
    {
        v_SubStepAdvance(nstep, m_time);
    }
}

#endif

