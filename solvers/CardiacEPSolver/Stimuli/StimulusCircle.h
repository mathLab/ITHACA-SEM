///////////////////////////////////////////////////////////////////////////////
//
// File StimulusCircle.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Circular stimulus header file.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUSCIRC
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUSCIRC

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
//#include <SpatialDomains/SpatialData.h>
#include <MultiRegions/ExpList.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdNodalTetExp.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>

namespace Nektar
{

    /// Protocol base class.
    class StimulusCirc: public Stimulus
    {
    public:
        /// Creates an instance of this class
        static StimulusSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField,
                const TiXmlElement* pXml)
        {
            return MemoryManager<StimulusCirc>
                ::AllocateSharedPtr(pSession, pField, pXml);
        }

        /// Name of class
        static std::string className;

        friend class MemoryManager<StimulusCirc>;

        virtual ~StimulusCirc() {}

        /// Initialise the protocol storage and set initial conditions
        void Initialise();

    protected:
        NekDouble m_px1;
        NekDouble m_py1;
        NekDouble m_pz1;
        NekDouble m_pr1;
        NekDouble m_pis;
        NekDouble m_strength;
        NekDouble v_amp;
        NekDouble m_chiCapMembrane;

        virtual void v_Update(Array<OneD, Array<OneD, NekDouble> >&outarray,
                              const NekDouble time);

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

    private:
        StimulusCirc(const LibUtilities::SessionReaderSharedPtr& pSession,
                     const MultiRegions::ExpListSharedPtr& pField,
                     const TiXmlElement* pXml);
    };

}

#endif
