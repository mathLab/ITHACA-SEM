///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.h
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
// Description: Body or FIeld forcing
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGBODY

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>


namespace Nektar
{
namespace SolverUtils
{
    class ForcingBody : public Forcing
    {
        public:

            friend class MemoryManager<ForcingBody>;

            /// Creates an instance of this class
            SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::weak_ptr<EquationSystem>      &pEquation,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce)
            {
                ForcingSharedPtr p = MemoryManager<ForcingBody>::
                                        AllocateSharedPtr(pSession, pEquation);
                p->InitObject(pFields, pNumForcingFields, pForce);
                return p;
            }

            ///Name of the class
            static std::string classNameBody;
            static std::string classNameField;

        protected:
            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce);

            SOLVER_UTILS_EXPORT virtual void v_Apply(
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                    const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble &time);
                    
            SOLVER_UTILS_EXPORT virtual void v_ApplyCoeff(
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >              &outarray,
                    const NekDouble                                   &time);
        private:
            std::string                     m_funcName;
            bool                            m_hasTimeFcnScaling;
            LibUtilities::EquationSharedPtr m_timeFcnEqn;
            bool                            m_transform;

            ForcingBody(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::weak_ptr<EquationSystem>      &pEquation);

            virtual ~ForcingBody(void){};

            void Update(
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const Array<OneD, Array<OneD, NekDouble> > &inarray,
                const NekDouble &time);
    };

}
}

#endif
