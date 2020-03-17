///////////////////////////////////////////////////////////////////////////////
//
// File: Forcing.h
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCING
#define NEKTAR_SOLVERUTILS_FORCING

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Core/SessionFunction.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace SolverUtils
{
    //  Forward declaration
    class Forcing;
    class EquationSystem;

    /// A shared pointer to an EquationSystem object
    SOLVER_UTILS_EXPORT typedef std::shared_ptr<Forcing> ForcingSharedPtr;

    /// Declaration of the forcing factory
    typedef LibUtilities::NekFactory<std::string, Forcing,
            const LibUtilities::SessionReaderSharedPtr&,
            const std::weak_ptr<EquationSystem>&,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&,
            const unsigned int&,
            const TiXmlElement*> ForcingFactory;

    /// Declaration of the forcing factory singleton
    SOLVER_UTILS_EXPORT ForcingFactory& GetForcingFactory();

    /**
     * @class Forcing
     * @brief Defines a forcing term to be explicitly applied.
     */
    class Forcing
    {
        public:
            SOLVER_UTILS_EXPORT virtual ~Forcing() {}

            /// Initialise the forcing object
            SOLVER_UTILS_EXPORT void InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&       pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce);

            /// Apply the forcing
            SOLVER_UTILS_EXPORT void Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const NekDouble                                   &time);

            SOLVER_UTILS_EXPORT static std::vector<ForcingSharedPtr> Load(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const std::weak_ptr<EquationSystem>      &pEquation,
                        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                        const unsigned int& pNumForcingFields = 0);

            SOLVER_UTILS_EXPORT const Array<OneD, const Array<OneD, NekDouble> >
                &GetForces();

            SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, NekDouble> > &UpdateForces();

        protected:
            /// Session reader
            LibUtilities::SessionReaderSharedPtr m_session;
            /// Weak pointer to equation system using this forcing
            const std::weak_ptr<EquationSystem> m_equ;
            /// Evaluated forcing function
            Array<OneD, Array<OneD, NekDouble> > m_Forcing;
            /// Number of variables
            int m_NumVariable;
            /// Map of known SessionFunctions
            std::map<std::string, SolverUtils::SessionFunctionSharedPtr> m_sessionFunctions;

            /// Constructor
            SOLVER_UTILS_EXPORT Forcing(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::weak_ptr<EquationSystem>      &pEquation);

            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&       pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce) = 0;

            SOLVER_UTILS_EXPORT virtual void v_Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const NekDouble &time)=0;

                /// Get a SessionFunction by name
            SOLVER_UTILS_EXPORT SessionFunctionSharedPtr GetFunction(
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &pFields,
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                std::string                                         pName,
                bool                                                pCache = false);

            SOLVER_UTILS_EXPORT void EvaluateTimeFunction(
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    std::string                                 pFieldName, 
                    Array<OneD, NekDouble>&                     pArray,
                    const std::string&                          pFunctionName,
                    NekDouble pTime = NekDouble(0));


            SOLVER_UTILS_EXPORT void EvaluateTimeFunction(
                    const NekDouble                          pTime,
                    const LibUtilities::EquationSharedPtr   &pEqn,
                    Array<OneD, NekDouble>&                  pArray);
    };
}
}

#endif
