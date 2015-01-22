///////////////////////////////////////////////////////////////////////////////
//
// File FilterThresholdMax.h
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
// Description: Outputs time when solution first exceeds a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERTHRESHOLDMAX_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERTHRESHOLDMAX_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class FilterThresholdMax : public Filter
        {
        public:
            friend class MemoryManager<FilterThresholdMax>;

            /// Creates an instance of this class
            static FilterSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams) {
                FilterSharedPtr p = MemoryManager<FilterThresholdMax>::AllocateSharedPtr(pSession, pParams);
                //p->InitObject();
                return p;
            }

            ///Name of the class
            static std::string className;

            SOLVER_UTILS_EXPORT FilterThresholdMax(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams);
            SOLVER_UTILS_EXPORT virtual ~FilterThresholdMax();

        protected:
            virtual void v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual bool v_IsTimeDependent();

        private:
            Array<OneD, NekDouble> m_threshold;
            NekDouble m_startTime;
            NekDouble m_thresholdValue;
            int m_thresholdVar;
            NekDouble m_initialValue;
            std::string m_outputFile;
            LibUtilities::FieldIOSharedPtr m_fld;
        };
    }
}

#endif /* FILTERTHRESHOLDMAX_H_ */
