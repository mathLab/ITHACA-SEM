///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForces.h
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERMOVINGBODY_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERMOVINGBODY_H

#include <SolverUtils/Filters/Filter.h>
#include <LocalRegions/Expansion3D.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class FilterMovingBody : public Filter
        {
        public:
            friend class MemoryManager<FilterMovingBody>;

            /// Creates an instance of this class
            static FilterSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams)
            {
                FilterSharedPtr p = MemoryManager<FilterMovingBody>::
                                        AllocateSharedPtr(pSession, pParams);
                //p->InitObject();
                return p;
            }

            ///Name of the class
            static std::string className;

            SOLVER_UTILS_EXPORT FilterMovingBody(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams);
            SOLVER_UTILS_EXPORT ~FilterMovingBody();

        protected:
            virtual void v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual bool v_IsTimeDependent();
            virtual void v_GetMotionVars(const Array<OneD, NekDouble> &inArray);
            virtual Array<OneD, NekDouble> v_GetAeroForces()
            {
                return m_Aeroforces;
            }
        private:
            /// ID's of boundary regions where we want the forces
            vector<unsigned int>            m_boundaryRegionsIdList;
            /// Determines if a given Boundary Region is in
            /// m_boundaryRegionsIdList
            vector<bool>                    m_boundaryRegionIsInList;
            unsigned int                    m_index;
            unsigned int                    m_outputFrequency;
            /// plane to take history point from if using a homogeneous1D
            /// expansion
            unsigned int                    m_outputPlane;
            bool                            m_isHomogeneous1D;
            LibUtilities::BasisSharedPtr    m_homogeneousBasis;
            std::string                     m_BoundaryString;
            /// number of planes for homogeneous1D expansion
            int                             m_planes;
            Array<OneD, NekDouble>          m_Aeroforces;
            Array<OneD, NekDouble>          m_MotionVars;
            Array<OneD, std::ofstream>      m_outputStream;
            std::string                        m_outputFile_fce;
            std::string                        m_outputFile_mot;
        };
    }
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
