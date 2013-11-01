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
// Description: Moving Body (Wavyness and acceleration)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY

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
    class ForcingMovingBody : public Forcing
    {
        public:

            friend class MemoryManager<ForcingMovingBody>;

            /// Creates an instance of this class
            SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce)
            {
                ForcingSharedPtr p = MemoryManager<ForcingMovingBody>::
                                                AllocateSharedPtr(pSession);
                p->InitObject(pFields, pNumForcingFields, pForce);
                return p;
            }

            ///Name of the class
            static std::string className;

        protected:
		
            SOLVER_UTILS_EXPORT virtual void v_InitObject(const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                                          const unsigned int& pNumForcingFields,
                                                          const TiXmlElement* pForce);

            SOLVER_UTILS_EXPORT virtual void v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                                     const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                     Array<OneD, Array<OneD, NekDouble> > &outarray,
                                                     NekDouble time);

        private:
        
            ForcingMovingBody(const LibUtilities::SessionReaderSharedPtr& pSession);
		
            void CheckIsFromFile();
		
            void UpdateMotion(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                              NekDouble time);
		
            void CalculateForcing(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);
        
            Array<OneD, Array<OneD, NekDouble> >    m_zeta;
            Array<OneD, Array<OneD, NekDouble> >    m_eta;
        
            Array<OneD, std::string> m_funcName;     // [0] is displacements -- [1] is velocities -- [2] is accelerations
            Array<OneD, std::string> m_motion;       // motion direction: [0] is 'x' and [1] is 'y'
            Array<OneD, bool>        m_IsFromFile;   // do determine if the the body motion come from an extern file
            
            NekDouble m_kinvis;
    };

}
}

#endif
