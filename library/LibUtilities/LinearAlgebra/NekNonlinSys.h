///////////////////////////////////////////////////////////////////////////////
//
// File  NekNonlinSys.h
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
// Description: NekNonlinSys header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H

#include <LibUtilities/LinearAlgebra/NekSys.h>

namespace Nektar
{
    namespace LibUtilities
    {
        class NekNonlinSys;

        typedef std::shared_ptr<NekNonlinSys> NekNonlinSysSharedPtr;
        
		typedef LibUtilities::NekFactory<
            std::string,
            NekNonlinSys,
            const LibUtilities::SessionReaderSharedPtr  &,
            const LibUtilities::CommSharedPtr           &,
            const int                                   > NekNonlinSysFactory;
        LIB_UTILITIES_EXPORT NekNonlinSysFactory& GetNekNonlinSysFactory();

        class  NekNonlinSys : public NekSys 
        {
            public:
                friend class MemoryManager<NekNonlinSys>;
                LIB_UTILITIES_EXPORT static NekNonlinSysSharedPtr CreateInstance
                (
                    const LibUtilities::SessionReaderSharedPtr  &pSession,
                    const LibUtilities::CommSharedPtr           &vComm,
                    const int                                   nDimen)
                {
                    NekNonlinSysSharedPtr p = MemoryManager<
                        NekNonlinSys>::AllocateSharedPtr(pSession,
                                                         vComm, nDimen);
                    return p;
                }
                LIB_UTILITIES_EXPORT NekNonlinSys(
                    const LibUtilities::SessionReaderSharedPtr  &pSession,
                    const LibUtilities::CommSharedPtr           &vComm,
                    const int                                   nDimen);
                LIB_UTILITIES_EXPORT ~NekNonlinSys();

                LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> 
                & GetRefSolution() const
                {
                    return m_Solution;
                }

                LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> 
                & GetRefResidual() const
                {
                    return m_Residual;
                }
                
            protected:
                /// Maximum iterations
                int                                       m_maxiter;
                /// Tolerance of iterative solver.
                NekDouble                                 m_tolerance;

                int                                       m_totalIterations = 0;

                Array<OneD, NekDouble>  m_Solution;
                Array<OneD, NekDouble>  m_Residual;
                Array<OneD, NekDouble>  m_DeltSltn;

                virtual void v_InitObject();
            private:
        };
    }
}
#endif
