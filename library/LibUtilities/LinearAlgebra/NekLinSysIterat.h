///////////////////////////////////////////////////////////////////////////////
//
// File  NekLinSysIterat.h
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
// Description: NekLinSysIterat header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_H

#include <LibUtilities/LinearAlgebra/NekNonlinLinSys.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
namespace Nektar
{
    namespace LibUtilities
    {

        class  NekLinSysIterat;

        typedef std::shared_ptr<NekLinSysIterat> NekLinSysIteratSharedPtr;

        typedef LibUtilities::NekFactory<
            std::string,
            NekLinSysIterat,
            const LibUtilities::SessionReaderSharedPtr  &,
            const LibUtilities::CommSharedPtr           &,
            const int                                    > NekLinSysIteratFactory;
        LIB_UTILITIES_EXPORT NekLinSysIteratFactory& GetNekLinSysIteratFactory();
        
        class  NekLinSysIterat : public NekNonlinLinSys
        {
        public:

            /// Support creation through MemoryManager.
            friend class MemoryManager<NekLinSysIterat>;
            /**
             * @brief Creates an instance of the SessionReader class.
             *
             * This function should be used by an application to instantiate the
             * session reader. It should be called at the very beginning of the
             * application before any other processing of command-line
             * arguments. After instantiating the class and setting up any
             * parallel communication, it also calls the main initialisation
             * of the object.
             */
            LIB_UTILITIES_EXPORT static NekLinSysIteratSharedPtr CreateInstance(
                const LibUtilities::SessionReaderSharedPtr  &pSession,
                const LibUtilities::CommSharedPtr           &vComm,
                const int                                   nDimen)
            {
                NekLinSysIteratSharedPtr p = MemoryManager<
                    NekLinSysIterat>::AllocateSharedPtr(pSession, vComm,nDimen);
                return p;
            }
            LIB_UTILITIES_EXPORT NekLinSysIterat(
                const LibUtilities::SessionReaderSharedPtr  &pSession,
                const LibUtilities::CommSharedPtr           &vComm,
                const int                                   nDimen);
            LIB_UTILITIES_EXPORT ~NekLinSysIterat();
            
            LIB_UTILITIES_EXPORT void setUniversalUniqueMap(Array<OneD, int> &map);
            LIB_UTILITIES_EXPORT void setRhsMagnitude(const NekDouble mag)
            {
                m_rhs_magnitude = mag;
            }
        protected:
                        /// Global to universal unique map
            Array<OneD, int>                            m_map;

            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_rhs_magnitude = NekConstants::kNekUnsetDouble;
                        /// maximum iterations
            int                                         m_maxiter;
            /// Tolerance of iterative solver.
            NekDouble                                   m_tolerance;

            int                                         m_totalIterations = 0;
            /// cnt to how many times rhs_magnitude is called 
            NekDouble                                   m_rhs_mag_sm = 0.9; 

            NekDouble                                   m_prec_factor = 1.0;

            void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);
            void setUniversalUniqueMap();
            
            virtual void v_InitObject();
        private:
        };

        static NekLinSysIteratSharedPtr NullNekLinSysIteratSharedPtr;

    }
}
    
#endif
