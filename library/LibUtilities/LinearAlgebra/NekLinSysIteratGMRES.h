///////////////////////////////////////////////////////////////////////////////
//
// File  NekLinSysIteratGMRES.h
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
// Description: NekLinSysIteratGMRES header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_GMRES_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_GMRES_H

#include <LibUtilities/LinearAlgebra/NekLinSysIterat.h>
namespace Nektar
{
    namespace LibUtilities
    {   
        /// A global linear system.
        class  NekLinSysIteratGMRES;

        typedef std::shared_ptr<NekLinSysIteratGMRES> NekLinSysIteratGMRESSharedPtr;
        
        class  NekLinSysIteratGMRES : public NekLinSysIterat
        {
        public:

            /// Support creation through MemoryManager.
            friend class MemoryManager<NekLinSysIteratGMRES>;
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

            LIB_UTILITIES_EXPORT static NekLinSysIteratSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr  &pSession,
                const LibUtilities::CommSharedPtr           &vComm,
                const int                                   nDimen)
            {
                NekLinSysIteratGMRESSharedPtr p = MemoryManager<
                    NekLinSysIteratGMRES>::AllocateSharedPtr(pSession, vComm, nDimen);
                p->InitObject();
                return p;
            }
            static std::string className;
            /// Constructor for full direct matrix solve.
            LIB_UTILITIES_EXPORT NekLinSysIteratGMRES(
                const LibUtilities::SessionReaderSharedPtr  &pSession,
                const LibUtilities::CommSharedPtr           &vComm,
                const int                                   nDimen);
            LIB_UTILITIES_EXPORT ~NekLinSysIteratGMRES();
            
            LIB_UTILITIES_EXPORT int GetMaxLinIte()
            {
                return (m_maxrestart*m_maxstorage);
            }
        
        protected:
            /// maximum gmres restart iteration
            int                                         m_maxrestart;
            /// maximum gmres search directions for one restart(determines the max storage usage)
            int                                         m_maxstorage;
            /// maximum bandwidth of Hessenburg matrix if use truncted Gmres(m)
            int                                         m_maxhesband;
            int                                         m_totalIterations;
            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_prec_factor;
            /// cnt to how many times rhs_magnitude is called
            NekDouble                                   m_rhs_mag_sm;
            /// Whether to apply projection technique
            bool                                        m_useProjection;

            bool                                        m_flag_LeftPrecond   = false;
            bool                                        m_flag_RightPrecond  = true;

            bool                                        m_DifferenceFlag0  = false;
            bool                                        m_DifferenceFlag1  = false;

            void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);

            virtual void v_InitObject();

            virtual int v_SolveSystem(
                const int                           nGlobal,
                const Array<OneD, const NekDouble>  &pInput,
                Array<OneD,      NekDouble>         &pOutput,
                const int                           nDir,
                const NekDouble                     tol    ,
                const NekDouble                     factor );
            
            virtual bool v_ConvergenceCheck(
                const int                           nIteration,
                const Array<OneD, const NekDouble>  &Residual,
                const NekDouble                     tol         );
        private:
            
            /// Actual iterative solve-GMRS
            int DoGMRES(
                const int pNumRows,
                const Array<OneD, const NekDouble> &pInput,
                    Array<OneD,       NekDouble> &pOutput,
                const int pNumDir);
            /// Actual iterative gmres solver for one restart
            NekDouble DoGmresRestart(
                const bool                         restarted,
                const bool                         truncted,
                const int                          nGlobal,
                const Array<OneD, const NekDouble> &pInput,
                Array<OneD,      NekDouble> &pOutput,
                const int                          nDir);
            
            // Arnoldi process
            void DoArnoldi(
                const int starttem,
                const int endtem,
                const int nGlobal,
                const int nDir,
                // V_total(:,1:nd) total search directions
                Array<OneD, Array<OneD,  NekDouble> > &V_local,
                // V[nd] current search direction
                Array<OneD,  NekDouble> &Vsingle1,
                // V[nd+1] new search direction
                Array<OneD, NekDouble> &Vsingle2,
                // One line of Hessenburg matrix
                Array<OneD, NekDouble> &hsingle
            );
            // QR fatorization through Givens rotation
            void DoGivensRotation(
                const int starttem,
                const int endtem,
                const int nGlobal,
                const int nDir,
                Array<OneD, NekDouble> &c,
                Array<OneD, NekDouble> &s,
                Array<OneD, NekDouble> &hsingle,
                Array<OneD, NekDouble> &eta
            );
            // Backward calculation to calculate coeficients of least square problem
            // To notice, Hessenburg's columnns and rows are reverse
            void DoBackward(
                const int  number,
                Array<OneD, Array<OneD, NekDouble> > &A,
                const Array<OneD, const NekDouble> &b,
                Array <OneD, NekDouble> &y
            );
            static std::string lookupIds[];
            static std::string def;
        };
    }
}
    
#endif
