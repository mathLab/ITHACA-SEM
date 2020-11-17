///////////////////////////////////////////////////////////////////////////////
//
// File  NekLinSysIterative.h
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
// Description: NekLinSysIterative header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERATIVE_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERATIVE_H

#include <boost/circular_buffer.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{

    // =====================================================================
    // ==== DEFINITION OF THE CLASS  LinSysOperators
    // =====================================================================
    class LinSysOperators
    {
    public:
        // typedef const Array<OneD, const TensorOfArray1D<NekDouble> > InArrayType;
        // typedef       Array<OneD,       Array<OneD, NekDouble> > OutArrayType;

        typedef const TensorOfArray1D<NekDouble> InArrayType;
        typedef       Array<OneD, NekDouble> OutArrayType;
        
        typedef std::function< void (InArrayType&, OutArrayType&, const bool&) >                                   FunctorType1;
        typedef std::function< void (InArrayType&, OutArrayType&, const NekDouble) >                  FunctorType2;
        typedef std::function< void (InArrayType&, OutArrayType&, const NekDouble, const NekDouble) > FunctorType3;
        typedef const FunctorType1& ConstFunctorType1Ref;
        typedef const FunctorType2& ConstFunctorType2Ref;
        typedef const FunctorType3& ConstFunctorType3Ref;
        typedef Array<OneD, FunctorType1> FunctorType1Array;
        typedef Array<OneD, FunctorType2> FunctorType2Array;
        typedef Array<OneD, FunctorType3> FunctorType3Array;
        
        LinSysOperators(void):
        m_functors1(3)
        {
        }
        LinSysOperators(LinSysOperators &in):
        m_functors1(3)
        {
            for (int i = 0; i < 3; i++)
            {
                m_functors1[i] = in.m_functors1[i];
            }
        }
        template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineLinSysRhs(FuncPointerT func, ObjectPointerT obj)
        {
            m_functors1[0] =  std::bind(
                func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }
        template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMatrixMultiply(FuncPointerT func, ObjectPointerT obj)
        {
            m_functors1[1] =  std::bind(
                func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }
        template<typename FuncPointerT, typename ObjectPointerT> 
            void DefinePrecond(FuncPointerT func, ObjectPointerT obj)
        {
            m_functors1[2] =  std::bind(
                func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
        }
        
        inline void DoLinSysRhs(InArrayType     &inarray, 
                                OutArrayType    &outarray,
                                const  bool      &flag = false) const
        {
            ASSERTL1(m_functors1[0],"DoLinSysRhs should be defined for this time integration scheme");
            m_functors1[0](inarray,outarray,flag);
        }
        
        inline void DoMatrixMultiply(InArrayType     &inarray, 
                                     OutArrayType    &outarray,
                                    const  bool      &flag = false) const
        {
            ASSERTL1(m_functors1[1],"DoMatrixMultiply should be defined for this time integration scheme");
            m_functors1[1](inarray,outarray,flag);
        }
        
        inline void DoPrecond(  InArrayType     &inarray, 
                                OutArrayType    &outarray,
                                const  bool      &flag = false) const
        {
            ASSERTL1(m_functors1[2],"DoPrecond should be defined for this time integration scheme");
            m_functors1[2](inarray,outarray,flag);
        }
    protected:
        FunctorType1Array m_functors1;
    private:
    };
 /// A global linear system.
    class  NekLinSysIterative;

    typedef std::shared_ptr<NekLinSysIterative> NekLinSysIterativeSharedPtr;
    
    class  NekLinSysIterative
    {
    public:

        /// Support creation through MemoryManager.
        friend class MemoryManager<NekLinSysIterative>;
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
        LIB_UTILITIES_EXPORT static NekLinSysIterativeSharedPtr CreateInstance(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::CommSharedPtr &vComm)
        {
            NekLinSysIterativeSharedPtr p = MemoryManager<
                NekLinSysIterative>::AllocateSharedPtr(pSession, vComm);
            return p;
        }
        /// Constructor for full direct matrix solve.
         NekLinSysIterative(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::CommSharedPtr &vComm);
         ~NekLinSysIterative();
        
        inline void setLinSysOperators(LinSysOperators &in)
        {
            m_oprtor = LinSysOperators(in);
        }
            
        int SolveLinearSystem(
            const int nGlobal,
            const TensorOfArray1D<NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const int nDir,
            const NekDouble  tol    =   1.0E-7,
            const NekDouble  factor =   1.0);
        
        int GetMaxLinIte()
        {
            return (m_maxrestart*m_maxstorage);
        }
     
 protected:
        /// Global to universal unique map
        Array<OneD, int>                            m_map;
        /// maximum gmres restart iteration
        int                                         m_maxrestart;
        /// maximum gmres search directions for one restart(determines the max storage usage)
        int                                         m_maxstorage;
        /// maximum bandwidth of Hessenburg matrix if use truncted Gmres(m)
        int                                         m_maxhesband;
        /// maximum iterations (for gmres m_maxiter =  m_maxrestart*m_maxdirction)
        // int                                         m_maxiter;
        /// Tolerance of iterative solver.
        NekDouble                                   m_tolerance;
        /// dot product of rhs to normalise stopping criterion
        NekDouble                                   m_rhs_magnitude;
        /// dot product of rhs to normalise stopping criterion
        NekDouble                                   m_prec_factor;
        /// cnt to how many times rhs_magnitude is called
        NekDouble                                   m_rhs_mag_sm;
        // MultiRegions::PreconditionerSharedPtr                     m_precon;
        // MultiRegions::PreconditionerType            m_precontype;
        // MultiRegions::IterativeMethodType           m_IteraterType;
        LibUtilities::CommSharedPtr                 m_Comm;
        
        int                                         m_totalIterations;
        /// Whether to apply projection technique
        bool                                        m_useProjection;
        /// Whether the iteration has been converged
        bool                                        m_converged;
        /// Root if parallel
        bool                                        m_root;
        /// verbose
        bool                                        m_verbose;

        bool                                        m_flag_LeftPrecond   = false;
        bool                                        m_flag_RightPrecond  = true;

        bool                                        m_DifferenceFlag0  = false;
        bool                                        m_DifferenceFlag1  = false;

        /// Storage for solutions to previous linear problems
        boost::circular_buffer<Array<OneD, NekDouble> > m_prevLinSol;
        /// Total counter of previous solutions
        int m_numPrevSols;
        /// A-conjugate projection technique

        LinSysOperators                             m_oprtor;
        
        void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);
 private:
        /// Actual iterative solve-GMRS
        int DoGMRES(
            const int pNumRows,
            const TensorOfArray1D<NekDouble> &pInput,
                  Array<OneD,       NekDouble> &pOutput,
            const int pNumDir);
        void UpdateKnownSolutions(
            const int pGlobalBndDofs,
            const TensorOfArray1D<NekDouble> &pSolution,
            const int pNumDirBndDofs);
        NekDouble CalculateAnorm(
            const int nGlobal,
            const TensorOfArray1D<NekDouble> &in,
            const int nDir);
        
        /// Actual iterative gmres solver for one restart
        NekDouble DoGmresRestart(
            const bool                         restarted,
            const bool                         truncted,
            const int                          nGlobal,
            const TensorOfArray1D<NekDouble> &pInput,
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
            const TensorOfArray1D<NekDouble> &b,
            Array <OneD, NekDouble> &y
        );
        static std::string lookupIds[];
        static std::string def;
    };
}
    
#endif
