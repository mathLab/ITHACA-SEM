///////////////////////////////////////////////////////////////////////////////
//
// File  NekSys.h
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
// Description: NekSys header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NekSys_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NekSys_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <iomanip>
namespace Nektar
{
namespace LibUtilities
{

// =====================================================================
// ==== DEFINITION OF THE CLASS  NekSysOperators
// ==== Defines operators needed by iterative solver of nonlinear or
// ==== linear system
// =====================================================================
class NekSysOperators
{
public:
    typedef const Array<OneD, NekDouble> InArrayType;
    typedef Array<OneD, NekDouble> OutArrayType;

    typedef std::function<void(InArrayType &, OutArrayType &, const bool &)>
        FunctorType1;
    typedef std::function<void(InArrayType &, InArrayType &, OutArrayType &,
                               const bool &)>
        FunctorType2;
    typedef Array<OneD, FunctorType1> FunctorType1Array;
    typedef Array<OneD, FunctorType2> FunctorType2Array;
    static const int nfunctor1 = 3;
    static const int nfunctor2 = 1;

    NekSysOperators(void) : m_functors1(nfunctor1), m_functors2(nfunctor2)
    {
    }
    NekSysOperators(const NekSysOperators &in)
        : m_functors1(nfunctor1), m_functors2(nfunctor2)
    {
        for (int i = 0; i < nfunctor1; ++i)
        {
            m_functors1[i] = in.m_functors1[i];
        }
        for (int i = 0; i < nfunctor2; ++i)
        {
            m_functors2[i] = in.m_functors2[i];
        }
    }

    NekSysOperators &operator=(const NekSysOperators &in)
    {
        for (int i = 0; i < nfunctor1; ++i)
        {
            m_functors1[i] = in.m_functors1[i];
        }
        for (int i = 0; i < nfunctor2; ++i)
        {
            m_functors2[i] = in.m_functors2[i];
        }

        return *this;
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineNekSysRhsEval(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[0] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }
    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineNekSysLhsEval(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[1] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }
    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineNekSysPrecond(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[2] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }
    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineNekSysFixPointIte(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors2[0] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4);
    }

    inline void DoNekSysRhsEval(InArrayType &inarray, OutArrayType &outarray,
                                const bool &flag = false) const
    {
        ASSERTL1(m_functors1[0], "DoNekSysRhsEval should be defined");
        m_functors1[0](inarray, outarray, flag);
    }

    inline void DoNekSysLhsEval(InArrayType &inarray, OutArrayType &outarray,
                                const bool &flag = false) const
    {
        ASSERTL1(m_functors1[1], "DoNekSysLhsEval should be defined");
        m_functors1[1](inarray, outarray, flag);
    }

    inline void DoNekSysPrecond(InArrayType &inarray, OutArrayType &outarray,
                                const bool &flag = false) const
    {
        if (m_functors1[2])
        {
            m_functors1[2](inarray, outarray, flag);
        }
        else
        {
            Vmath::Vcopy(outarray.size(), inarray, 1, outarray, 1);
        }
    }

    inline void DoNekSysFixPointIte(InArrayType &rhs, InArrayType &xn,
                                    OutArrayType &xn1,
                                    const bool &flag = false) const
    {
        ASSERTL1(m_functors2[0], "DoNekSysFixPointIte should be defined");
        m_functors2[0](rhs, xn, xn1, flag);
    }

protected:
    /* Defines three operators
        DoNekSysRhsEval   :
            evaluations the RHS of the Nonlinear/Linear system.
            May not be used for linear system.
        DoNekSysLhsEval   :
            evaluations the LHS of the Nonlinear/Linear system (Ax),
            where A is the matrix x is solution vector.
            For linear system A is the coefficient matrix;
            For nonlinear system A is the coefficient matrix in
            each nonlinear iterations, for example A is the
            Jacobian matrix for Newton method;
        DoNekSysPrecond      :
            Preconditioning operator of the system.
        DoNekSysFixPointIte  :
            Operator to calculate RHS of fix point iterations
            (x^{n+1}=M^{-1}(b-N*x^{n}), with M+N=A).
    */
    FunctorType1Array m_functors1;
    FunctorType2Array m_functors2;
};

class NekSys;

typedef std::shared_ptr<NekSys> NekSysSharedPtr;

class NekSys : public std::enable_shared_from_this<NekSys>
{
public:
    /// Support creation through MemoryManager.
    friend class MemoryManager<NekSys>;

    LIB_UTILITIES_EXPORT static NekSysSharedPtr CreateInstance(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen)
    {
        NekSysSharedPtr p =
            MemoryManager<NekSys>::AllocateSharedPtr(pSession, vComm, nDimen);
        return p;
    }
    LIB_UTILITIES_EXPORT NekSys(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen);
    LIB_UTILITIES_EXPORT void InitObject()
    {
        v_InitObject();
    }
    LIB_UTILITIES_EXPORT virtual ~NekSys();

    LIB_UTILITIES_EXPORT inline void setSysOperators(const NekSysOperators &in)
    {
        m_operator = in;
    }

    LIB_UTILITIES_EXPORT int SolveSystem(
        const int nGlobal, const Array<OneD, const NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput, const int nDir,
        const NekDouble tol = 1.0E-7, const NekDouble factor = 1.0)
    {
        return v_SolveSystem(nGlobal, pInput, pOutput, nDir, tol, factor);
    }

    LIB_UTILITIES_EXPORT bool ConvergenceCheck(
        const int nIteration, const Array<OneD, const NekDouble> &Residual,
        const NekDouble tol = 1.0E-7)
    {
        return v_ConvergenceCheck(nIteration, Residual, tol);
    }

protected:
    /// Maximum iterations
    int m_maxiter;
    /// Tolerance of iterative solver.
    NekDouble m_tolerance;
    /// Communicate.
    LibUtilities::CommSharedPtr m_Comm;
    /// Whether the iteration has been converged
    bool m_converged;
    /// Root if parallel
    bool m_root;
    /// Verbose
    bool m_verbose;
    /// Operators
    NekSysOperators m_operator;
    /// The dimension of the system
    int m_SysDimen;

    virtual void v_InitObject()
    {
    }

    virtual int v_SolveSystem(const int nGlobal,
                              const Array<OneD, const NekDouble> &pInput,
                              Array<OneD, NekDouble> &pOutput, const int nDir,
                              const NekDouble tol, const NekDouble factor)
    {
        boost::ignore_unused(nGlobal, pInput, pOutput, nDir, tol, factor);
        ASSERTL0(false, "LinSysIterSovler NOT CORRECT.");
        return 0;
    }

    virtual bool v_ConvergenceCheck(
        const int nIteration, const Array<OneD, const NekDouble> &Residual,
        const NekDouble tol);
};
} // namespace LibUtilities
} // namespace Nektar
#endif
