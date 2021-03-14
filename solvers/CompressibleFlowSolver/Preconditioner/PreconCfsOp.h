///////////////////////////////////////////////////////////////////////////////
//
// File PreconCfsOp.h
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
// Description: PreconCfsOp header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSOP
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSOP

#include <CompressibleFlowSolver/Preconditioner/PreconCfs.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
// =====================================================================
// ==== DEFINITION OF THE CLASS  PreconCfsOps
// ==== Defines operators needed by operator based preconditioner in the
// ==== CFS solver
// =====================================================================
class NekPreconCfsOperators
{
public:
    typedef const Array<OneD, NekDouble> InArrayType;
    typedef Array<OneD, NekDouble> OutArrayType;

    typedef std::function<void(
        const Array<OneD, const Array<OneD, NekDouble>> &,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &, SNekBlkMatSharedPtr &,
        Array<OneD, SNekBlkMatSharedPtr> &, Array<OneD, SNekBlkMatSharedPtr> &,
        Array<OneD, Array<OneD, NekSingle>> &, TensorOfArray4D<NekSingle> &,
        TensorOfArray4D<NekSingle> &, TensorOfArray5D<NekSingle> &,
        TensorOfArray4D<NekSingle> &, TensorOfArray5D<NekSingle> &)>
        FunctorType1;

    typedef std::function<void(InArrayType &, InArrayType &, OutArrayType &,
                               const bool &)>
        FunctorType2;
    typedef Array<OneD, FunctorType1> FunctorType1Array;
    typedef Array<OneD, FunctorType2> FunctorType2Array;
    static const int nfunctor1 = 1;
    static const int nfunctor2 = 0;

    NekPreconCfsOperators(void) : m_functors1(nfunctor1), m_functors2(nfunctor2)
    {
    }
    NekPreconCfsOperators(const NekPreconCfsOperators &in)
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

    NekPreconCfsOperators &operator=(const NekPreconCfsOperators &in)
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
    void DefineCalcPreconMatBRJCoeff(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[0] = std::bind(
            func, obj, std::placeholders::_1, std::placeholders::_2,
            std::placeholders::_3, std::placeholders::_4, std::placeholders::_5,
            std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
            std::placeholders::_9, std::placeholders::_10,
            std::placeholders::_11);
    }

    inline void DoCalcPreconMatBRJCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &gmtxarray,
        SNekBlkMatSharedPtr &gmtVar, Array<OneD, SNekBlkMatSharedPtr> &TraceJac,
        Array<OneD, SNekBlkMatSharedPtr> &TraceJacDeriv,
        Array<OneD, Array<OneD, NekSingle>> &TraceJacDerivSign,
        TensorOfArray4D<NekSingle> &TraceJacArray,
        TensorOfArray4D<NekSingle> &TraceJacDerivArray,
        TensorOfArray5D<NekSingle> &TraceIPSymJacArray,
        TensorOfArray4D<NekSingle> &StdMatDataDBB,
        TensorOfArray5D<NekSingle> &StdMatDataDBDB)
    {
        ASSERTL1(m_functors1[0], "DoNekSysResEval should be defined");
        m_functors1[0](inarray, gmtxarray, gmtVar, TraceJac, TraceJacDeriv,
                       TraceJacDerivSign, TraceJacArray, TraceJacDerivArray,
                       TraceIPSymJacArray, StdMatDataDBB, StdMatDataDBDB);
    }

protected:
    /* Defines three operators
        DoNekSysResEval   :
            evaluations the residual of the Nonlinear/Linear system
            ie. the residual b-Ax and N(x) for linear and
            nonlinear systems, respectively
            May not be used for linear system.
        DoNekSysLhsEval   :
            evaluations the LHS of the Nonlinear/Linear system (Ax),
            where A is the matrix and x is solution vector.
            For linear system A is the coefficient matrix;
            For nonlinear system A is the coefficient matrix in
            each nonlinear iterations, for example A is the
            Jacobian matrix for Newton method;
        DoNekSysPrecon      :
            Preconditioning operator of the system.
        DoNekSysFixPointIte  :
            Operator to calculate RHS of fixed point iterations
            (x^{n+1}=M^{-1}(b-N*x^{n}), with M+N=A).
    */
    FunctorType1Array m_functors1;
    FunctorType2Array m_functors2;
};

//  Forward declaration
class PreconCfsOp;

typedef std::shared_ptr<PreconCfsOp> PreconCfsOpSharedPtr;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<
    std::string, PreconCfsOp,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &,
    const LibUtilities::SessionReaderSharedPtr &,
    const LibUtilities::CommSharedPtr &>
    PreconCfsOpFactory;

/// Declaration of the boundary condition factory singleton
PreconCfsOpFactory &GetPreconCfsOpFactory();

/**
 * High level abstraction of operator based preconditioner.
 * In some cases, the operators is not presented in matrix form
 * (matrix-free implementation for example), thus operator classes needs to
 * be provided to perform the operators needed in preconditioning.
 *
 * A brother class, PreconCfsMat for example, which stores the system
 * coefficient matrix, could be designed. With the matrix, preconditioning
 * matrix could be formed directly within this class
 * (for example the ILU factorization).
 * The operator defined in this class is not necessary in PreconCfsMat.
 *
 */
class PreconCfsOp : public PreconCfs
{
public:
    PreconCfsOp(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::CommSharedPtr &vComm);
    virtual ~PreconCfsOp()
    {
    }

    inline void SetOperators(const NekPreconCfsOperators &in)
    {
        m_operator = in;
    }

protected:
    NekPreconCfsOperators m_operator;

    virtual void v_InitObject();

private:
    void NullPreconCfsOp(void);

    virtual void v_DoPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
        const bool &flag);

    virtual void v_BuildPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, const Array<OneD, NekDouble>> &intmp,
        const NekDouble time, const NekDouble lambda);

    static std::string lookupIds[];
    static std::string def;
};
typedef std::shared_ptr<PreconCfsOp> PreconCfsOpSharedPtr;
} // namespace Nektar

#endif
