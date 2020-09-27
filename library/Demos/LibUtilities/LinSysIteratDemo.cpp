///////////////////////////////////////////////////////////////////////////////
//
// File: LinSysDemo.cpp
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
// Description: Demo for testing functionality of NodalUtil classes
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSysIterat.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

class LinSysDemo
{
    typedef const Array<OneD, NekDouble> InArrayType;
    typedef       Array<OneD, NekDouble> OutArrayType;
    public:
        LinSysDemo(
            const LibUtilities::SessionReaderSharedPtr        &pSession,
            const LibUtilities::CommSharedPtr                 &pComm     )
            : m_session(pSession),
              m_comm(pComm)
        {
            AllocateInitMatrix();

            std::string LinIteratSovlerType = "FixedpointJacobi";
            if (pSession->DefinesSolverInfo("LinIteratSovler"))
            {
              LinIteratSovlerType = pSession->GetSolverInfo("LinIteratSovler");
            }
            
            ASSERTL0(LibUtilities::GetNekLinSysIteratFactory().
                ModuleExists(LinIteratSovlerType), "NekLinSysIterat '" + 
                LinIteratSovlerType + "' is not defined.\n");
            m_linsol = LibUtilities::GetNekLinSysIteratFactory().CreateInstance(
                            LinIteratSovlerType, m_session, m_comm, m_matDim);

            m_LinSysOprtors.DefineNonlinLinSysLhsEval(&LinSysDemo::DoLhs, this);
            m_LinSysOprtors.DefineNonlinLinFixPointIte(&LinSysDemo::
                DoFixedPoint, this);
            m_linsol->setSysOperators(m_LinSysOprtors);
            UniqueMap();
            m_linsol->setUniversalUniqueMap(m_map);
        }
        ~LinSysDemo()
        {
        }

        void DoSolve()
        {
            Array<OneD, NekDouble> pOutput(m_matDim, 0.0);

            int ntmpIts = m_linsol->SolveSystem(m_matDim, m_SysRhs,
                                                pOutput, 0,1.0E-9);
            boost::ignore_unused(ntmpIts);
            // The number of sigificant digits
            int ndigits     = 9;  
            // Extra width to place -, E, and power
            int nothers     = 10; 
            // The second value determines the number of sigificant digits
            int nwidthcolm  = nothers + ndigits - 1;
            cout << " ntmpIts= " << ntmpIts << endl
                 << std::scientific << std::setw(nwidthcolm)
                 << std::setprecision(ndigits - 1);

            string vars = "uvwx";
            for (int i = 0;i < m_matDim; ++i)
            {
                cout << "L 2 error (variable " << vars[i] << ") : " <<
                        pOutput[i] << endl;
            }
        }

        void AllocateInitMatrix()
        {
            m_matDim = 4; 
            m_matrix = MemoryManager<DNekMat>::AllocateSharedPtr
                                            (m_matDim, m_matDim, 0.0);
            m_matDat = m_matrix->GetPtr();
            m_matDat[0]  = 10.0;
            m_matDat[1]  = -1.0;
            m_matDat[2]  =  2.0;
            m_matDat[3]  =  0.0;
            m_matDat[4]  = -1.0;
            m_matDat[5]  = 11.0;
            m_matDat[6]  = -1.0;
            m_matDat[7]  =  3.0;
            m_matDat[8]  =  2.0;
            m_matDat[9]  = -1.0;
            m_matDat[10] = 10.0;
            m_matDat[11] = -1.0;
            m_matDat[12] =  0.0;
            m_matDat[13] =  3.0;
            m_matDat[14] = -1.0;
            m_matDat[15] =  8.0;

            m_SysRhs = Array<OneD, NekDouble>(m_matDim);
            m_SysRhs[0]  =   6.0;
            m_SysRhs[1]  =  25.0;
            m_SysRhs[2]  = -11.0;
            m_SysRhs[3]  =  15.0;
        }

        void DoFixedPoint(
                    InArrayType     &rhs, 
                    InArrayType     &inarray, 
                    OutArrayType    &outarray,
                    const  bool     &flag = false)
        {
            boost::ignore_unused(flag);

            int ntmp = inarray.size();
            ASSERTL0(m_matDim == ntmp, "m_matDim == ntmp not true");
            NekVector<NekDouble> vecInn (m_matDim, inarray, eWrapper);
            NekVector<NekDouble> vecOut (m_matDim, outarray, eWrapper);
            vecOut = (*m_matrix) * vecInn;

            Vmath::Vsub(m_matDim, rhs, 1, outarray, 1, outarray, 1);
            for (int i = 0; i < m_matDim; ++i)
            {
                outarray[i] = outarray[i] / (*m_matrix)(i,i);
            }
            Vmath::Vadd(m_matDim, inarray, 1, outarray, 1, outarray, 1);

        }

        void DoLhs(
                    InArrayType     &inarray, 
                    OutArrayType    &outarray,
                    const  bool     &flag = false)
        {
            boost::ignore_unused(flag);
            int ntmp = inarray.size();
            ASSERTL0(m_matDim == ntmp,"m_matDim == ntmp not true");
            NekVector<NekDouble> vecInn (m_matDim, inarray, eWrapper);
            NekVector<NekDouble> vecOut (m_matDim, outarray, eWrapper);
            vecOut = (*m_matrix) * vecInn;
        }

        void UniqueMap()
        {
            m_map = Array<OneD, int>(m_matDim, 1);
        }

    protected:
        int                                         m_matDim;
        DNekMatSharedPtr                            m_matrix;
        Array<OneD, NekDouble>                      m_matDat;
        Array<OneD, NekDouble>                      m_SysRhs;
        NekLinSysIteratSharedPtr                    m_linsol;
        LibUtilities::NonlinLinSysOperators         m_LinSysOprtors;
        LibUtilities::SessionReaderSharedPtr        m_session;
        LibUtilities::CommSharedPtr                 m_comm;
        Array<OneD, int>                            m_map;
};

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;

    session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    session->InitSession();

    LinSysDemo linsys(session,session->GetComm());

    linsys.DoSolve();

    // Finalise communications
    session->Finalise();
    return 0;
}
