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
// #include <SolverUtils/Driver.h>
// #include <SolverUtils/EquationSystem.h>
// #include <LibUtilities/BasicUtils/SessionReader.h>
// #include <SpatialDomains/MeshGraph.h>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekNonlinSys.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
// using namespace Nektar::SolverUtils;

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

            std::string SovlerType = "Newton";
            // if(pSession->DefinesSolverInfo("NonlinIteratSovler"))
            // {
            //     SovlerType = pSession->GetSolverInfo("NonlinIteratSovler");
            // }
            ASSERTL0(LibUtilities::GetNekNonlinSysFactory().
                     ModuleExists(SovlerType),
                     "NekNonlinSys '" + SovlerType + "' is not defined.\n");
            m_nonlinsol = LibUtilities::GetNekNonlinSysFactory().CreateInstance(
                                SovlerType, m_session,m_comm, m_matDim);

            m_LinSysOprtors.DefineNonlinLinSysRhsEval(&LinSysDemo::DoRhs, this);
            m_LinSysOprtors.DefineNonlinLinSysLhsEval(&LinSysDemo::DoLhs, this);
            m_nonlinsol->setSysOperators(m_LinSysOprtors);
        }
        ~LinSysDemo()
        {
        }

        void DoSolve()
        {
            Array<OneD, NekDouble> pOutput(m_matDim, 0.9);

            int ntmpIts =  m_nonlinsol->SolveSystem(m_matDim, pOutput, 
                                                    pOutput, 0, 1.0E-9);

            int ndigits     = 9;  // the number of sigificant digits
            int nothers     = 10; // extra width to place -, E, and power
            int nwidthcolm  = nothers+ndigits - 1; // the second value determines the number of sigificant digits
            cout    << "ntmpIts = " << ntmpIts << endl
                    << std::scientific << std::setw(nwidthcolm) <<
                       std::setprecision(ndigits - 1);

            string vars = "uvwx";
            for (int i = 0;i < m_matDim; ++i)
            {
                cout << "L 2 error (variable " << vars[i] << ") : " <<
                         pOutput[i] << endl;
            }
        }

        void AllocateInitMatrix()
        {
            m_matDim = 2; 
        }

        void DoLhs(
                    InArrayType     &inarray, 
                    OutArrayType    &outarray,
                    const  bool     &flag = false)
        {
            boost::ignore_unused(flag);
            const Array<OneD, const NekDouble> refsol = 
                                               m_nonlinsol->GetRefSolution();

            NekDouble x = refsol[0];
            NekDouble y = refsol[1];

            NekDouble f1 = 3.0 * x * x * inarray[0] + inarray[1];
            NekDouble f2 = 3.0 * y * y * inarray[1] - inarray[0];

            outarray[0] =  f1;
            outarray[1] =  f2;
        }

        void DoRhs(
                    InArrayType     &inarray, 
                    OutArrayType    &outarray,
                    const  bool     &flag = false)
        {
            boost::ignore_unused(flag);
            int ntmp = inarray.size();
            ASSERTL0(m_matDim == ntmp, "m_matDim == ntmp not true");
            NekDouble x = inarray[0];
            NekDouble y = inarray[1];
            NekDouble f1 = x * x * x + y - 1.0;
            NekDouble f2 = -x + y * y * y + 1.0;
            outarray[0]  = f1;
            outarray[1]  = f2;
        }

    protected:
        int                                         m_matDim;
        DNekMatSharedPtr                            m_matrix;
        Array<OneD, NekDouble>                      m_matDat;
        Array<OneD, NekDouble>                      m_SysRhs;
        NekNonlinSysSharedPtr                       m_nonlinsol;
        LibUtilities::NonlinLinSysOperators         m_LinSysOprtors;
        LibUtilities::SessionReaderSharedPtr        m_session;
        LibUtilities::CommSharedPtr                 m_comm;
        Array<OneD, int>                            m_map;
};

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    // SpatialDomains::MeshGraphSharedPtr graph;

    session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    session->InitSession();
    // graph = SpatialDomains::MeshGraph::Read(session);

    LinSysDemo linsys(session, session->GetComm());

    linsys.DoSolve();

    // Finalise communications
    session->Finalise();
    return 0;
}
