///////////////////////////////////////////////////////////////////////////////
//
// File DriverStandard.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/DriverStandard.h>

namespace Nektar
{
    namespace SolverUtils
    {
        string DriverStandard::className = GetDriverFactory().RegisterCreatorFunction("Standard", DriverStandard::create);
        string DriverStandard::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","Standard",0);

        /**
	 *
         */
        DriverStandard::DriverStandard(const LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
        {
        }
    
    
        /**
         *
         */
        DriverStandard:: ~DriverStandard()
        {
        }
    
    
        /**
         *
         */
        void DriverStandard::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
		}
    
    
        void DriverStandard::v_Execute(ostream &out)
        
        {
            m_equ[0]->PrintSummary(out);
            m_equ[0]->DoInitialise();
            m_equ[0]->DoSolve();
            m_equ[0]->Output();
        
            // Evaluate and output computation time and solution accuracy.
            // The specific format of the error output is essential for the
            // regression tests to work.
            // Evaluate L2 Error
            for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
            {
                NekDouble vL2Error = m_equ[0]->L2Error(i,false);
                NekDouble vLinfError = m_equ[0]->LinfError(i);
                if (m_comm->GetRank() == 0)
                {
                    out << "L 2 error (variable " << m_equ[0]->GetVariable(i) << ") : " << vL2Error << endl;
                    out << "L inf error (variable " << m_equ[0]->GetVariable(i) << ") : " << vLinfError << endl;
                }
            }
        }
    }
}

