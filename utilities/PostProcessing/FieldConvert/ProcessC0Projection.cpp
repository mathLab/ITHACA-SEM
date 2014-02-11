////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessC0Projection.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes C0 projection.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessC0Projection.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessC0Projection::className =
		GetModuleFactory().RegisterCreatorFunction(
												   ModuleKey(eProcessModule, "C0Projection"),
												   ProcessC0Projection::create, "Computes C0 projection.");
		
        ProcessC0Projection::ProcessC0Projection(FieldSharedPtr f) : ProcessModule(f)
        {
			//f->m_declareExpansionAsContField = true;
        }
		
        ProcessC0Projection::~ProcessC0Projection()
        {
        }
		
        void ProcessC0Projection::Process(po::variables_map &vm)
        {
            if (m_f->m_verbose)
            {
                cout << "ProcessC0Projection: Projects fiels into C0 space..." << endl;
            }
            
            int i;
            int nfields = m_f->m_fielddef[0]->m_fields.size();
						
			for (i = 0; i < nfields; ++i)
			{
				m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
										m_f->m_exp[i]->UpdatePhys());
			}
			
			for (i = 0; i < nfields; ++i)
			{
				m_f->m_exp[i]->FwdTrans(m_f->m_exp[i]->GetPhys(),
										m_f->m_exp[i]->UpdateCoeffs());
			}
			
        }
    }
}
