////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFld.cpp
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
//  Description: FLD file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFld.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputFld::m_className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "fld"), OutputFld::create,
                "Writes a FLD file.");

        OutputFld::OutputFld(FieldSharedPtr f) : OutputModule(f)
        {
        }

        OutputFld::~OutputFld()
        {
        }
        
        void OutputFld::Process(po::variables_map &vm)
        {
            int i, j;
            if (m_f->m_verbose)
            {
                cout << "OutputFld: Writing file..." << endl;
            }
            
            // Extract the output filename and extension
            string filename = m_config["outfile"].as<string>();
            // amend for parallel output if required 
            if(m_f->m_session->GetComm()->GetSize() != 1)
            {
                int    dot  = filename.find_last_of('.');
                string ext = filename.substr(dot,filename.length()-dot);
                string procId = "_P" + boost::lexical_cast<std::string>(m_f->m_session->GetComm()->GetRank());
                string start = filename.substr(0,dot);
                filename = start + procId + ext;
            }
            
            // Write the output file
            LibUtilities::Write(filename, m_f->m_fielddef, m_f->m_data);
        }        
    }
}
