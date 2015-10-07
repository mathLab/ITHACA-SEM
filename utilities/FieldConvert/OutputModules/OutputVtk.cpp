////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.cpp
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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputVtk.h"

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputVtk::m_className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "vtu"), OutputVtk::create,
        "Writes a VTU file.");

OutputVtk::OutputVtk(FieldSharedPtr f) : OutputModule(f)
{
    m_requireEquiSpaced = true;
}

OutputVtk::~OutputVtk()
{
}

void OutputVtk::Process(po::variables_map &vm)
{
    if(!m_f->m_exp.size()) // do nothing if no expansion defined
    {
        return;
    }

    int i, j;
    if (m_f->m_verbose)
    {
        cout << "OutputVtk: Writing file..." << endl;
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    // amend for parallel output if required
    if(m_f->m_session->GetComm()->GetSize() != 1)
    {
        int    dot  = filename.find_last_of('.');
        string ext = filename.substr(dot,filename.length()-dot);
        string procId = "_P" + boost::lexical_cast<std::string>(
            m_f->m_session->GetComm()->GetRank());
        string start = filename.substr(0,dot);
        filename = start + procId + ext;
    }

    // Write solution.
    ofstream outfile(filename.c_str());
    m_f->m_exp[0]->WriteVtkHeader(outfile);

    int nfields, nstrips;
    if (m_f->m_fielddef.size() == 0)
    {
        nfields = 0;
    }
    else
    {
        nfields = m_f->m_fielddef[0]->m_fields.size();
    }
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    // Homogeneous strip variant
    for(int s = 0; s < nstrips; ++s)
    {
        // For each field write out field data for each expansion.
        for (i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            m_f->m_exp[0]->WriteVtkPieceHeader(outfile,i,s);

            // For this expansion write out each field.
            for (j = 0; j < nfields; ++j)
            {
                m_f->m_exp[s*nfields+j]->WriteVtkPieceData(
                    outfile, i, m_f->m_fielddef[0]->m_fields[j]);
            }
            m_f->m_exp[0]->WriteVtkPieceFooter(outfile, i);
        }
    }
    m_f->m_exp[0]->WriteVtkFooter(outfile);
    cout << "Written file: " << filename << endl;
}

}
}
