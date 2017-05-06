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
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputFld::m_className[2] = {
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "fld"),
                                               OutputFld::create,
                                               "Writes a Fld file."),
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "chk"),
                                               OutputFld::create,
                                               "Writes a Fld file."),
};

OutputFld::OutputFld(FieldSharedPtr f) : OutputFileBase(f)
{
    m_config["format"] = ConfigOption(
        false, "Xml", "Output format of field file");
}

OutputFld::~OutputFld()
{
}

void OutputFld::OutputFromPts(po::variables_map &vm)
{
    //ASSERTL0(false, "OutputFld can't write using Pts information.");
    OutputFromExp(vm);
}

void OutputFld::OutputFromExp(po::variables_map &vm)
{
    ASSERTL0(m_f->m_variables.size(),
            "OutputFld: need input data.")

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    // Set up communicator and FieldIO object.
    LibUtilities::CommSharedPtr c = m_f->m_session ? m_f->m_session->GetComm() :
        LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::GetFieldIOFactory().CreateInstance(
            GetIOFormat(), c, true);

    int i, j, s;
    int nfields = m_f->m_variables.size();
    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    if(m_f->m_exp[0]->GetNumElmts() != 0)
    {
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            m_f->m_exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        for (s = 0; s < nstrips; ++s)
        {
            for (j = 0; j < nfields; ++j)
            {
                for (i = 0; i < FieldDef.size() / nstrips; ++i)
                {
                    int n = s * FieldDef.size() / nstrips + i;

                    FieldDef[n]->m_fields.push_back(m_f->m_variables[j]);
                    m_f->m_exp[s * nfields + j]->AppendFieldData(
                        FieldDef[n], FieldData[n]);
                }
            }
        }
        fld->Write(filename, FieldDef, FieldData, m_f->m_fieldMetaDataMap);
    }
    else
    {
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            std::vector<LibUtilities::FieldDefinitionsSharedPtr>();
        std::vector<std::vector<NekDouble> > FieldData =
            std::vector<std::vector<NekDouble> >();
        fld->Write(filename, FieldDef, FieldData, m_f->m_fieldMetaDataMap);
    }
}

void OutputFld::OutputFromData(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    // Set up communicator and FieldIO object.
    LibUtilities::CommSharedPtr c = m_f->m_session ? m_f->m_session->GetComm() :
        LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::GetFieldIOFactory().CreateInstance(
            GetIOFormat(), c, true);

    fld->Write(filename, m_f->m_fielddef, m_f->m_data,
                   m_f->m_fieldMetaDataMap);
}

std::string OutputFld::GetIOFormat()
{
    std::string iofmt("Xml");
    if(m_f->m_session)
    {
        if (m_f->m_session->DefinesSolverInfo("IOFormat"))
        {
            iofmt = m_f->m_session->GetSolverInfo("IOFormat");
        }
        if (m_f->m_session->DefinesCmdLineArgument("io-format"))
        {
            iofmt =
                m_f->m_session->GetCmdLineArgument<std::string>("io-format");
        }
    }
    if(m_config["format"].m_beenSet)
    {
        iofmt = m_config["format"].as<string>();
    }
    return iofmt;
}

}
}
