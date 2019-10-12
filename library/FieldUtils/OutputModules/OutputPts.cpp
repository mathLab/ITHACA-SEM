////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputPts.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2017 Kilian Lackhove
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: pts file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>

#include "OutputPts.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputPts::m_className[5] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "pts"), OutputPts::create, "Writes a pts file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "csv"), OutputPts::create, "Writes a csv file."),
};


OutputPts::OutputPts(FieldSharedPtr f) : OutputFileBase(f)
{
}

OutputPts::~OutputPts()
{
}

void OutputPts::OutputFromPts(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    if (boost::filesystem::path(filename).extension() == ".csv")
    {
        LibUtilities::CsvIO csvIO(m_f->m_comm);
        csvIO.Write(filename, m_f->m_fieldPts);
    }
    else
    {
        LibUtilities::PtsIO ptsIO(m_f->m_comm);
        ptsIO.Write(filename, m_f->m_fieldPts);
    }
}

void OutputPts::OutputFromExp(po::variables_map &vm)
{
    Array<OneD, Array<OneD, NekDouble> > tmp(
        m_f->m_exp[0]->GetCoordim(0) +
        m_f->m_variables.size());

    switch (m_f->m_exp[0]->GetCoordim(0))
    {
        case 1:
            tmp[0] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            m_f->m_exp[0]->GetCoords(tmp[0]);
            break;

        case 2:
            tmp[1] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            tmp[0] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            m_f->m_exp[0]->GetCoords(tmp[0], tmp[1]);
            break;

        case 3:
            tmp[2] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            tmp[1] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            tmp[0] = Array<OneD, NekDouble>(m_f->m_exp[0]->GetTotPoints());
            m_f->m_exp[0]->GetCoords(tmp[0], tmp[1], tmp[2]);
            break;
    }

    for (int i = 0; i < m_f->m_variables.size(); ++i)
    {
        tmp[i + m_f->m_exp[0]->GetCoordim(0)] =
            m_f->m_exp[i]->GetPhys();
    }
    m_f->m_fieldPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            m_f->m_exp[0]->GetCoordim(0),
            m_f->m_variables,
            tmp);

    OutputFromPts(vm);
}

void OutputPts::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal,
             "OutputPts can't write using only FieldData.");
}

fs::path OutputPts::GetPath(std::string &filename,
                            po::variables_map &vm)
{
    boost::ignore_unused(vm);
    return   fs::path(filename);
}

fs::path OutputPts::GetFullOutName(std::string &filename,
                                po::variables_map &vm)
{
    boost::ignore_unused(vm);
    return   fs::path(filename);
}

}
}

