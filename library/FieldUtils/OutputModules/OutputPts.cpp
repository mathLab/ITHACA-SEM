////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputPts.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2016 Kilian Lackhove
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
//  Description: pts file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputPts.h"
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputPts::m_className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "pts"), OutputPts::create, "Writes a pts file.");

OutputPts::OutputPts(FieldSharedPtr f) : OutputModule(f)
{
}

OutputPts::~OutputPts()
{
}

void OutputPts::Process(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "OutputPts: Writing file..." << endl;
        }
    }

    fs::path writefile(filename);
    int writepts = 1;
    if (fs::exists(writefile) && (vm.count("forceoutput") == 0))
    {
        LibUtilities::CommSharedPtr comm = m_f->m_comm;
        int rank                         = comm->GetRank();
        writepts = 0; // set to zero for reduce all to be correct.

        if (rank == 0)
        {
            string answer;
            cout << "Did you wish to overwrite " << filename << " (y/n)? ";
            getline(cin, answer);
            if (answer.compare("y") == 0)
            {
                writepts = 1;
            }
            else
            {
                cout << "Not writing file " << filename
                     << " because it already exists" << endl;
            }
        }

        comm->AllReduce(writepts, LibUtilities::ReduceSum);
    }

    if (writepts)
    {
        LibUtilities::PtsIO ptsIO(m_f->m_comm);
        LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;
        if(m_f->m_fieldPts == LibUtilities::NullPtsField)
        {
            Array<OneD, Array<OneD, NekDouble> > tmp(
                m_f->m_exp[0]->GetCoordim(0) +
                m_f->m_fielddef[0]->m_fields.size());

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

            for (int i = 0; i < m_f->m_fielddef[0]->m_fields.size(); ++i)
            {
                tmp[i + m_f->m_exp[0]->GetCoordim(0)] =
                    m_f->m_exp[i]->GetPhys();
            }
            fPts =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                    m_f->m_exp[0]->GetCoordim(0),
                    m_f->m_fielddef[0]->m_fields,
                    tmp);
        }
        ptsIO.Write(filename, fPts);
    }
}
}
}
