////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFileBase.cpp
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
//  Description: Base class for outputting to a file
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFileBase.h"
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <boost/format.hpp>
#include <iomanip>

namespace Nektar
{
namespace FieldUtils
{

OutputFileBase::OutputFileBase(FieldSharedPtr f) : OutputModule(f)
{
    m_requireEquiSpaced = false;
}

OutputFileBase::~OutputFileBase()
{
}

void OutputFileBase::Process(po::variables_map &vm)
{
    string filename = m_config["outfile"].as<string>();
    fs::path outFile(filename);
    int writeFile = 1;
    if (fs::exists(outFile) && (vm.count("forceoutput") == 0))
    {
        LibUtilities::CommSharedPtr comm;
        int rank = 0;
        if (m_f->m_session)
        {
            comm = m_f->m_session->GetComm();
            rank = comm->GetRank();
        }
        else
        {
            comm = LibUtilities::GetCommFactory().CreateInstance(
                "Serial", 0, 0);
        }

        writeFile = 0; // set to zero for reduce all to be correct.

        if (rank == 0)
        {
            string answer;
            cout << "Did you wish to overwrite " << filename << " (y/n)? ";
            getline(cin, answer);
            if (answer.compare("y") == 0)
            {
                writeFile = 1;
            }
            else
            {
                cout << "Not writing file " << filename
                     << " because it already exists" << endl;
            }
        }

        comm->AllReduce(writeFile, LibUtilities::ReduceSum);
    }
    if(writeFile)
    {
        if(m_f->m_fieldPts != LibUtilities::NullPtsField)
        {
            OutputFromPts(vm);
        }
        else if(m_f->m_exp.size())
        {
            // reset expansion definition to use equispaced points if required.
            if (m_requireEquiSpaced && (vm.count("noequispaced") == 0 ) )
            {
                // Information to create new expansion
                int numFields   = m_f->m_exp.size();
                m_f->m_fielddef = m_f->m_exp[0]->GetFieldDefinitions();

                // Set points to equispaced
                int nPointsNew  = 0;
                if (vm.count("output-points"))
                {
                    nPointsNew = vm["output-points"].as<int>();
                }
                m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);

                // Save original expansion
                vector<MultiRegions::ExpListSharedPtr> expOld = m_f->m_exp;
                // Create new expansion
                m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_numHomogeneousDir,
                                                true);
                for(int i = 1; i < numFields; ++i)
                {
                    m_f->m_exp[i] =
                            m_f->AppendExpList(m_f->m_numHomogeneousDir);
                }
                // Extract result to new expansion
                for(int i = 0; i < numFields; ++i)
                {
                    m_f->m_exp[i]->ExtractCoeffsToCoeffs(
                            expOld[i],
                            expOld[i]->GetCoeffs(),
                            m_f->m_exp[i]->UpdateCoeffs());
                    m_f->m_exp[i]->BwdTrans(
                            m_f->m_exp[i]->GetCoeffs(),
                            m_f->m_exp[i]->UpdatePhys());
                }
            }

            OutputFromExp(vm);
        }
        else if(m_f->m_data.size())
        {
            OutputFromData(vm);
        }
    }
}
}
}
