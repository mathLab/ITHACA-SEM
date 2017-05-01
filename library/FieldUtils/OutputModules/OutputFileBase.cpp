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
        LibUtilities::CommSharedPtr comm = m_f->m_comm;
        int rank                         = comm->GetRank();
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
