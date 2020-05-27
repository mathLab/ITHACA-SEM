///////////////////////////////////////////////////////////////////////////////
//
// File: PartitionAnalyse.cpp
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
// Description: Small utility to export histogram of partition sizes.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshPartition.h>

#include <iostream>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::SpatialDomains;

class FauxComm : public CommSerial
{
public:
    FauxComm(int argc, char* argv[], int size)
        : CommSerial(argc, argv)
    {
        m_size = size;
        m_type = "Faux parallel";
    }
    FauxComm(int size) : CommSerial(0, NULL)
    {
        m_size = size;
        m_type = "Faux parallel";
    }
    virtual ~FauxComm() {}
    void v_SplitComm(int pRows, int pColumns)
    {
        m_commRow    = std::shared_ptr<FauxComm>(new FauxComm(pColumns));
        m_commColumn = std::shared_ptr<FauxComm>(new FauxComm(pRows));
    }
};

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: PartitionAnalyse <nproc> <xml file1> [xml file 2..n]"
             << endl;
        return 1;
    }

    int nParts = atoi(argv[1]);
    vector<string> filenames(argv + 2, argv + argc);
    
    CommSharedPtr vComm = std::shared_ptr<FauxComm>(
        new FauxComm(argc, argv, nParts));

    char **new_argv = new char*[argc];
    new_argv[0] = strdup("PartitionAnalyse");
    new_argv[1] = strdup("--part-info");
    for (int i = 0; i < argc-2; ++i)
    {
        new_argv[i+2] = strdup(filenames[i].c_str());
    }

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(
            argc, new_argv, filenames, vComm);

    return 0;
}

