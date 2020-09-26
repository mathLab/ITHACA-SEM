////////////////////////////////////////////////////////////////////////////////
//
// File: CsvIO.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Csv IO
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/CsvIO.h>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <vector>
#include <fstream>

#include <boost/format.hpp>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#include "ErrorUtil.hpp"
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace LibUtilities
{


CsvIO::CsvIO(CommSharedPtr pComm, bool sharedFilesystem)
    : PtsIO(pComm, sharedFilesystem)
{
}

/**
 * @brief Save a pts field to a file
 *
 * @param outFile    filename of the file
 * @param ptsField  the pts field
 */
void CsvIO::Write(const std::string &outFile,
                  const Nektar::LibUtilities::PtsFieldSharedPtr &ptsField,
                  const bool backup)
{
    size_t nTotvars = ptsField->GetNFields() + ptsField->GetDim();
    size_t np = ptsField->GetNpoints();

    std::string filename = SetUpOutput(outFile, true, backup);
    SetUpFieldMetaData(outFile);

    std::ofstream ptsFile;
    ptsFile.open(filename.c_str());

    std::vector<std::string> xyz;
    xyz.push_back("x");
    xyz.push_back("y");
    xyz.push_back("z");
    xyz.resize(ptsField->GetDim());

    std::string fn = boost::algorithm::join(xyz, ",");
    ptsFile << "# " << fn << ",";
    fn = boost::algorithm::join(ptsField->GetFieldNames(), ",");
    ptsFile << fn;
    ptsFile << std::endl;

    Array<OneD, Array<OneD, NekDouble> > pts;
    ptsField->GetPts(pts);
    for (size_t i = 0; i < np; ++i)
    {
        ptsFile << pts[0][i];
        for (size_t j = 1; j < nTotvars; ++j)
        {
            ptsFile << "," << pts[j][i];
        }
        ptsFile << std::endl;
    }

    ptsFile.close();
}


void CsvIO::v_ImportFieldData(const std::string inFile, PtsFieldSharedPtr& ptsField)
{
    std::stringstream errstr;
    errstr << "Unable to load file: " << inFile << std::endl;
    std::ifstream in(inFile.c_str());
    ASSERTL0(in.is_open(), errstr.str());

    std::string line;
    std::getline(in, line);
    boost::erase_first(line, "#");

    std::vector<std::string> fieldNames;
    bool valid = ParseUtils::GenerateVector(line, fieldNames);
    ASSERTL0(valid, "Unable to process list of fields from line: " + line);

    int dim = 0;
    for (auto &it : fieldNames)
    {
        if (it == "x" || it == "y" || it == "z")
        {
            dim++;
        }
    }

    size_t totvars = fieldNames.size();

    std::vector<NekDouble> ptsSerial;
    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;
    Tokenizer tok(line);
    while (getline(in, line))
    {
        tok.assign(line);

        ASSERTL0(std::distance(tok.begin(), tok.end()) ==
                    std::iterator_traits<Tokenizer::iterator>::difference_type(totvars),
                 "wrong number of columns in line: " + line);

        for (auto &it : tok)
        {
            try
            {
                ptsSerial.push_back(
                    boost::lexical_cast<NekDouble>(
                        boost::trim_copy(std::string(it))));
            }
            catch(const boost::bad_lexical_cast &)
            {
                NEKERROR(ErrorUtil::efatal, "could not convert line: " + line);
            }
        }
    }

    size_t npts = ptsSerial.size() / totvars;

    Array<OneD, Array<OneD, NekDouble> > pts(totvars);
    for (size_t i = 0; i < totvars; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts);
    }

    for (size_t i = 0; i < npts; ++i)
    {
        for (size_t j = 0; j < totvars; ++j)
        {
            pts[j][i] = ptsSerial[i * totvars + j];
        }
    }

    // reorder pts to make x,y,z the first columns
    std::vector<std::string> dimNames = {"x", "y", "z"};
    for (int i = 0; i < dim; ++i)
    {
        auto p = std::find(fieldNames.begin(), fieldNames.end(), dimNames[i]);
        if (p != fieldNames.end())
        {
            auto j = std::distance(fieldNames.begin(), p);

            if (i == j)
            {
                continue;
            }

            Array<OneD, NekDouble> tmp = pts[i];
            pts[i] = pts[j];
            pts[j] = tmp;

            std::string tmp2 = fieldNames[i];
            fieldNames[i] = fieldNames[j];
            fieldNames[j] = tmp2;
        }
    }
    fieldNames.erase(fieldNames.begin(), fieldNames.begin() + dim);

    ptsField = MemoryManager<PtsField>::AllocateSharedPtr(dim, fieldNames, pts);
}

}
}
