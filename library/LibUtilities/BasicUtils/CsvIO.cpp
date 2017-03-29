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
// License for the specific language governing rights and limitations under
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

#include <fstream>

#include <boost/format.hpp>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#include "ErrorUtil.hpp"
#include <LibUtilities/BasicUtils/FileSystem.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{


CsvIO::CsvIO(CommSharedPtr pComm, bool sharedFilesystem)
    : FieldIOXml(pComm, sharedFilesystem)
{
}

/**
 * @brief Import a pts field from file
 *
 * @param inFile    filename of the file to read
 * @param ptsField  the resulting pts field.
 */
void CsvIO::Import(const string &inFile,
                   PtsFieldSharedPtr &ptsField,
                   FieldMetaDataMap &fieldmetadatamap)
{
    std::string infile = inFile;

    fs::path pinfilename(infile);

    // check to see that infile is a directory
    if (fs::is_directory(pinfilename))
    {
        fs::path infofile("Info.xml");
        fs::path fullpath = pinfilename / infofile;
        infile = PortablePath(fullpath);

        std::vector<std::string> filenames;
        std::vector<std::vector<unsigned int> > elementIDs_OnPartitions;

        ImportMultiFldFileIDs(
            infile, filenames, elementIDs_OnPartitions, fieldmetadatamap);

        // Load metadata
        ImportFieldMetaData(infile, fieldmetadatamap);

        if (filenames.size() == m_comm->GetSize())
        {
            // only load the file that matches this rank
            filenames.clear();
            boost::format pad("P%1$07d.%2$s");
            pad % m_comm->GetRank() % GetFileEnding();
            filenames.push_back(pad.str());
        }

        for (int i = 0; i < filenames.size(); ++i)
        {
            fs::path pfilename(filenames[i]);
            fullpath = pinfilename / pfilename;
            string fname = PortablePath(fullpath);

            if (i == 0)
            {
                ImportFieldData(fname, ptsField);
            }
            else
            {
                LibUtilities::PtsFieldSharedPtr newPtsField;
                ImportFieldData(fname, newPtsField);
                Array<OneD, Array<OneD, NekDouble> > pts;
                newPtsField->GetPts(pts);
                ptsField->AddPoints(pts);
            }
        }
    }
    else
    {
        ImportFieldData(infile, ptsField);
    }
}

/**
 * @brief Save a pts field to a file
 *
 * @param outFile    filename of the file
 * @param ptsField  the pts field
 */
void CsvIO::Write(const string &outFile,
                  const Nektar::LibUtilities::PtsFieldSharedPtr &ptsField,
                  const bool backup)
{
    int nTotvars = ptsField->GetNFields() + ptsField->GetDim();
    int np = ptsField->GetNpoints();

    std::string filename = SetUpOutput(outFile, true, backup);
    SetUpFieldMetaData(outFile);

    std::ofstream ptsFile;
    ptsFile.open(filename.c_str());

    vector<string> xyz;
    xyz.push_back("x");
    xyz.push_back("y");
    xyz.push_back("z");
    xyz.resize(ptsField->GetDim());

    string fn = boost::algorithm::join(xyz, ",");
    ptsFile << "# " << fn << ",";
    fn = boost::algorithm::join(ptsField->GetFieldNames(), ",");
    ptsFile << fn;
    ptsFile << endl;

    Array<OneD, Array<OneD, NekDouble> > pts;
    ptsField->GetPts(pts);
    for (int i = 0; i < np; ++i)
    {
        ptsFile << pts[0][i];
        for (int j = 1; j < nTotvars; ++j)
        {
            ptsFile << "," << pts[j][i];
        }
        ptsFile << endl;
    }

    ptsFile.close();
}

void CsvIO::ImportFieldData(const std::string inFile, PtsFieldSharedPtr& ptsField)
{
    std::stringstream errstr;
    errstr << "Unable to load file: " << inFile << std::endl;
    ifstream in(inFile.c_str());
    ASSERTL0(in.is_open(), errstr.str());

    string line;
    getline(in, line);
    boost::erase_first(line, "#");

    vector<string> fieldNames;
    bool valid = ParseUtils::GenerateOrderedStringVector(line.c_str(), fieldNames);
    ASSERTL0(valid, "Unable to process list of fields" + line);

    int dim = 0;
    for (vector<string>::iterator it = fieldNames.begin(); it != fieldNames.end(); ++it)
    {
        if (*it == "x" || *it == "y" || *it == "z")
        {
            dim++;
        }
    }
    fieldNames.erase(fieldNames.begin(), fieldNames.begin() + dim);

    int nfields = fieldNames.size();
    int totvars = dim + nfields;

    vector<NekDouble> ptsSerial;
    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;
    Tokenizer tok(line);
    while (getline(in, line))
    {
        tok.assign(line);

        ASSERTL0(distance(tok.begin(), tok.end()) == totvars, "wrong number of columns: " + line);

        for (Tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
        {
            try
            {
                ptsSerial.push_back(boost::lexical_cast<NekDouble>(boost::trim_copy(string(*it))));
            }
            catch(const boost::bad_lexical_cast &)
            {
                ASSERTL0(false, "could not convert line: " + line);
            }
        }
    }

    int npts = ptsSerial.size() / totvars;

    Array<OneD, Array<OneD, NekDouble> > pts(totvars);
    for (int i = 0; i < totvars; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts);
    }

    for (int i = 0; i < npts; ++i)
    {
        for (int j = 0; j < totvars; ++j)
        {
            pts[j][i] = ptsSerial[i * totvars + j];
        }
    }

    ptsField = MemoryManager<PtsField>::AllocateSharedPtr(dim, fieldNames, pts);
}

void CsvIO::SetUpFieldMetaData(const string outname)
{
    ASSERTL0(!outname.empty(), "Empty path given to SetUpFieldMetaData()");

    int nprocs = m_comm->GetSize();
    int rank = m_comm->GetRank();

    fs::path specPath(outname);

    // Collate per-process element lists on root process to generate
    // the info file.
    if (rank == 0)
    {
        // Set up output names
        std::vector<std::string> filenames;
        std::vector<std::vector<unsigned int> > ElementIDs;
        for (int i = 0; i < nprocs; ++i)
        {
            boost::format pad("P%1$07d.%2$s");
            pad % i % GetFileEnding();
            filenames.push_back(pad.str());

            std::vector<unsigned int> tmp;
            tmp.push_back(0);
            ElementIDs.push_back(tmp);
        }

        // Write the Info.xml file
        string infofile =
            LibUtilities::PortablePath(specPath / fs::path("Info.xml"));

        cout << "Writing: " << specPath << endl;

        const FieldMetaDataMap fieldmetadatamap;
        WriteMultiFldFileIDs(infofile, filenames, ElementIDs, fieldmetadatamap);
    }
}
}
}
