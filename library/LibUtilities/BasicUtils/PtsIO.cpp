////////////////////////////////////////////////////////////////////////////////
//
// File: PtsIO.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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
// Description: I/O routines relating to Pts data
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/PtsIO.h>

#include <LibUtilities/BasicUtils/FileSystem.h>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>

#include <sstream>
#include <iostream>
#include <fstream>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{

void Import(const string &inFile, PtsFieldSharedPtr &ptsField)
{
#ifdef NEKTAR_USE_MPI
    int size;
    int init;
    MPI_Initialized(&init);

    // If MPI has been initialised we can check the number of processes
    // and, if > 1, tell the user he should not be running this
    // function in parallel. If it is not initialised, we do not
    // initialise it here, and assume the user knows what they are
    // doing.
    if (init)
    {
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        ASSERTL0(size == 1,
                "This static function is not available in parallel. Please"
                "instantiate a FieldIO object for parallel use.");
    }
#endif
    CommSharedPtr c = GetCommFactory().CreateInstance("Serial", 0, 0);
    PtsIO p(c);
    p.Import(inFile, ptsField);
}


void Write(const string &outFile, const PtsFieldSharedPtr &ptsField)
{
#ifdef NEKTAR_USE_MPI
    int size;
    int init;
    MPI_Initialized(&init);

    // If MPI has been initialised we can check the number of processes
    // and, if > 1, tell the user he should not be running this
    // function in parallel. If it is not initialised, we do not
    // initialise it here, and assume the user knows what they are
    // doing.
    if (init)
    {
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        ASSERTL0(size == 1,
                "This static function is not available in parallel. Please"
                "instantiate a FieldIO object for parallel use.");
    }
#endif
    CommSharedPtr c = GetCommFactory().CreateInstance("Serial", 0, 0);
    PtsIO p(c);
    p.Write(outFile, ptsField);
}



/**
 * @brief Import a pts field from file
 *
 * @param inFile    filename of the file to read
 * @param ptsField  the resulting pts field.
 */
void PtsIO::Import(const string &inFile, PtsFieldSharedPtr &ptsField)
{
    TiXmlDocument docInput;
    ASSERTL0(docInput.LoadFile(inFile), "Unable to open file '" + inFile + "'.");

    TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");
    TiXmlElement *points = nektar->FirstChildElement("POINTS");
    int dim;
    int err = points->QueryIntAttribute("DIM", &dim);

    ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute DIM.");

    std::string fields = points->Attribute("FIELDS");

    vector<string> fieldNames;
    if (!fields.empty())
    {
        bool valid = ParseUtils::GenerateOrderedStringVector(
                         fields.c_str(), fieldNames);
        ASSERTL0(valid,
                 "Unable to process list of field variable in  FIELDS attribute:  " + fields);
    }

    int nfields = fieldNames.size();
    int totvars = dim + nfields;

    TiXmlNode *pointsBody = points->FirstChild();

    std::istringstream pointsDataStrm(pointsBody->ToText()->Value());

    vector<NekDouble> ptsSerial;
    Array<OneD,  Array<OneD,  NekDouble> > pts(totvars);

    try
    {
        NekDouble      ptsStream;
        while (!pointsDataStrm.fail())
        {
            pointsDataStrm >> ptsStream;

            ptsSerial.push_back(ptsStream);
        }
    }
    catch (...)
    {
        ASSERTL0(false, "Unable to read Points data.");
    }

    int npts = ptsSerial.size() / totvars;

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


/**
 * @brief Save a pts field to a file
 *
 * @param outFile    filename of the file
 * @param ptsField  the pts field
 */
void PtsIO::Write(const string &outFile,
                  const Nektar::LibUtilities::PtsFieldSharedPtr &ptsField)
{
    int nTotvars = ptsField->GetNFields() + ptsField->GetDim();
    int np = ptsField->GetNpoints();

    std::string filename = SetUpOutput(outFile);

    // until tinyxml gains support for line break, write the xml manually
    std::ofstream ptsFile;
    ptsFile.open(filename.c_str());

    ptsFile << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
    ptsFile << "<NEKTAR>" << endl;
    ptsFile << "  <POINTS ";
    ptsFile << "DIM=\"" << ptsField->GetDim() << "\" ";
    string fn = boost::algorithm::join(ptsField->GetFieldNames(), ",");
    ptsFile << "FIELDS=\"" << fn << "\" ";
    ptsFile << ">" << endl;

    Array <OneD, Array <OneD, NekDouble > > pts;
    ptsField->GetPts(pts);
    for (int i = 0; i < np; ++i)
    {
        ptsFile << "    ";
        ptsFile << pts[0][i];
        for (int j = 1; j < nTotvars; ++j)
        {
            ptsFile << " " << pts[j][i];
        }
        ptsFile << endl;
    }
    ptsFile << "  </POINTS>" << endl;
    ptsFile << "</NEKTAR>" << endl;

    ptsFile.close();

    /*
    // Create the file (partition)
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    TiXmlElement *pointsTag = new TiXmlElement("POINTS");
    root->LinkEndChild(pointsTag);

    pointsTag->SetAttribute("DIM", ptsField->GetDim());

    string fn = boost::algorithm::join(ptsField->GetFieldNames(), ",");
    pointsTag->SetAttribute("FIELDS", fn);

    Array <OneD, Array <OneD, NekDouble > > pts;
    ptsField->GetPts(pts);
    ostringstream os;
    for (int i = 0; i < np; ++i)
    {
        os << pts[0][i];
        for (int j = 1; j < nTotvars; ++j)
        {
            os << " " << pts[j][i];
        }
        os << " ";
    }

    pointsTag->LinkEndChild(new TiXmlText(os.str()));

    doc.SaveFile(filename);
    */
}

std::string PtsIO::SetUpOutput(const std::string outname)
{
    ASSERTL0(!outname.empty(), "Empty path given to SetUpOutput()");

    int nprocs = m_comm->GetSize();

    // Directory name if in parallel, regular filename if in serial
    fs::path specPath(outname);

    // Remove any existing file which is in the way
    if (m_comm->RemoveExistingFiles())
    {
        try
        {
            fs::remove_all(specPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                     "Filesystem error: " + string(e.what()));
        }
    }

    // serial processing just add ending.
    if (nprocs == 1)
    {
        cout << "Writing: " << specPath << endl;
        return LibUtilities::PortablePath(specPath);
    }

    // Create the destination directory
    try
    {
        fs::create_directory(specPath);
    }
    catch (fs::filesystem_error &e)
    {
        ASSERTL0(false, "Filesystem error: " + string(e.what()));
    }

    // Pad rank to 8char filenames, e.g. P0000000.fld
    boost::format pad("P%1$07d.fld");
    pad % m_comm->GetRank();

    // Generate full path name
    fs::path poutfile(pad.str());
    fs::path fulloutname = specPath / poutfile;

    // Return the full path to the partition for this process
    return LibUtilities::PortablePath(fulloutname);
}


}
}
