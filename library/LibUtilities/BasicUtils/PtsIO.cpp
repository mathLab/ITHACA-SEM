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

#include <fstream>

#include <boost/format.hpp>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#include <LibUtilities/BasicUtils/FileSystem.h>


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
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERTL0(size == 1,
                 "This static function is not available in parallel. Please "
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
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERTL0(size == 1,
                 "This static function is not available in parallel. Please "
                 "instantiate a FieldIO object for parallel use.");
    }
#endif
    CommSharedPtr c = GetCommFactory().CreateInstance("Serial", 0, 0);
    PtsIO p(c);
    p.Write(outFile, ptsField);
}

PtsIO::PtsIO(CommSharedPtr pComm, bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}

/**
 * @brief Import a pts field from file
 *
 * @param inFile    filename of the file to read
 * @param ptsField  the resulting pts field.
 */
void PtsIO::Import(const string &inFile,
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

        // TODO: This currently only loads the filename matching our rank.
        filenames.clear();
        boost::format pad("P%1$07d.%2$s");
        pad % m_comm->GetRank() % GetFileEnding();
        filenames.push_back(pad.str());

        for (int i = 0; i < filenames.size(); ++i)
        {
            fs::path pfilename(filenames[i]);
            fullpath = pinfilename / pfilename;
            string fname = PortablePath(fullpath);

            TiXmlDocument doc1(fname);
            bool loadOkay1 = doc1.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << fname << std::endl;
            errstr << "Reason: " << doc1.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc1.ErrorRow() << ", Column "
                   << doc1.ErrorCol() << std::endl;
            ASSERTL0(loadOkay1, errstr.str());

            ImportFieldData(doc1, ptsField);
        }
    }
    else
    {
        TiXmlDocument doc(infile);
        bool loadOkay = doc.LoadFile();

        std::stringstream errstr;
        errstr << "Unable to load file: " << infile << std::endl;
        errstr << "Reason: " << doc.ErrorDesc() << std::endl;
        errstr << "Position: Line " << doc.ErrorRow() << ", Column "
               << doc.ErrorCol() << std::endl;
        ASSERTL0(loadOkay, errstr.str());

        ImportFieldData(doc, ptsField);
    }
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
    SetUpFieldMetaData(outFile);

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

    Array<OneD, Array<OneD, NekDouble> > pts;
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

    // this is what the above cpart would read if tinyxml
    // supported line breaks
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

void PtsIO::ImportFieldData(TiXmlDocument docInput, PtsFieldSharedPtr &ptsField)
{
    TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");
    TiXmlElement *points = nektar->FirstChildElement("POINTS");
    int dim;
    int err = points->QueryIntAttribute("DIM", &dim);

    ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute DIM.");

    std::string fields = points->Attribute("FIELDS");

    vector<string> fieldNames;
    if (!fields.empty())
    {
        bool valid =
            ParseUtils::GenerateOrderedStringVector(fields.c_str(), fieldNames);
        ASSERTL0(
            valid,
            "Unable to process list of field variable in  FIELDS attribute:  " +
                fields);
    }

    int nfields = fieldNames.size();
    int totvars = dim + nfields;

    TiXmlNode *pointsBody = points->FirstChild();

    std::istringstream pointsDataStrm(pointsBody->ToText()->Value());

    vector<NekDouble> ptsSerial;
    Array<OneD, Array<OneD, NekDouble> > pts(totvars);

    try
    {
        NekDouble ptsStream;
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

void PtsIO::SetUpFieldMetaData(const string outname)
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
