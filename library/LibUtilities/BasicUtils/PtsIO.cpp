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
#include <LibUtilities/BasicConst/GitRevision.h>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/asio/ip/host_name.hpp>

#include <sstream>
#include <iostream>
#include <fstream>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace ptime = boost::posix_time;
namespace ip = boost::asio::ip;
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
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        ASSERTL0(size == 1,
                "This static function is not available in parallel. Please "
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
void PtsIO::Import(const string &inFile, PtsFieldSharedPtr &ptsField, FieldMetaDataMap &fieldmetadatamap)
{
    std::string infile = inFile;

    fs::path pinfilename(infile);

    if(fs::is_directory(pinfilename)) // check to see that infile is a directory
    {
        fs::path infofile("Info.xml");
        fs::path fullpath = pinfilename / infofile;
        infile = PortablePath(fullpath);

        std::vector<std::string> filenames;
        std::vector<std::vector<unsigned int> > elementIDs_OnPartitions;


        ImportMultiFldFileIDs(infile,filenames, fieldmetadatamap);

        // Load metadata
        ImportFieldMetaData(infile, fieldmetadatamap);

        //HACK: only load the filename matching our rank.
        filenames.clear();
        boost::format pad("P%1$07d.pts");
        pad % m_comm->GetRank();
        filenames.push_back(pad.str());

        for(int i = 0; i < filenames.size(); ++i)
        {
            fs::path pfilename(filenames[i]);
            fullpath = pinfilename / pfilename;
            string fname = PortablePath(fullpath);

            TiXmlDocument doc1(fname);
            bool loadOkay1 = doc1.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << fname << std::endl;
            errstr << "Reason: " << doc1.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc1.ErrorRow() << ", Column " << doc1.ErrorCol() << std::endl;
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
        errstr << "Position: Line " << doc.ErrorRow() << ", Column " <<
            doc.ErrorCol() << std::endl;
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

    cout << "writing to " << filename << endl;

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
    int rank   = m_comm->GetRank();

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
    m_comm->Block();

    // serial processing just add ending.
    if (nprocs == 1)
    {
        cout << "Writing: " << specPath << endl;
        return LibUtilities::PortablePath(specPath);
    }

    // Collate per-process element lists on root process to generate
    // the info file.
    if (rank == 0)
    {
        // Create the destination directory
        try
        {
            fs::create_directory(specPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(false, "Filesystem error: " + string(e.what()));
        }

        // Set up output names
        std::vector<std::string> filenames;
        for (int i = 0; i < nprocs; ++i)
        {
            boost::format pad("P%1$07d.pts");
            pad % i;
            filenames.push_back(pad.str());
        }

        // Write the Info.xml file
        string infofile = LibUtilities::PortablePath(
                              specPath / fs::path("Info.xml"));

        cout << "Writing: " << specPath << endl;
        WriteMultiFldFileIDs(infofile, filenames);
    }
    m_comm->Block();

    // Pad rank to 8char filenames, e.g. P0000000.pts
    boost::format pad("P%1$07d.pts");
    pad % m_comm->GetRank();

    // Generate full path name
    fs::path poutfile(pad.str());
    fs::path fulloutname = specPath / poutfile;

    // Return the full path to the partition for this process
    return LibUtilities::PortablePath(fulloutname);
}

// TODO: copied from FieldIO, unify
void PtsIO::WriteMultiFldFileIDs(const std::string &outFile,
                                 const std::vector<std::string> fileNames,
                                 const PtsMetaDataMap &fieldmetadatamap)
{
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    AddInfoTag(root, fieldmetadatamap);

    for (int t = 0; t < fileNames.size(); ++t)
    {
        TiXmlElement *elemIDs = new TiXmlElement("Partition");
        root->LinkEndChild(elemIDs);
        elemIDs->SetAttribute("FileName", fileNames[t]);
    }

    doc.SaveFile(outFile);
}


void PtsIO::ImportMultiFldFileIDs(const std::string &inFile,
                                    std::vector<std::string> &fileNames,
                                    FieldMetaDataMap &fieldmetadatamap)
{
    TiXmlDocument doc(inFile);
    bool loadOkay = doc.LoadFile();


    std::stringstream errstr;
    errstr << "Unable to load file: " << inFile<< std::endl;
    errstr << "Reason: " << doc.ErrorDesc() << std::endl;
    errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
    ASSERTL0(loadOkay, errstr.str());

    // Handle on XML document
    TiXmlHandle docHandle(&doc);

    // Retrieve main NEKTAR tag - XML specification states one
    // top-level element tag per file.
    TiXmlElement* master = doc.FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");

    // Partition element tag name
    std::string strPartition = "Partition";

    // First attempt to get the first Partition element
    TiXmlElement* fldfileIDs = master->FirstChildElement(strPartition.c_str());
    if (!fldfileIDs)
    {
        // If this files try previous name
        strPartition = "MultipleFldFiles";
        fldfileIDs = master->FirstChildElement("MultipleFldFiles");
    }
    ASSERTL0(fldfileIDs,
            "Unable to find 'Partition' or 'MultipleFldFiles' tag "
            "within nektar tag.");

    while (fldfileIDs)
    {
    // Read file name of partition file
    const char *attr = fldfileIDs->Attribute("FileName");
    ASSERTL0(attr, "'FileName' not provided as an attribute of '"
                    + strPartition + "' tag.");
    fileNames.push_back(std::string(attr));

//     const char* elementIDs = fldfileIDs->GetText();
//     ASSERTL0(elementIDs, "Element IDs not specified.");
//
//     std::string elementIDsStr(elementIDs);
//
//     std::vector<unsigned int> idvec;
//     ParseUtils::GenerateSeqVector(elementIDsStr.c_str(),idvec);
//
//     elementList.push_back(idvec);

    fldfileIDs = fldfileIDs->NextSiblingElement(strPartition.c_str());
    }
}


void PtsIO::ImportFieldMetaData(std::string filename,
                                 FieldMetaDataMap &fieldmetadatamap)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();

    std::stringstream errstr;
    errstr << "Unable to load file: " << filename << std::endl;
    errstr << "Reason: " << doc.ErrorDesc() << std::endl;
    errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
    ASSERTL0(loadOkay, errstr.str());

    ImportFieldMetaData(doc,fieldmetadatamap);
}


void PtsIO::ImportFieldMetaData(TiXmlDocument &doc,
                        FieldMetaDataMap &fieldmetadatamap)
{

    TiXmlHandle docHandle(&doc);
    TiXmlElement* master = 0;    // Master tag within which all data is contained.
    TiXmlElement* metadata = 0;

    master = doc.FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");
    std::string strLoop = "NEKTAR";

    // Retain original metadata structure for backwards compatibility
    // TODO: Remove old metadata format
    metadata = master->FirstChildElement("FIELDMETADATA");
    if(metadata)
    {
        TiXmlElement *param = metadata->FirstChildElement("P");

        while (param)
        {
            TiXmlAttribute *paramAttr = param->FirstAttribute();
            std::string attrName(paramAttr->Name());
            std::string paramString;

            if(attrName == "PARAM")
            {
                paramString.insert(0,paramAttr->Value());
            }
            else
            {
                ASSERTL0(false,"PARAM not provided as an attribute in FIELDMETADATA section");
            }

            // Now read body of param
            std::string paramBodyStr;

            TiXmlNode *paramBody = param->FirstChild();

            paramBodyStr += paramBody->ToText()->Value();

            fieldmetadatamap[paramString] = paramBodyStr;
            param = param->NextSiblingElement("P");
        }
    }

    // New metadata format
    metadata = master->FirstChildElement("Metadata");
    if(metadata)
    {
        TiXmlElement *param = metadata->FirstChildElement();

        while (param)
        {
            std::string paramString = param->Value();
            if (paramString != "Provenance")
            {
                // Now read body of param
                TiXmlNode *paramBody = param->FirstChild();
                std::string paramBodyStr = paramBody->ToText()->Value();

                fieldmetadatamap[paramString] = paramBodyStr;
            }
            param = param->NextSiblingElement();
        }
    }

}


// TODO: copied from FieldIO, unify
/**
 * \brief add information about provenance and fieldmetadata
 */
void PtsIO::AddInfoTag(TiXmlElement *root,
                         const PtsMetaDataMap &fieldmetadatamap)
{
    PtsMetaDataMap ProvenanceMap;

    // Nektar++ release version from VERSION file
    ProvenanceMap["NektarVersion"] = string(NEKTAR_VERSION);

    // Date/time stamp
    ptime::time_facet *facet = new ptime::time_facet("%d-%b-%Y %H:%M:%S");
    std::stringstream wss;
    wss.imbue(locale(wss.getloc(), facet));
    wss << ptime::second_clock::local_time();
    ProvenanceMap["Timestamp"] = wss.str();

    // Hostname
    boost::system::error_code ec;
    ProvenanceMap["Hostname"] = ip::host_name(ec);

    // Git information
    // If built from a distributed package, do not include this
    if (NekConstants::kGitSha1 != "GITDIR-NOTFOUND")
    {
        ProvenanceMap["GitSHA1"]   = NekConstants::kGitSha1;
        ProvenanceMap["GitBranch"] = NekConstants::kGitBranch;
    }

    TiXmlElement *infoTag = new TiXmlElement("Metadata");
    root->LinkEndChild(infoTag);

    TiXmlElement *v;
    PtsMetaDataMap::const_iterator infoit;

    TiXmlElement *provTag = new TiXmlElement("Provenance");
    infoTag->LinkEndChild(provTag);
    for (infoit = ProvenanceMap.begin(); infoit != ProvenanceMap.end(); ++infoit)
    {
        v = new TiXmlElement((infoit->first).c_str());
        v->LinkEndChild(new TiXmlText((infoit->second).c_str()));
        provTag->LinkEndChild(v);
    }

    //---------------------------------------------
    // write field info section
    if (fieldmetadatamap != NullPtsMetaDataMap)
    {
        for (infoit = fieldmetadatamap.begin(); infoit != fieldmetadatamap.end(); ++infoit)
        {
            v = new TiXmlElement((infoit->first).c_str());
            v->LinkEndChild(new TiXmlText((infoit->second).c_str()));
            infoTag->LinkEndChild(v);
        }
    }
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


}
}
