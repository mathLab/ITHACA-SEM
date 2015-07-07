////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIO.cpp
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
//  Description: I/O routines relating to Fields
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/make_shared.hpp>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicConst/GitRevision.h>

#include <loki/Singleton.h>

#include "zlib.h"
#include <set>
#include <fstream>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

// Buffer size for zlib compression/decompression
#define CHUNK 16384

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
        FieldIOFactory& GetFieldIOFactory()
        {
            typedef Loki::SingletonHolder<FieldIOFactory, Loki::CreateUsingNew,
                    Loki::NoDestroy> Type;
            return Type::Instance();
        }

        const std::string FieldIO::GetFileType(const std::string& filename,
                CommSharedPtr comm)
        {
            // We'll use 0 => XML
            // and 1 => HDF5
            int code = 0;

            int size = comm->GetSize();
            int rank = comm->GetRank();

            if (size == 1 || rank == 0)
            {
                std::string datafilename;
                if (fs::is_directory(filename)) // check to see that infile is a directory
                {
                    fs::path p0file("P0000000.fld");
                    fs::path fullpath = filename / p0file;
                    datafilename = PortablePath(fullpath);
                }
                else
                {
                    datafilename = filename;
                }
                // Read first 8 bytes
                // If they are (in hex) 89  48  44  46  0d  0a  1a  0a
                // then it's an HDF5 file.

                // XML is potentially a nightmare with all the different encodings
                // so we'll just assume it's OK if it's not HDF
                const char magic[8] =
                { 0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a };

                std::ifstream datafile(datafilename.c_str(), ios_base::binary);
                char filedata[8];
                datafile.read(filedata, 8);

                code = 1;
                for (unsigned i = 0; i < 8; ++i)
                {
                    if (filedata[i] != magic[i])
                    {
                        code = 0;
                        break;
                    }
                }
            }

            if (size > 1)
                comm->Bcast(code, 0);

            std::string iofmt;
            if (code == 0)
                iofmt = "Xml";
            else if (code == 1)
                iofmt = "Hdf5";
            else
            // Error
            ASSERTL0(false, "Unknown file format");
            return iofmt;
        }

        FieldIOSharedPtr MakeFieldIOForFile(
                const LibUtilities::SessionReaderSharedPtr session,
                const std::string& filename)
        {
            const std::string iofmt = FieldIO::GetFileType(filename,
                    session->GetComm());
            return GetFieldIOFactory().CreateInstance(iofmt, session->GetComm());
        }
        /**
         * This function allows for data to be written to an FLD file when a
         * session and/or communicator is not instantiated. Typically used in
         * utilities which do not take XML input and operate in serial only.
         */
        void Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldinfomap)
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
            FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance("Xml", c);
            f->Write(outFile, fielddefs, fielddata, fieldinfomap);
        }

        /**
         * This function allows for data to be imported from an FLD file when
         * a session and/or communicator is not instantiated. Typically used in
         * utilities which only operate in serial.
         */
        LIB_UTILITIES_EXPORT void Import(const std::string& infilename,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                FieldMetaDataMap &fieldinfomap,
                const Array<OneD, int> ElementiDs)
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
            FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance("Xml", c);
            f->Import(infilename, fielddefs, fielddata, fieldinfomap,
                    ElementiDs);
        }

        /**
         *
         */
        FieldIO::FieldIO(LibUtilities::CommSharedPtr pComm) :
                m_comm(pComm)
        {
        }

        /**
         *
         */
        void FieldIO::WriteMultiFldFileIDs(const std::string &outFile,
                const std::vector<std::string> fileNames,
                std::vector<std::vector<unsigned int> > &elementList,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            TiXmlDocument doc;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            ASSERTL0(fileNames.size() == elementList.size(),
                    "Outfile names and list of elements ids does not match");

            TiXmlElement * root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            AddInfoTag(root, fieldmetadatamap);

            for (int t = 0; t < fileNames.size(); ++t)
            {
                if (elementList[t].size())
                {
                    TiXmlElement * elemIDs = new TiXmlElement("Partition");
                    root->LinkEndChild(elemIDs);

                    elemIDs->SetAttribute("FileName", fileNames[t]);

                    string IDstring;

                    GenerateSeqString(elementList[t], IDstring);

                    elemIDs->LinkEndChild(new TiXmlText(IDstring));
                }
            }

            doc.SaveFile(outFile);
        }

        void FieldIO::Import(const std::string& infilename,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                FieldMetaDataMap &fieldinfomap,
                const Array<OneD, int> ElementIDs)
        {
            std::string infile = infilename;

            fs::path pinfilename(infilename);

            if (fs::is_directory(pinfilename)) // check to see that infile is a directory
            {
                fs::path infofile("Info.xml");
                fs::path fullpath = pinfilename / infofile;
                infile = PortablePath(fullpath);

                std::vector < std::string > filenames;
                std::vector < std::vector<unsigned int>
                        > elementIDs_OnPartitions;

                ImportMultiFldFileIDs(infile, filenames,
                        elementIDs_OnPartitions, fieldinfomap);

                // Load metadata
                if (GetClassName() == "Xml")
                {
                    ImportFieldMetaData(infile, fieldinfomap);
                }
                else
                {
                    FieldIOSharedPtr infoReader =
                            GetFieldIOFactory().CreateInstance("Xml", m_comm);
                    infoReader->ImportFieldMetaData(infile, fieldinfomap);
                }

                if (ElementIDs == NullInt1DArray) //load all fields
                {
                    for (int i = 0; i < filenames.size(); ++i)
                    {
                        fs::path pfilename(filenames[i]);
                        fullpath = pinfilename / pfilename;
                        string fname = PortablePath(fullpath);
                        v_ImportFile(fname, fielddefs, fielddata,
                                DataSourceSharedPtr());
                    }

                }
                else // only load relevant partitions
                {
                    int i, j;
                    map<int, vector<int> > FileIDs;
                    map<int, vector<int> >::iterator it;
                    set<int> LoadFile;

                    for (i = 0; i < elementIDs_OnPartitions.size(); ++i)
                    {
                        for (j = 0; j < elementIDs_OnPartitions[i].size(); ++j)
                        {
                            FileIDs[elementIDs_OnPartitions[i][j]].push_back(i);
                        }
                    }

                    for (i = 0; i < ElementIDs.num_elements(); ++i)
                    {
                        it = FileIDs.find(ElementIDs[i]);
                        if (it != FileIDs.end())
                        {
                            for (j = 0; j < it->second.size(); ++j)
                            {
                                LoadFile.insert(it->second[j]);
                            }
                        }
                    }

                    set<int>::iterator iter;
                    for (iter = LoadFile.begin(); iter != LoadFile.end();
                            ++iter)
                    {
                        fs::path pfilename(filenames[*iter]);
                        fullpath = pinfilename / pfilename;
                        string fname = PortablePath(fullpath);
                        v_ImportFile(fname, fielddefs, fielddata,
                                DataSourceSharedPtr());
                    }
                }
            }
            else // serial format case
            {
                DataSourceSharedPtr dfile = ImportFieldMetaData(infile,
                        fieldinfomap);
                v_ImportFile(infile, fielddefs, fielddata, dfile);
            }
        }

        /** 
         *
         */
        void FieldIO::ImportMultiFldFileIDs(const std::string &inFile,
                std::vector<std::string> &fileNames,
                std::vector<std::vector<unsigned int> > &elementList,
                FieldMetaDataMap &fieldmetadatamap)
        {
            TiXmlDocument doc(inFile);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << inFile << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column "
                    << doc.ErrorCol() << std::endl;
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
            TiXmlElement* fldfileIDs = master->FirstChildElement(
                    strPartition.c_str());
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
                ASSERTL0(attr,
                        "'FileName' not provided as an attribute of '"
                                + strPartition + "' tag.");
                fileNames.push_back(std::string(attr));

                const char* elementIDs = fldfileIDs->GetText();
                ASSERTL0(elementIDs, "Element IDs not specified.");

                std::string elementIDsStr(elementIDs);

                std::vector<unsigned int> idvec;
                ParseUtils::GenerateSeqVector(elementIDsStr.c_str(), idvec);

                elementList.push_back(idvec);

                fldfileIDs = fldfileIDs->NextSiblingElement(
                        strPartition.c_str());
            }
        }

        /**
         * \brief add information about provenance and fieldmetadata
         */
        void FieldIO::AddInfoTag(TiXmlElement * root,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            TagWriterSharedPtr w = boost::make_shared < XmlTagWriter > (root);
            AddInfoTag(w, fieldmetadatamap);
        }

        TagWriter::~TagWriter()
        {
        }

        XmlTagWriter::XmlTagWriter(TiXmlElement* elem) :
                m_El(elem)
        {
        }
        TagWriterSharedPtr XmlTagWriter::AddChild(const std::string& name)
        {
            TiXmlElement* child = new TiXmlElement(name.c_str());
            m_El->LinkEndChild(child);
            return TagWriterSharedPtr(new XmlTagWriter(child));
        }

        void XmlTagWriter::SetAttr(const std::string& key,
                const std::string& val)
        {
            TiXmlElement* child = new TiXmlElement(key.c_str());
            child->LinkEndChild(new TiXmlText(val.c_str()));
            m_El->LinkEndChild(child);
        }

        void FieldIO::AddInfoTag(TagWriterSharedPtr root,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            FieldMetaDataMap ProvenanceMap;

            // Nektar++ release version from VERSION file
            ProvenanceMap["NektarVersion"] = string(NEKTAR_VERSION);

            // Date/time stamp
            ptime::time_facet *facet = new ptime::time_facet(
                    "%d-%b-%Y %H:%M:%S");
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
                ProvenanceMap["GitSHA1"] = NekConstants::kGitSha1;
                ProvenanceMap["GitBranch"] = NekConstants::kGitBranch;
            }

            TagWriterSharedPtr infoTag = root->AddChild("Metadata");

            FieldMetaDataMap::const_iterator infoit;

            TagWriterSharedPtr provTag = infoTag->AddChild("Provenance");
            for (infoit = ProvenanceMap.begin(); infoit != ProvenanceMap.end();
                    ++infoit)
            {
                provTag->SetAttr(infoit->first, infoit->second);
            }

            //---------------------------------------------
            // write field info section
            if (fieldmetadatamap != NullFieldMetaDataMap)
            {
                for (infoit = fieldmetadatamap.begin();
                        infoit != fieldmetadatamap.end(); ++infoit)
                {
                    infoTag->SetAttr(infoit->first, infoit->second);
                }
            }
        }

        /**
         *
         */
        void FieldIO::GenerateSeqString(
                const std::vector<unsigned int> &elmtids, std::string &idString)
        {
            std::stringstream idStringStream;
            bool setdash = true;
            unsigned int endval;

            idStringStream << elmtids[0];
            for (int i = 1; i < elmtids.size(); ++i)
            {
                if (elmtids[i] == elmtids[i - 1] + 1)
                {
                    if (setdash)
                    {
                        idStringStream << "-";
                        setdash = false;
                    }

                    if (i == elmtids.size() - 1) // last element
                    {
                        idStringStream << elmtids[i];
                    }
                    else
                    {
                        endval = elmtids[i];
                    }
                }
                else
                {
                    if (setdash == false) // finish off previous dash sequence
                    {
                        idStringStream << endval;
                        setdash = true;
                    }

                    idStringStream << "," << elmtids[i];
                }
            }
            idString = idStringStream.str();
        }

        /**
         *
         */
        std::string FieldIO::SetUpOutput(const std::string outname,
                const std::vector<FieldDefinitionsSharedPtr>& fielddefs,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            ASSERTL0(!outname.empty(), "Empty path given to SetUpOutput()");

            int nprocs = m_comm->GetSize();
            int rank = m_comm->GetRank();

            // Directory name if in parallel, regular filename if in serial
            fs::path specPath(outname);

            // Remove any existing file which is in the way
            if (m_comm->RemoveExistingFiles())
            {
                if (rank == 0)
                {
                    try
                    {
                        fs::remove_all(specPath);
                    } catch (fs::filesystem_error& e)
                    {
                        ASSERTL0(
                                e.code().value()
                                        == berrc::no_such_file_or_directory,
                                "Filesystem error: " + string(e.what()));
                    }
                }
            }

            // serial processing just add ending.
            if (nprocs == 1)
            {
                cout << "Writing: " << specPath << endl;
                return LibUtilities::PortablePath(specPath);
            }

            // Compute number of elements on this process and share with other
            // processes. Also construct list of elements on this process from
            // available vector of field definitions.
            std::vector<unsigned int> elmtnums(nprocs, 0);
            std::vector<unsigned int> idlist;
            int i;
            for (i = 0; i < fielddefs.size(); ++i)
            {
                elmtnums[rank] += fielddefs[i]->m_elementIDs.size();
                idlist.insert(idlist.end(), fielddefs[i]->m_elementIDs.begin(),
                        fielddefs[i]->m_elementIDs.end());
            }
            m_comm->AllReduce(elmtnums, LibUtilities::ReduceMax);

            // Create the destination directory
            try
            {
                fs::create_directory(specPath);
            } catch (fs::filesystem_error& e)
            {
                ASSERTL0(false, "Filesystem error: " + string(e.what()));
            }

            // Collate per-process element lists on root process to generate
            // the info file.
            if (rank == 0)
            {
                std::vector < std::vector<unsigned int> > ElementIDs(nprocs);

                // Populate the list of element ID lists from all processes
                ElementIDs[0] = idlist;
                for (i = 1; i < nprocs; ++i)
                {
                    std::vector<unsigned int> tmp(elmtnums[i]);
                    m_comm->Recv(i, tmp);
                    ElementIDs[i] = tmp;
                }

                // Set up output names
                std::vector < std::string > filenames;
                for (int i = 0; i < nprocs; ++i)
                {
                    boost::format pad("P%1$07d.fld");
                    pad % i;
                    filenames.push_back(pad.str());
                }

                // Write the Info.xml file
                string infofile = LibUtilities::PortablePath(
                        specPath / fs::path("Info.xml"));

                cout << "Writing: " << specPath << endl;
                WriteMultiFldFileIDs(infofile, filenames, ElementIDs,
                        fieldmetadatamap);
            }
            else
            {
                // Send this process's ID list to the root process
                m_comm->Send(0, idlist);
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

        /**
         *
         */
        int FieldIO::CheckFieldDefinition(
                const FieldDefinitionsSharedPtr &fielddefs)
        {
            int i;

            if (fielddefs->m_elementIDs.size() == 0) // empty partition
            {
                return 0;
            }
            //ASSERTL0(fielddefs->m_elementIDs.size() > 0, "Fielddefs vector must contain at least one element of data .");

            unsigned int numbasis = 0;

            // Determine nummodes vector lists are correct length
            switch (fielddefs->m_shapeType)
            {
                case eSegment:
                    numbasis = 1;
                    if (fielddefs->m_numHomogeneousDir)
                    {
                        numbasis += fielddefs->m_numHomogeneousDir;
                    }

                    break;
                case eTriangle:
                case eQuadrilateral:
                    if (fielddefs->m_numHomogeneousDir)
                    {
                        numbasis = 3;
                    }
                    else
                    {
                        numbasis = 2;
                    }
                    break;
                case eTetrahedron:
                case ePyramid:
                case ePrism:
                case eHexahedron:
                    numbasis = 3;
                    break;
                default:
                    ASSERTL0(false, "Unsupported shape type.")
                    ;
                    break;
            }

            unsigned int datasize = 0;

            ASSERTL0(fielddefs->m_basis.size() == numbasis,
                    "Length of basis vector is incorrect");

            if (fielddefs->m_uniOrder == true)
            {
                unsigned int cnt = 0;
                // calculate datasize
                switch (fielddefs->m_shapeType)
                {
                    case eSegment:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        if (fielddefs->m_numHomogeneousDir == 1)
                        {
                            datasize += l * fielddefs->m_numModes[cnt++];
                        }
                        else if (fielddefs->m_numHomogeneousDir == 2)
                        {
                            int m = fielddefs->m_numModes[cnt++];
                            datasize += l * m * fielddefs->m_numModes[cnt++];
                        }
                        else
                        {
                            datasize += l;
                        }
                    }
                        break;
                    case eTriangle:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];

                        if (fielddefs->m_numHomogeneousDir == 1)
                        {
                            datasize += StdTriData::getNumberOfCoefficients(l,
                                    m) * fielddefs->m_homogeneousZIDs.size();
                        }
                        else
                        {
                            datasize += StdTriData::getNumberOfCoefficients(l,
                                    m);
                        }
                    }
                        break;
                    case eQuadrilateral:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        if (fielddefs->m_numHomogeneousDir == 1)
                        {
                            datasize += l * m
                                    * fielddefs->m_homogeneousZIDs.size();
                        }
                        else
                        {
                            datasize += l * m;
                        }
                    }
                        break;
                    case eTetrahedron:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdTetData::getNumberOfCoefficients(l, m,
                                n);
                    }
                        break;
                    case ePyramid:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdPyrData::getNumberOfCoefficients(l, m,
                                n);
                    }
                        break;
                    case ePrism:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdPrismData::getNumberOfCoefficients(l, m,
                                n);
                    }
                        break;
                    case eHexahedron:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += l * m * n;
                    }
                        break;
                    default:
                        ASSERTL0(false, "Unsupported shape type.")
                        ;
                        break;
                }

                datasize *= fielddefs->m_elementIDs.size();
            }
            else
            {
                unsigned int cnt = 0;
                // calculate data length
                for (i = 0; i < fielddefs->m_elementIDs.size(); ++i)
                {
                    switch (fielddefs->m_shapeType)
                    {
                        case eSegment:
                            datasize += fielddefs->m_numModes[cnt++];
                            break;
                        case eTriangle:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            datasize += StdTriData::getNumberOfCoefficients(l,
                                    m);
                        }
                            break;
                        case eQuadrilateral:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            datasize += l * m;
                        }
                            break;
                        case eTetrahedron:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdTetData::getNumberOfCoefficients(l,
                                    m, n);
                        }
                            break;
                        case ePyramid:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdPyrData::getNumberOfCoefficients(l,
                                    m, n);
                        }
                            break;
                        case ePrism:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdPrismData::getNumberOfCoefficients(l,
                                    m, n);
                        }
                            break;
                        case eHexahedron:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += l * m * n;
                        }
                            break;
                        default:
                            ASSERTL0(false, "Unsupported shape type.")
                            ;
                            break;
                    }
                }
            }

            return datasize;
        }
    }
}
