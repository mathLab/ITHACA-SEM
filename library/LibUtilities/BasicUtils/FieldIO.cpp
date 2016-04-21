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

#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/make_shared.hpp>

#include <LibUtilities/BasicConst/GitRevision.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

#include <loki/Singleton.h>

#include <fstream>
#include <set>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace ptime = boost::posix_time;
namespace ip    = boost::asio::ip;

namespace Nektar
{
namespace LibUtilities
{
FieldIOFactory &GetFieldIOFactory()
{
    typedef Loki::
        SingletonHolder<FieldIOFactory, Loki::CreateUsingNew, Loki::NoDestroy>
            Type;
    return Type::Instance();
}

/**
 * @brief Determine file type of given input file.
 *
 * This method attempts to identify the file type of a given input file @p
 * filename. It returns a string corresponding to GetFieldIOFactory() or throws
 * an assertion if it cannot be identified.
 *
 * @param filename  Input filename
 * @param comm      Communicator for parallel runs
 *
 * @return FieldIO format of @p filename.
 */
const std::string FieldIO::GetFileType(const std::string &filename,
                                       CommSharedPtr comm)
{
    // We'll use 0 => XML and 1 => HDF5.
    int code = 0;
    int size = comm->GetSize();
    int rank = comm->GetRank();

    if (size == 1 || rank == 0)
    {
        std::string datafilename;

        // If input is a directory, check for root processor file.
        if (fs::is_directory(filename))
        {
            fs::path p0file("P0000000.fld");
            fs::path fullpath = filename / p0file;
            datafilename      = PortablePath(fullpath);
        }
        else
        {
            datafilename = filename;
        }

        // Read first 8 bytes. If they correspond with magic bytes below it's an
        // HDF5 file. XML is potentially a nightmare with all the different
        // encodings so we'll just assume it's OK if it's not HDF.
        const char magic[8] = {0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a};

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
    {
        comm->Bcast(code, 0);
    }

    std::string iofmt;
    if (code == 0)
    {
        iofmt = "Xml";
    }
    else if (code == 1)
    {
        iofmt = "Hdf5";
    }
    else
    {
        // Error
        ASSERTL0(false, "Unknown file format");
    }

    return iofmt;
}

/**
 * @brief Construct a FieldIO object for a given input filename.
 *
 * This is a convenience function that takes an input filename and constructs
 * the appropriate FieldIO subclass, using FieldIO::GetFileType.
 *
 * @param session   Session reader
 * @param filename  Input filename
 *
 * @return FieldIO class reader for @p filename.
 */
FieldIOSharedPtr MakeFieldIOForFile(
    const LibUtilities::SessionReaderSharedPtr session,
    const std::string &filename)
{
    const std::string iofmt =
        FieldIO::GetFileType(filename, session->GetComm());
    return GetFieldIOFactory().CreateInstance(
        iofmt,
        session->GetComm(),
        session->DefinesCmdLineArgument("shared-filesystem"));
}

/**
 * @brief This function allows for data to be written to an FLD file when a
 * session and/or communicator is not instantiated. Typically used in utilities
 * which do not take XML input and operate in serial only.
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
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERTL0(size == 1,
                 "This static function is not available in parallel. Please"
                 "instantiate a FieldIO object for parallel use.");
    }
#endif
    CommSharedPtr c    = GetCommFactory().CreateInstance("Serial", 0, 0);
    FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance("Xml", c, false);
    f->Write(outFile, fielddefs, fielddata, fieldinfomap);
}

/**
 * @brief This function allows for data to be imported from an FLD file when a
 * session and/or communicator is not instantiated. Typically used in utilities
 * which only operate in serial.
 */
LIB_UTILITIES_EXPORT void Import(
    const std::string &infilename,
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
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERTL0(size == 1,
                 "This static function is not available in parallel. Please"
                 "instantiate a FieldIO object for parallel use.");
    }
#endif
    CommSharedPtr c    = GetCommFactory().CreateInstance("Serial", 0, 0);
    FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance("Xml", c, false);
    f->Import(infilename, fielddefs, fielddata, fieldinfomap, ElementiDs);
}

/**
 * @brief Constructor for FieldIO base class.
 */
FieldIO::FieldIO(LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    : m_comm(pComm), m_sharedFilesystem(sharedFilesystem)
{
}

void FieldIO::AddInfoTag(TagWriterSharedPtr root,
                         const FieldMetaDataMap &fieldmetadatamap)
{
    FieldMetaDataMap ProvenanceMap;

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
             infoit != fieldmetadatamap.end();
             ++infoit)
        {
            infoTag->SetAttr(infoit->first, infoit->second);
        }
    }
}

/**
 *
 */
int FieldIO::CheckFieldDefinition(const FieldDefinitionsSharedPtr &fielddefs)
{
    int i;

    if (fielddefs->m_elementIDs.size() == 0) // empty partition
    {
        return 0;
    }
    // ASSERTL0(fielddefs->m_elementIDs.size() > 0, "Fielddefs vector must
    // contain at least one element of data .");

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
            ASSERTL0(false, "Unsupported shape type.");
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
                    datasize += StdTriData::getNumberOfCoefficients(l, m) *
                                fielddefs->m_homogeneousZIDs.size();
                }
                else
                {
                    datasize += StdTriData::getNumberOfCoefficients(l, m);
                }
            }
            break;
            case eQuadrilateral:
            {
                int l = fielddefs->m_numModes[cnt++];
                int m = fielddefs->m_numModes[cnt++];
                if (fielddefs->m_numHomogeneousDir == 1)
                {
                    datasize += l * m * fielddefs->m_homogeneousZIDs.size();
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
                datasize += StdTetData::getNumberOfCoefficients(l, m, n);
            }
            break;
            case ePyramid:
            {
                int l = fielddefs->m_numModes[cnt++];
                int m = fielddefs->m_numModes[cnt++];
                int n = fielddefs->m_numModes[cnt++];
                datasize += StdPyrData::getNumberOfCoefficients(l, m, n);
            }
            break;
            case ePrism:
            {
                int l = fielddefs->m_numModes[cnt++];
                int m = fielddefs->m_numModes[cnt++];
                int n = fielddefs->m_numModes[cnt++];
                datasize += StdPrismData::getNumberOfCoefficients(l, m, n);
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
                ASSERTL0(false, "Unsupported shape type.");
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
                    datasize += StdTriData::getNumberOfCoefficients(l, m);
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
                    datasize += StdTetData::getNumberOfCoefficients(l, m, n);
                }
                break;
                case ePyramid:
                {
                    int l = fielddefs->m_numModes[cnt++];
                    int m = fielddefs->m_numModes[cnt++];
                    int n = fielddefs->m_numModes[cnt++];
                    datasize += StdPyrData::getNumberOfCoefficients(l, m, n);
                }
                break;
                case ePrism:
                {
                    int l = fielddefs->m_numModes[cnt++];
                    int m = fielddefs->m_numModes[cnt++];
                    int n = fielddefs->m_numModes[cnt++];
                    datasize += StdPrismData::getNumberOfCoefficients(l, m, n);
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
                    ASSERTL0(false, "Unsupported shape type.");
                    break;
            }
        }
    }

    return datasize;
}
}
}
