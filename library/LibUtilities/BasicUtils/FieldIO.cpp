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
#include <boost/format.hpp>

#include <LibUtilities/BasicConst/GitRevision.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

#include <chrono>
#include <ctime>
#include <ios>
#include <iomanip>
#include <fstream>
#include <set>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace berrc = boost::system::errc;
namespace ip    = boost::asio::ip;

namespace Nektar
{
namespace LibUtilities
{

std::string fldCmdFormat = SessionReader::RegisterCmdLineArgument(
    "io-format", "i", "Default input/output format (e.g. Xml, Hdf5)");

/**
 * @brief Returns the FieldIO factory.
 */
FieldIOFactory &GetFieldIOFactory()
{
    static FieldIOFactory instance;
    return instance;
}

/// Enumerator for auto-detection of FieldIO types.
enum FieldIOType {
    eXML,
    eHDF5
};


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
    FieldIOType ioType = eXML;
    int size  = comm->GetSize();
    bool root = comm->TreatAsRankZero();

    if (size == 1 || root)
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
        const unsigned char magic[8] = {
            0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a};

        std::ifstream datafile(datafilename.c_str(), std::ios_base::binary);
        ASSERTL0(datafile.good(), "Unable to open file: " + filename);

        ioType = eHDF5;
        for (unsigned i = 0; i < 8 && datafile.good(); ++i)
        {
            int byte = datafile.get();
            if (byte != magic[i])
            {
                ioType = eXML;
                break;
            }
        }
    }

    if (size > 1)
    {
        int code = (int)ioType;
        comm->Bcast(code, 0);
        ioType = (FieldIOType)code;
    }

    std::string iofmt;
    if (ioType == eXML)
    {
        iofmt = "Xml";
    }
    else if (ioType == eHDF5)
    {
        iofmt = "Hdf5";
    }
    else
    {
        // Error
        NEKERROR(ErrorUtil::efatal, "Unknown file format");
    }

    return iofmt;
}

/**
 * @brief Returns an object for the default FieldIO method.
 *
 * This function returns a FieldIO class as determined by the hard-coded default
 * (XML), which can be overridden by changing the session reader SOLVERINFO
 * variable FieldIOFormat.
 *
 * @param session  Session reader
 *
 * @return FieldIO object
 */
FieldIOSharedPtr FieldIO::CreateDefault(
    const LibUtilities::SessionReaderSharedPtr session)
{
    std::string iofmt("Xml");
    if (session->DefinesSolverInfo("IOFormat"))
    {
        iofmt = session->GetSolverInfo("IOFormat");
    }

    if (session->DefinesCmdLineArgument("io-format"))
    {
        iofmt = session->GetCmdLineArgument<std::string>("io-format");
    }

    return GetFieldIOFactory().CreateInstance(
        iofmt,
        session->GetComm(),
        session->GetSharedFilesystem());
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
FieldIOSharedPtr FieldIO::CreateForFile(
    const LibUtilities::SessionReaderSharedPtr session,
    const std::string &filename)
{
    const std::string iofmt =
        FieldIO::GetFileType(filename, session->GetComm());
    return GetFieldIOFactory().CreateInstance(
        iofmt,
        session->GetComm(),
        session->GetSharedFilesystem());
}

/**
 * @brief This function allows for data to be written to an FLD file when a
 * session and/or communicator is not instantiated. Typically used in utilities
 * which do not take XML input and operate in serial only.
 *
 * @param outFile       Output filename
 * @param fielddefs     Field definitions that define the output
 * @param fielddata     Binary field data that stores the output corresponding
 *                      to @p fielddefs.
 * @param fieldinfomap  Associated field metadata map.
 */
void Write(const std::string &outFile,
           std::vector<FieldDefinitionsSharedPtr> &fielddefs,
           std::vector<std::vector<NekDouble> > &fielddata,
           const FieldMetaDataMap &fieldinfomap,
           const bool backup)
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
    f->Write(outFile, fielddefs, fielddata, fieldinfomap, backup);
}

/**
 * @brief This function allows for data to be imported from an FLD file when a
 * session and/or communicator is not instantiated. Typically used in utilities
 * which only operate in serial.
 *
 * @param infilename    Input filename (or directory if parallel format)
 * @param fielddefs     On return contains field definitions as read from the
 *                      input.
 * @param fielddata     On return, contains binary field data that stores the
 *                      input corresponding to @p fielddefs.
 * @param fieldinfo     On returnm, contains the associated field metadata map.
 * @param ElementIDs    Element IDs that lie on this processor, which can be
 *                      optionally supplied to avoid reading the entire file on
 *                      each processor.
 */
LIB_UTILITIES_EXPORT void Import(
    const std::string &infilename,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    FieldMetaDataMap &fieldinfomap,
    const Array<OneD, int> &ElementIDs)
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
    const std::string iofmt = FieldIO::GetFileType(infilename, c);
    FieldIOSharedPtr f = GetFieldIOFactory().CreateInstance(iofmt, c, false);
    f->Import(infilename, fielddefs, fielddata, fieldinfomap, ElementIDs);
}

/**
 * @brief Constructor for FieldIO base class.
 */
FieldIO::FieldIO(LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    : m_comm(pComm), m_sharedFilesystem(sharedFilesystem)
{
}

/**
 * @brief Add provenance information to the field metadata map.
 *
 * This routine adds some basic provenance information to the field metadata to
 * enable better tracking of version information:
 *
 *   - Nektar++ version
 *   - Date and time of simulation
 *   - Hostname of the machine the simulation was performed on
 *   - git SHA1 and branch name, if Nektar++ was compiled from git and not
 *     e.g. a tarball.
 *
 * @param root              Root tag, which is encapsulated using the TagWriter
 *                          structure to enable multi-file format support.
 * @param fieldmetadatamap  Any existing field metadata.
 */
void FieldIO::AddInfoTag(TagWriterSharedPtr root,
                         const FieldMetaDataMap &fieldmetadatamap)
{
    FieldMetaDataMap ProvenanceMap;

    // Nektar++ release version from VERSION file
    ProvenanceMap["NektarVersion"] = std::string(NEKTAR_VERSION);

    // Date/time stamp
    auto now = std::chrono::system_clock::now();
    auto now_t = std::chrono::system_clock::to_time_t(now);
    auto now_tm = *std::localtime(&now_t);
    char buffer[128];
    strftime(buffer, sizeof(buffer), "%d-%b-%Y %H:%M:%S", &now_tm);
    ProvenanceMap["Timestamp"] = buffer;

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

    TagWriterSharedPtr provTag = infoTag->AddChild("Provenance");
    for (auto &infoit : ProvenanceMap)
    {
        provTag->SetAttr(infoit.first, infoit.second);
    }

    //---------------------------------------------
    // write field info section
    if (fieldmetadatamap != NullFieldMetaDataMap)
    {
        for (auto &infoit : fieldmetadatamap)
        {
            infoTag->SetAttr(infoit.first, infoit.second);
        }
    }
}

/**
 * @brief Set up the filesystem ready for output.
 *
 * This function sets up the output given an output filename @p outname. This
 * will therefore remove any file or directory with the desired output filename
 * and return the absolute path to the output.
 *
 * If @p perRank is set, a new directory will be created to contain one file per
 * process rank.
 *
 * @param outname   Desired output filename.
 * @param perRank   True if one file-per-rank output is required.
 *
 * @return Absolute path to resulting file.
 */
std::string FieldIO::SetUpOutput(const std::string outname, bool perRank, bool backup)
{
    ASSERTL0(!outname.empty(), "Empty path given to SetUpOutput()");

    int nprocs = m_comm->GetSize();
    bool root  = m_comm->TreatAsRankZero();

    // Path to output: will be directory if parallel, normal file if
    // serial.
    fs::path specPath(outname), fulloutname;

    // in case we are rank 0 or not on a shared filesystem, check if the specPath already exists
    if (backup && (root || !m_sharedFilesystem) && fs::exists(specPath))
    {
        // rename. foo/bar_123.chk -> foo/bar_123.bak0.chk and in case
        // foo/bar_123.bak0.chk already exists, foo/bar_123.chk -> foo/bar_123.bak1.chk
        fs::path bakPath = specPath;
        int cnt = 0;
        while (fs::exists(bakPath))
        {
            bakPath = specPath.parent_path();
            bakPath += specPath.stem();
            bakPath += fs::path(".bak" + std::to_string(cnt++));
            bakPath += specPath.extension();
        }
        std::cout << "renaming " << specPath << " -> " << bakPath << std::endl;
        try
        {
            fs::rename(specPath, bakPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                     "Filesystem error: " + std::string(e.what()));
        }
    }

    // wait until rank 0 has moved the old specPath and the changes
    // have propagated through the filesystem
    if (backup)
    {
        m_comm->Block();
        int exists = 1;
        while (exists && perRank)
        {
            exists = fs::exists(specPath);
            m_comm->AllReduce(exists, ReduceMax);
        }
    }

    if (nprocs == 1)
    {
        fulloutname = specPath;
    }
    else
    {
        // Guess at filename that might belong to this process.
        boost::format pad("P%1$07d.%2$s");
        pad % m_comm->GetRank() % GetFileEnding();

        // Generate full path name
        fs::path poutfile(pad.str());
        fulloutname = specPath / poutfile;
    }

    // Remove any existing file which is in the way
    if (m_comm->RemoveExistingFiles() && !backup)
    {
        if (m_sharedFilesystem)
        {
            // First, each process clears up its .fld file. This might or might
            // not be there (we might have changed numbers of processors between
            // runs, for example), but we can try anyway.
            try
            {
                fs::remove_all(fulloutname);
            }
            catch (fs::filesystem_error &e)
            {
                ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                         "Filesystem error: " + std::string(e.what()));
            }
        }

        m_comm->Block();

        // Now get rank 0 processor to tidy everything else up.
        if (root || !m_sharedFilesystem)
        {
            try
            {
                fs::remove_all(specPath);
            }
            catch (fs::filesystem_error &e)
            {
                ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                         "Filesystem error: " + std::string(e.what()));
            }
        }

        // wait until rank 0 has removed specPath and the changes
        // have propagated through the filesystem
        m_comm->Block();
        int exists = 1;
        while (exists && perRank)
        {
            exists = fs::exists(specPath);
            m_comm->AllReduce(exists, ReduceMax);
        }
    }

    if (root)
    {
        std::cout << "Writing: " << specPath;
    }

    // serial processing just add ending.
    if (nprocs == 1)
    {
        return LibUtilities::PortablePath(specPath);
    }

    // Create the destination directory
    if (perRank)
    {
        try
        {
            if (root || !m_sharedFilesystem)
            {
                fs::create_directory(specPath);
            }
        }
        catch (fs::filesystem_error &e)
        {
            NEKERROR(ErrorUtil::efatal,
                     "Filesystem error: " + std::string(e.what()));
        }

        m_comm->Block();

        // Sit in a loop and make sure target directory has been created
        int created = 0;
        while (!created)
        {
            created = fs::is_directory(specPath);
            m_comm->AllReduce(created, ReduceMin);
        }
    }
    else
    {
        fulloutname = specPath;
    }

    // Return the full path to the partition for this process
    return LibUtilities::PortablePath(fulloutname);
}

/**
 * @brief Check field definitions for correctness and return storage size.
 *
 * @param fielddefs  Field definitions to check.
 */
int FieldIO::CheckFieldDefinition(const FieldDefinitionsSharedPtr &fielddefs)
{
    int i;

    if (fielddefs->m_elementIDs.size() == 0) // empty partition
    {
        return 0;
    }

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
            NEKERROR(ErrorUtil::efatal, "Unsupported shape type.");
            break;
    }

    size_t datasize = 0;

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
                    datasize += l * fielddefs->m_homogeneousZIDs.size();
                    cnt++;
                }
                else if (fielddefs->m_numHomogeneousDir == 2)
                {
                    datasize += l * fielddefs->m_homogeneousYIDs.size();
                    cnt += 2;
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
                NEKERROR(ErrorUtil::efatal, "Unsupported shape type.");
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
                {
                    int l = fielddefs->m_numModes[cnt++];
                    if (fielddefs->m_numHomogeneousDir == 1)
                    {
                        datasize += l * fielddefs->m_homogeneousZIDs.size();
                        cnt++;
                    }
                    else if (fielddefs->m_numHomogeneousDir == 2)
                    {
                        datasize += l * fielddefs->m_homogeneousYIDs.size();
                        cnt += 2;
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
                        cnt++;
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
                        cnt++;
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
                    NEKERROR(ErrorUtil::efatal, "Unsupported shape type.");
                    break;
            }
        }
    }

    return (int)datasize;
}
}
}
