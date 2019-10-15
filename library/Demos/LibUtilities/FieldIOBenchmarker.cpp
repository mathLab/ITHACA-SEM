////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOBenchmarker.cpp
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
//  Description: Measure the performance of FieldIO XML and HDF5 classes.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/Communication/CommMpi.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <unordered_set>
#include <string>

typedef std::unordered_set<int> IntSet;

using namespace Nektar;
using namespace LibUtilities;

namespace po    = boost::program_options;
namespace berrc = boost::system::errc;

typedef std::vector<FieldDefinitionsSharedPtr> DefVec;
typedef std::vector<std::vector<NekDouble> > DatVec;

/// Struct that contains experimental setup.
struct Experiment
{
    /// Test read (false) or write (true)
    bool write;
    /// Test the HDF5 (true) or XML (false) reader/writer.
    bool hdf;
    /// If true, print additional debugging information.
    bool verbose;
    /// Number of writes to perform.
    int n;
    /// Input file source name
    std::string dataSource;
    /// Output filename
    std::string dataDest;
    /// Communicator for the experiment.
    CommSharedPtr comm;
};

typedef std::vector<NekDouble> Results;

Results TestRead(Experiment &exp);
Results TestWrite(Experiment &exp);
void PrintResults(Experiment &exp, Results &results);

int main(int argc, char *argv[])
{
    Experiment exp;
    exp.write   = false;
    exp.hdf     = false;
    exp.verbose = false;
    exp.n       = 3;
    exp.comm    = GetCommFactory().CreateInstance("ParallelMPI", argc, argv);

    po::options_description desc("Available options");
    desc.add_options()("help,h", "Produce this help message.")(
        "mode,m", po::value<char>(),
        "Choose r[ead] (default), x[ml write] or h[df5 write]")(
        "number,n", po::value<unsigned>(),
        "Number of iterations to perform, default 3")("verbose,v",
                                                      "Enable verbose mode.");

    po::options_description hidden("Hidden options");
    hidden.add_options()("input-file", po::value<std::string>(),
                         "Input filename")("output-file",
                                           po::value<std::string>());

    po::options_description cmdline_options;
    cmdline_options.add(hidden).add(desc);

    po::options_description visible("Allowed options");
    visible.add(desc);

    po::positional_options_description p;
    p.add("input-file", 1).add("output-file", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv)
                      .options(cmdline_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << desc;
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1)
    {
        std::cerr
            << "Usage: FieldIOBenchmarker [options] inputfile [outputfile]"
            << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    ASSERTL0(vm.count("input-file"), "Must specify input.");

    exp.dataSource = vm["input-file"].as<std::string>();

    if (vm.count("verbose") && exp.comm->GetRank() == 0)
    {
        exp.verbose = true;
    }

    if (vm.count("number"))
    {
        exp.n = vm["number"].as<unsigned>();
    }

    if (vm.count("mode"))
    {
        char mode = vm["mode"].as<char>();
        switch (mode)
        {
            case 'r':
                exp.write = false;
                break;
            case 'x':
                exp.write = true;
                exp.hdf   = false;
                break;
            case 'h':
                exp.write = true;
                exp.hdf   = true;
                break;
            default:
                std::cout << "Unrecognised mode: " << mode << std::endl;
                std::cout << desc << std::endl;
                return 1;
                break;
        }
    }

    Results res;
    if (exp.write)
    {
        if (vm.count("output-file"))
        {
            exp.dataDest = vm["output-file"].as<std::string>();
        }
        else
        {
            exp.dataDest = exp.dataSource + ".tmp";
        }

        res = TestWrite(exp);
    }
    else
    {
        res = TestRead(exp);
    }

    PrintResults(exp, res);
    exp.comm->Finalise();
}

/**
 * @brief Read elemental IDs from XML field file format for this rank.
 *
 * Here we read the Info.xml to figure out in which file the elements are
 * stored.  The elements are then divided amongst the ranks (in a trivial
 * decomposition based on order in the Info.xml).
 *
 * @param exp  Experiment to be run.
 * @param fio  FieldIO object to perform test on.
 *
 * @return Array containing element IDs for this rank.
 */
Array<OneD, int> ReadIDsForThisRank(Experiment &exp, FieldIOSharedPtr fio)
{
    std::vector<std::string> fileNames;
    std::vector<std::vector<unsigned int> > elementList;
    FieldMetaDataMap fieldmetadatamap;

    std::string infoFile = exp.dataSource + "/Info.xml";

    std::shared_ptr<FieldIOXml> fioXml =
        std::dynamic_pointer_cast<FieldIOXml>(fio);
    fioXml->ImportMultiFldFileIDs(infoFile, fileNames, elementList,
                                  fieldmetadatamap);

    unsigned totalEls = 0;
    std::vector<unsigned> elStartFile(elementList.size(), 0);
    std::vector<unsigned> elStopFile(elementList.size(), 0);

    for (unsigned i = 0; i < elementList.size(); ++i)
    {
        elStartFile[i] = totalEls;
        totalEls += elementList[i].size();
        elStopFile[i] = totalEls;
    }
    NekDouble elemPerNode = totalEls / (NekDouble)exp.comm->GetSize();
    unsigned elStart      = elemPerNode * exp.comm->GetRank();
    unsigned elStop       = elemPerNode * (exp.comm->GetRank() + 1);
    unsigned nEls         = elStop - elStart;

    Array<OneD, int> ElementIDs(nEls);

    for (unsigned iFile = 0, iEl = elStart; iEl < elStop;)
    {
        // Find the index of the file that contains the element index we want
        while (!(elStartFile[iFile] <= iEl && iEl < elStopFile[iFile]))
            iFile++;

        unsigned startInFile = iEl - elStartFile[iFile];
        unsigned stopInFile;

        // Determine how much of the file we want
        if (elStop > elStopFile[iFile])
        {
            // Need some of the next one too
            // Copy to the end
            stopInFile = elStopFile[iFile] - elStartFile[iFile];
        }
        else
        {
            // This is the last file we need
            // Copy up to elStop
            stopInFile = elStop - elStartFile[iFile];
        }

        // Copy the chunk
        std::memcpy(&ElementIDs[iEl - elStart],
                    &elementList[iFile][startInFile],
                    (stopInFile - startInFile) * sizeof(int));

        iEl += stopInFile - startInFile;
    }
    return ElementIDs;
}

/**
 * @brief Construct data for this rank based on element IDs and input data.
 *
 * Extract from @p inFieldDefs and @p inFieldData those elements which are in @p
 * ElementIDs and return them in the parameters @p outFieldDefs and @p
 * outFieldData.
 *
 * @param inFieldDefs    Field definitions of input file
 * @param inFieldData    Input field data corresponding to the definitions
 * @param ElementIDs     Element IDs for this rank
 * @param outFieldDefs   Output field definitions containing this rank's data
 * @param outFieldData   Output field data corresponding to the definitions
 */
void FilterDataForThisRank(const DefVec &inFieldDefs,
                           const DatVec &inFieldData,
                           Array<OneD, int> ElementIDs,
                           DefVec &outFieldDefs,
                           DatVec &outFieldData)
{
    // Create a set with all the IDs
    IntSet IDs(ElementIDs.begin(), ElementIDs.end());

    // Clear the output vectors
    outFieldDefs.clear();
    outFieldData.clear();

    // Loop through all the loaded elements and copy over if in the requested
    // set
    DefVec::const_iterator inDefIt = inFieldDefs.begin();
    DatVec::const_iterator inDatIt = inFieldData.begin();
    for (; inDefIt != inFieldDefs.end(); ++inDefIt, ++inDatIt)
    {
        FieldDefinitionsSharedPtr inDef = *inDefIt;
        // Use list to avoid endless reallocation.
        std::list<unsigned int> elOut;
        std::list<NekDouble> datOut;

        unsigned dat_per_el = inDatIt->size() / inDef->m_elementIDs.size();

        std::vector<unsigned int>::const_iterator elIt =
            inDef->m_elementIDs.begin();
        std::vector<NekDouble>::const_iterator datIt = inDatIt->begin();
        for (; elIt != inDef->m_elementIDs.end(); ++elIt, datIt += dat_per_el)
        {
            if (IDs.find(*elIt) != IDs.end())
            {
                // Copy across element id
                elOut.push_back(*elIt);
                // and data
                datOut.insert(datOut.end(), datIt, datIt + dat_per_el);
            }
        }
        if (elOut.size())
        {
            // create the outFieldDefs
            FieldDefinitionsSharedPtr defOut =
                FieldDefinitionsSharedPtr(new FieldDefinitions(
                    inDef->m_shapeType,
                    std::vector<unsigned int>(elOut.begin(), elOut.end()),
                    inDef->m_basis, inDef->m_uniOrder, inDef->m_numModes,
                    inDef->m_fields, inDef->m_numHomogeneousDir,
                    inDef->m_homogeneousLengths, inDef->m_homoStrips,
                    inDef->m_homogeneousSIDs, inDef->m_homogeneousZIDs,
                    inDef->m_homogeneousYIDs, inDef->m_points,
                    inDef->m_pointsDef, inDef->m_numPoints,
                    inDef->m_numPointsDef));
            // Add to return
            outFieldDefs.push_back(defOut);
            // create the out data vector from our list
            outFieldData.push_back(
                std::vector<NekDouble>(datOut.begin(), datOut.end()));
        }
    }
}

/**
 * @brief Read all data in those files that this rank wants (and any other data
 * in them too).
 *
 * @param Exp  Experiment to be run
 * @param outFieldDefs   Output field definitions
 * @param outFieldData   Output data corresponding to the definitions
 */
void ReadWholeFilesForThisRank(Experiment &exp,
                               DefVec &outFieldDefs,
                               DatVec &outFieldData)
{
    std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    FieldIOSharedPtr fio =
        GetFieldIOFactory().CreateInstance(ft, exp.comm, true);

    if (fs::is_directory(exp.dataSource))
    {
        Array<OneD, int> ElementIDs = ReadIDsForThisRank(exp, fio);
        FieldMetaDataMap fieldmetadatamap;

        // Load all the data from files that contain any of the IDs we want.
        fio->Import(exp.dataSource, outFieldDefs, outFieldData,
                    fieldmetadatamap, ElementIDs);
    }
    else
    {
        fio->Import(exp.dataSource, outFieldDefs, outFieldData);
    }
}

/**
 * @brief Read only the data that this rank wants.
 *
 * @param exp           Experiment setup details
 * @param outFieldDefs  Resulting output field definitions for this rank.
 * @param outFieldData  Output field data for this rank.
 */
void ReadDecomposed(Experiment &exp, DefVec &outFieldDefs, DatVec &outFieldData)
{
    std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    FieldIOSharedPtr fio =
        GetFieldIOFactory().CreateInstance(ft, exp.comm, true);

    Array<OneD, int> ElementIDs = ReadIDsForThisRank(exp, fio);
    DefVec fileFieldDefs;
    DatVec fileFieldData;
    FieldMetaDataMap fieldmetadatamap;

    // Load all the data from files that contain any of the IDs we want.
    fio->Import(exp.dataSource, fileFieldDefs, fileFieldData, fieldmetadatamap,
                ElementIDs);
    // Filter it
    FilterDataForThisRank(fileFieldDefs, fileFieldData, ElementIDs,
                          outFieldDefs, outFieldData);
}

/**
 * @brief Test read speed for this experimental setup.
 *
 * This routine performs Experiment::n reads, timing the read times and
 * returning a Results struct containing the results.
 *
 * @param exp  Experimental setup.
 *
 * @return Resulting timings.
 */
Results TestRead(Experiment &exp)
{
    if (exp.verbose)
    {
        std::cout << "Beginning read experiment with " << exp.n << " loops."
                  << std::endl;
        std::cout << "Determining file type... ";
    }

    const std::string ft = FieldIO::GetFileType(exp.dataSource, exp.comm);
    if (exp.verbose)
    {
        std::cout << ft << std::endl;
    }

    Results res(exp.n, 0.0);
    for (unsigned i = 0; i < exp.n; ++i)
    {
        if (exp.verbose)
        {
            std::cout << "Test " << i << " of " << exp.n;
        }

        std::vector<FieldDefinitionsSharedPtr> fielddefs;
        std::vector<std::vector<NekDouble> > fielddata;
        // Synchronise
        exp.comm->Block();

        NekDouble t0 = MPI_Wtime();

        ReadWholeFilesForThisRank(exp, fielddefs, fielddata);

        NekDouble t1 = MPI_Wtime();
        t1 -= t0;

        if (exp.verbose)
        {
            std::cout << ": t = " << t1 << " s" << std::endl;
        }

        res[i] = t1;
    }
    return res;
}

/**
 * @brief Test write speed for this experimental setup.
 *
 * This routine performs Experiment::n writes, timing the write times and
 * returning a Results struct containing the results.
 *
 * @param exp  Experimental setup.
 *
 * @return Resulting timings.
 */
Results TestWrite(Experiment &exp)
{
    if (exp.verbose)
        std::cout << "Reading in input: " << exp.dataSource << std::endl;

    std::vector<FieldDefinitionsSharedPtr> fielddefs;
    std::vector<std::vector<NekDouble> > fielddata;
    ReadDecomposed(exp, fielddefs, fielddata);

    std::string outtype;
    if (exp.hdf)
    {
        outtype = "Hdf5";
    }
    else
    {
        outtype = "Xml";
    }

    if (exp.verbose)
    {
        std::cout << "Beginning write (" << outtype << ") experiment with "
                  << exp.n << " loops." << std::endl;
        std::cout << "Writing to temp file: " << exp.dataDest << std::endl;
    }

    Results res(exp.n, 0);
    for (unsigned i = 0; i < exp.n; ++i)
    {
        if (exp.verbose)
        {
            std::cout << "Test " << i << " of " << exp.n << std::endl;
        }

        // Synchronise - have to do this before removing any old data in case
        // any ranks haven't closed their file yet.
        exp.comm->Block();
        // Remove any existing files
        fs::path specPath(exp.dataDest);
        try
        {
            fs::remove_all(specPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                     "Filesystem error: " + std::string(e.what()));
        }

        // Synchronise to make sure we're all at the same point.
        exp.comm->Block();

        NekDouble t0 = MPI_Wtime();

        FieldIOSharedPtr fio =
            GetFieldIOFactory().CreateInstance(outtype, exp.comm, true);

        fio->Write(exp.dataDest, fielddefs, fielddata);

        NekDouble t1 = MPI_Wtime();
        t1 -= t0;

        if (exp.verbose)
        {
            std::cout << ": t = " << t1 << " s" << std::endl;
        }

        res[i] = t1;
    }
    return res;
}

/**
 * @brief Print out the results of the timing.
 *
 * Given an experimental setup @p exp and timing in @p results, this routine
 * prints the mean time per read/write.
 *
 * @param exp      Experimental setup.
 * @param results  Results of timing for @p exp
 */
void PrintResults(Experiment &exp, Results &results)
{
    NekDouble sum   = 0.0;
    NekDouble sumSq = 0.0;

    for (Results::const_iterator it = results.begin(); it != results.end();
         ++it)
    {
        NekDouble x = *it;
        sum += x;
        sumSq += x * x;
    }

    NekDouble mean = sum / exp.n;
    // NekDouble var = sumSq / exp.n - mean*mean;
    // NekDouble std = std::sqrt(var);

    if (exp.comm->GetSize() > 1)
    {
        // Use all version since reduce to specified rank isn't wrapped.
        exp.comm->AllReduce(mean, ReduceSum);
        mean /= exp.comm->GetSize();
    }

    if (exp.comm->GetRank() == 0)
    {
        std::cout << "Mean: " << mean << std::endl;
        // std::cout << "Std: " << std << std::endl;
    }
}
