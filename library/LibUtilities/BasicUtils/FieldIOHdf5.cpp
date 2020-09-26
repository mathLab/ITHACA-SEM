////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOHdf5.cpp
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
//  Description: I/O routines relating to Fields into HDF
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIOHdf5.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <unordered_set>
#include <functional>

namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{
namespace H5
{

template <> inline DataTypeSharedPtr DataTypeTraits<BasisType>::GetType()
{
    return PredefinedDataType::Native<int>();
}

}

std::string FieldIOHdf5::className =
    GetFieldIOFactory().RegisterCreatorFunction(
        "Hdf5", FieldIOHdf5::create, "HDF5-based output of field data.");

/// Version of the Nektar++ HDF5 format, which is embedded into the main NEKTAR
/// group as an attribute.
const unsigned int FieldIOHdf5::FORMAT_VERSION = 1;

// The following definitions allow us to consistently refer to indexes pulled
// out of the various datasets.

/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of elements in decomposition (i.e. field definition).
const unsigned int FieldIOHdf5::ELEM_DCMP_IDX  = 0;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of data points in decomposition (i.e. field
/// definition).
const unsigned int FieldIOHdf5::VAL_DCMP_IDX   = 1;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of elements multiplied by the dimension of the
/// element, giving number of modes when variable polynomial order is defined.
const unsigned int FieldIOHdf5::ORDER_DCMP_IDX = 2;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of y-planes for homogeneous
/// simulations.
const unsigned int FieldIOHdf5::HOMY_DCMP_IDX  = 3;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of z-planes for homogeneous
/// simulations.
const unsigned int FieldIOHdf5::HOMZ_DCMP_IDX  = 4;
/// A helper for FieldIOHdf5::v_Write and FieldIOHdf5::v_Import. Describes the
/// position of the number of the number of strips for homogeneous simulations.
const unsigned int FieldIOHdf5::HOMS_DCMP_IDX  = 5;
/// The hash of the field definition information, which defines the name of the
/// attribute containing the field definition itself.
const unsigned int FieldIOHdf5::HASH_DCMP_IDX  = 6;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the decomposition per field definition.
const unsigned int FieldIOHdf5::MAX_DCMPS      = FieldIOHdf5::HASH_DCMP_IDX + 1;

/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// elements in the cnt array.
const unsigned int FieldIOHdf5::ELEM_CNT_IDX   = 0;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// data points in the cnt array.
const unsigned int FieldIOHdf5::VAL_CNT_IDX    = 1;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// order points in the cnt array.
const unsigned int FieldIOHdf5::ORDER_CNT_IDX  = 2;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous y-planes in the cnt array.
const unsigned int FieldIOHdf5::HOMY_CNT_IDX   = 3;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous z-planes in the cnt array.
const unsigned int FieldIOHdf5::HOMZ_CNT_IDX   = 4;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous strips in the cnt array.
const unsigned int FieldIOHdf5::HOMS_CNT_IDX   = 5;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the cnt array per field definition.
const unsigned int FieldIOHdf5::MAX_CNTS       = FieldIOHdf5::HOMS_CNT_IDX + 1;

/// A helper for FieldIOHdf5::v_Write. Describes the position of the element IDs
/// within the indexing set.
const unsigned int FieldIOHdf5::IDS_IDX_IDX    = 0;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the data size
/// within the indexing set.
const unsigned int FieldIOHdf5::DATA_IDX_IDX   = 1;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the element
/// order within the indexing set.
const unsigned int FieldIOHdf5::ORDER_IDX_IDX  = 2;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// y-planes within the indexing set.
const unsigned int FieldIOHdf5::HOMY_IDX_IDX   = 3;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// z-planes within the indexing set.
const unsigned int FieldIOHdf5::HOMZ_IDX_IDX   = 4;
/// A helper for FieldIOHdf5::v_Write. Describes the position of the number of
/// homogeneous strips within the indexing set.
const unsigned int FieldIOHdf5::HOMS_IDX_IDX   = 5;
/// A helper for FieldIOHdf5::v_Write. Describes the maximum number of items in
/// the indexing set.
const unsigned int FieldIOHdf5::MAX_IDXS       = FieldIOHdf5::HOMS_IDX_IDX + 1;

/**
 * @brief Construct the FieldIO object for HDF5 output.
 *
 * @param pComm              Communicator object.
 * @param sharedFilesystem   True if this system has a shared filesystem.
 */
FieldIOHdf5::FieldIOHdf5(LibUtilities::CommSharedPtr pComm,
                         bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}

/**
 * @brief Write a HDF5 file to @p outFile given the field definitions @p
 * fielddefs, field data @p fielddata and metadata @p fieldmetadatamap.
 *
 * The writing strategy for HDF5 output is as follows:
 *
 *   - Each rank determines the amount of data needed to be written into each
 *     dataset.
 *   - Each rank communicates its decomposition information to the root process.
 *   - The root processor initialises the output structure, writes the
 *     decomposition dataset and all the field definition information.
 *   - Other ranks may have field definitions that do not belong to the root
 *     process, in which case they open the file and append this (since
 *     attributes cannot be written in parallel).
 *   - Each of the other ranks writes their data contributions to the rest of
 *     the set.
 *
 * @param outFile           Output filename.
 * @param fielddefs         Input field definitions.
 * @param fielddata         Input field data.
 * @param fieldmetadatamap  Field metadata.
 */
void FieldIOHdf5::v_Write(const std::string &outFile,
                          std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                          std::vector<std::vector<NekDouble> > &fielddata,
                          const FieldMetaDataMap &fieldmetadatamap,
                          const bool backup)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::v_Write(): ";
    double tm0 = 0.0, tm1 = 0.0;

    if (m_comm->GetRank() == 0)
    {
        tm0 = m_comm->Wtime();
    }

    SetUpOutput(outFile, false, backup);

    // We make a number of assumptions in this code:
    //   1. All element ids have the same type: unsigned int
    //   2. All elements within a given field have the same number of values
    //   3. All element values have the same type, NekDouble

    // Determine the root MPI process, i.e., the lowest ranked process handling
    // nMaxFields fields, that will be responsible for writing our file.
    ASSERTL1(fielddefs.size() == fielddata.size(),
             prfx.str() + "fielddefs and fielddata have incompatible lengths.");

    size_t nFields    = fielddefs.size();
    size_t nMaxFields = nFields;
    m_comm->AllReduce(nMaxFields, LibUtilities::ReduceMax);

    int root_rank = -1;
    bool amRoot = false;
    LibUtilities::CommSharedPtr max_fields_comm;

    if (m_comm->GetSize() > 1)
    {
        max_fields_comm = m_comm->CommCreateIf((nFields == nMaxFields) ? 1 : 0);
    }
    else
    {
        max_fields_comm = m_comm;
    }

    if (max_fields_comm)
    {
        int rank  = m_comm->GetRank();
        root_rank = rank;
        max_fields_comm->AllReduce(root_rank, LibUtilities::ReduceMin);
        amRoot = (rank == root_rank);
        if (!amRoot)
        {
            root_rank = -1;
        }
    }

    m_comm->AllReduce(root_rank, LibUtilities::ReduceMax);
    ASSERTL1(root_rank >= 0 && root_rank < m_comm->GetSize(),
             prfx.str() + "invalid root rank.");

    std::vector<uint64_t> decomps(nMaxFields * MAX_DCMPS, 0);
    std::vector<uint64_t> all_hashes(nMaxFields * m_comm->GetSize(), 0);
    std::vector<uint64_t> cnts(MAX_CNTS, 0);
    std::vector<std::string> fieldNames(nFields);
    std::vector<std::string> shapeStrings(nFields);
    std::vector<std::vector<NekDouble> > homoLengths(nFields);
    std::vector<std::vector<unsigned int> > homoSIDs(nFields),
        homoYIDs(nFields), homoZIDs(nFields);
    std::vector<std::vector<unsigned int> > numModesPerDirVar(nFields);
    std::vector<std::string> numModesPerDirUni(nFields);

    int homDim = -1;
    int varOrder = 0;

    for (int f = 0; f < nFields; ++f)
    {
        if (!fielddefs[f]->m_uniOrder)
        {
            varOrder = 1;
            break;
        }
    }

    m_comm->AllReduce(varOrder, LibUtilities::ReduceMax);

    // Calculate the total number of elements handled by this MPI process and
    // the total number of bytes required to store the elements. Base the name
    // of each field on the hash of the field definition.
    for (int f = 0; f < nFields; ++f)
    {
        ASSERTL1(fielddata[f].size() > 0,
                 prfx.str() +
                     "fielddata vector must contain at least one value.");
        ASSERTL1(fielddata[f].size() ==
                     fielddefs[f]->m_fields.size() *
                         CheckFieldDefinition(fielddefs[f]),
                 prfx.str() + "fielddata vector has invalid size.");

        std::size_t nFieldElems = fielddefs[f]->m_elementIDs.size();
        std::size_t nElemVals   = fielddata[f].size();

        decomps[f * MAX_DCMPS + ELEM_DCMP_IDX] = nFieldElems;
        decomps[f * MAX_DCMPS + VAL_DCMP_IDX]  = nElemVals;

        cnts[ELEM_CNT_IDX] += nFieldElems;
        cnts[VAL_CNT_IDX]  += nElemVals;

        // Hash the field specification
        std::stringstream hashStream;
        std::size_t nSubFields = fielddefs[f]->m_fields.size();
        for (int sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fielddefs[f]->m_fields[sf];
        }

        nSubFields = fielddefs[f]->m_basis.size();
        for (int sf = 0; sf < nSubFields; ++sf)
        {
            hashStream << fielddefs[f]->m_basis[sf];
        }

        // Determine SHAPE attribute
        std::stringstream shapeStringStream;
        shapeStringStream << ShapeTypeMap[fielddefs[f]->m_shapeType];

        if (fielddefs[f]->m_numHomogeneousDir > 0)
        {
            if (homDim == -1)
            {
                homDim = fielddefs[f]->m_numHomogeneousDir;
            }

            ASSERTL1(homDim == fielddefs[f]->m_numHomogeneousDir,
                     "HDF5 does not support variable homogeneous directions in "
                     "the same file.");

            shapeStringStream << "-HomogenousExp"
                              << fielddefs[f]->m_numHomogeneousDir << "D";
        }

        if (fielddefs[f]->m_homoStrips)
        {
            shapeStringStream << "-Strips";
        }

        shapeStrings[f] = shapeStringStream.str();
        hashStream << shapeStringStream.str();

        // Determine HOMOGENEOUS attributes
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            nSubFields = fielddefs[f]->m_homogeneousLengths.size();
            homoLengths[f].resize(nSubFields);
            for (int sf = 0; sf < nSubFields; ++sf)
            {
                NekDouble len = fielddefs[f]->m_homogeneousLengths[sf];
                hashStream << len;
                homoLengths[f][sf] = len;
            }

            nSubFields = fielddefs[f]->m_homogeneousYIDs.size();
            if (nSubFields > 0)
            {
                homoYIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMY_DCMP_IDX] = nSubFields;
                cnts[HOMY_CNT_IDX] += nSubFields;
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    homoYIDs[f][sf] = fielddefs[f]->m_homogeneousYIDs[sf];
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousZIDs.size();
            if (nSubFields > 0)
            {
                homoZIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMZ_DCMP_IDX] = nSubFields;
                cnts[HOMZ_CNT_IDX] += nSubFields;
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    homoZIDs[f][sf] = fielddefs[f]->m_homogeneousZIDs[sf];
                }
            }

            nSubFields = fielddefs[f]->m_homogeneousSIDs.size();
            if (nSubFields > 0)
            {
                homoSIDs[f].resize(nSubFields);
                decomps[f * MAX_DCMPS + HOMS_DCMP_IDX] = nSubFields;
                cnts[HOMS_CNT_IDX] += nSubFields;
                for (int sf = 0; sf < nSubFields; ++sf)
                {
                    homoSIDs[f][sf] = fielddefs[f]->m_homogeneousSIDs[sf];
                }
            }
        }

        if (fielddefs[f]->m_uniOrder)
        {
            std::vector<unsigned int> elemModes(fielddefs[f]->m_basis.size());

            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_basis.size(); ++i)
            {
                elemModes[i] = fielddefs[f]->m_numModes[i];
            }

            if (varOrder)
            {
                for (std::vector<int>::size_type i = 0; i < nFieldElems; ++i)
                {
                    std::copy(elemModes.begin(), elemModes.end(),
                              std::back_inserter(numModesPerDirVar[f]));
                }
                decomps[f * MAX_DCMPS + ORDER_DCMP_IDX] =
                    nFieldElems * elemModes.size();
                cnts[ORDER_CNT_IDX] += nFieldElems * elemModes.size();
            }
            else
            {
                std::stringstream numModesStringStream;
                numModesStringStream << "UNIORDER:";
                for (std::vector<int>::size_type i = 0;
                     i < elemModes.size(); i++)
                {
                    if (i > 0)
                    {
                        numModesStringStream << ",";
                    }
                    numModesStringStream << elemModes[i];
                }

                numModesPerDirUni[f] = numModesStringStream.str();
                hashStream << numModesPerDirUni[f];
            }
        }
        else
        {
            numModesPerDirVar[f] = fielddefs[f]->m_numModes;
            decomps[f * MAX_DCMPS + ORDER_DCMP_IDX] =
                fielddefs[f]->m_numModes.size();
            cnts[ORDER_CNT_IDX] += fielddefs[f]->m_numModes.size();
        }

        std::hash<std::string> string_hasher;
        std::stringstream fieldNameStream;
        uint64_t fieldDefHash = string_hasher(hashStream.str());

        decomps[f * MAX_DCMPS + HASH_DCMP_IDX]         = fieldDefHash;
        all_hashes[m_comm->GetRank() * nMaxFields + f] = fieldDefHash;

        fieldNameStream << fieldDefHash;
        fieldNames[f] = fieldNameStream.str();
    }

    // Gather information from all MPI processes
    std::vector<uint64_t> all_cnts = m_comm->Gather(root_rank, cnts);
    std::vector<uint64_t> all_idxs(m_comm->GetSize() * MAX_IDXS, 0);
    std::vector<uint64_t> all_decomps = m_comm->Gather(root_rank, decomps);
    std::vector<uint64_t> all_dsetsize(MAX_CNTS, 0);

    // The root rank creates the file layout from scratch
    if (amRoot)
    {
        H5::FileSharedPtr outfile = H5::File::Create(outFile, H5F_ACC_TRUNC);
        ASSERTL1(outfile, prfx.str() + "cannot create HDF5 file.");
        H5::GroupSharedPtr root = outfile->CreateGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot create root group.");
        TagWriterSharedPtr info_writer(new H5TagWriter(root));
        AddInfoTag(info_writer, fieldmetadatamap);

        // Record file format version as attribute in main group.
        root->SetAttribute("FORMAT_VERSION", FORMAT_VERSION);

        // Calculate the indexes to be used by each MPI process when reading the
        // IDS and DATA datasets
        std::size_t nTotElems = 0, nTotVals = 0, nTotOrder = 0;
        std::size_t nTotHomY = 0, nTotHomZ = 0, nTotHomS = 0;
        int nRanks = m_comm->GetSize();
        for (int r = 0; r < nRanks; ++r)
        {
            all_idxs[r * MAX_IDXS + IDS_IDX_IDX]   = nTotElems;
            all_idxs[r * MAX_IDXS + DATA_IDX_IDX]  = nTotVals;
            all_idxs[r * MAX_IDXS + ORDER_IDX_IDX] = nTotOrder;
            all_idxs[r * MAX_IDXS + HOMY_IDX_IDX]  = nTotHomY;
            all_idxs[r * MAX_IDXS + HOMZ_IDX_IDX]  = nTotHomZ;
            all_idxs[r * MAX_IDXS + HOMS_IDX_IDX]  = nTotHomS;

            nTotElems += all_cnts[r * MAX_CNTS + ELEM_CNT_IDX];
            nTotVals  += all_cnts[r * MAX_CNTS + VAL_CNT_IDX];
            nTotOrder += all_cnts[r * MAX_CNTS + ORDER_CNT_IDX];
            nTotHomY  += all_cnts[r * MAX_CNTS + HOMY_CNT_IDX];
            nTotHomZ  += all_cnts[r * MAX_CNTS + HOMZ_CNT_IDX];
            nTotHomS  += all_cnts[r * MAX_CNTS + HOMS_CNT_IDX];
        }

        all_dsetsize[ELEM_CNT_IDX ] = nTotElems;
        all_dsetsize[VAL_CNT_IDX  ] = nTotVals;
        all_dsetsize[ORDER_CNT_IDX] = nTotOrder;
        all_dsetsize[HOMY_CNT_IDX ] = nTotHomY;
        all_dsetsize[HOMZ_CNT_IDX ] = nTotHomZ;
        all_dsetsize[HOMS_CNT_IDX ] = nTotHomS;

        // Create DECOMPOSITION dataset: basic field info for each MPI process
        H5::DataTypeSharedPtr decomps_type =
            H5::DataType::OfObject(all_decomps[0]);
        H5::DataSpaceSharedPtr decomps_space =
            H5::DataSpace::OneD(all_decomps.size());
        H5::DataSetSharedPtr decomps_dset =
            root->CreateDataSet("DECOMPOSITION", decomps_type, decomps_space);
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot create DECOMPOSITION dataset.");

        // Create IDS dataset: element ids
        H5::DataTypeSharedPtr ids_type =
            H5::DataType::OfObject(fielddefs[0]->m_elementIDs[0]);
        H5::DataSpaceSharedPtr ids_space = H5::DataSpace::OneD(nTotElems);
        H5::DataSetSharedPtr ids_dset =
            root->CreateDataSet("ELEMENTIDS", ids_type, ids_space);
        ASSERTL1(ids_dset, prfx.str() + "cannot create ELEMENTIDS dataset.");

        // Create DATA dataset: element data
        H5::DataTypeSharedPtr data_type =
            H5::DataType::OfObject(fielddata[0][0]);
        H5::DataSpaceSharedPtr data_space = H5::DataSpace::OneD(nTotVals);
        H5::DataSetSharedPtr data_dset =
            root->CreateDataSet("DATA", data_type, data_space);
        ASSERTL1(data_dset, prfx.str() + "cannot create DATA dataset.");

        // Create HOMOGENEOUSYIDS dataset: homogeneous y-plane IDs
        if (nTotHomY > 0)
        {
            H5::DataTypeSharedPtr homy_type =
                H5::DataType::OfObject(homoYIDs[0][0]);
            H5::DataSpaceSharedPtr homy_space = H5::DataSpace::OneD(nTotHomY);
            H5::DataSetSharedPtr homy_dset =
                root->CreateDataSet("HOMOGENEOUSYIDS", homy_type, homy_space);
            ASSERTL1(homy_dset,
                     prfx.str() + "cannot create HOMOGENEOUSYIDS dataset.");
        }

        // Create HOMOGENEOUSYIDS dataset: homogeneous z-plane IDs
        if (nTotHomZ > 0)
        {
            H5::DataTypeSharedPtr homz_type =
                H5::DataType::OfObject(homoZIDs[0][0]);
            H5::DataSpaceSharedPtr homz_space = H5::DataSpace::OneD(nTotHomZ);
            H5::DataSetSharedPtr homz_dset =
                root->CreateDataSet("HOMOGENEOUSZIDS", homz_type, homz_space);
            ASSERTL1(homz_dset,
                     prfx.str() + "cannot create HOMOGENEOUSZIDS dataset.");
        }

        // Create HOMOGENEOUSSIDS dataset: homogeneous strip IDs
        if (nTotHomS > 0)
        {
            H5::DataTypeSharedPtr homs_type =
                H5::DataType::OfObject(homoSIDs[0][0]);
            H5::DataSpaceSharedPtr homs_space = H5::DataSpace::OneD(nTotHomS);
            H5::DataSetSharedPtr homs_dset =
                root->CreateDataSet("HOMOGENEOUSSIDS", homs_type, homs_space);
            ASSERTL1(homs_dset,
                     prfx.str() + "cannot create HOMOGENEOUSSIDS dataset.");
        }

        // Create POLYORDERS dataset: elemental polynomial orders
        if (varOrder)
        {
            H5::DataTypeSharedPtr order_type =
                H5::DataType::OfObject(numModesPerDirVar[0][0]);
            H5::DataSpaceSharedPtr order_space = H5::DataSpace::OneD(nTotOrder);
            H5::DataSetSharedPtr order_dset =
                root->CreateDataSet("POLYORDERS", order_type, order_space);
            ASSERTL1(order_dset,
                     prfx.str() + "cannot create POLYORDERS dataset.");
        }
    }

    m_comm->Bcast(all_dsetsize, root_rank);

    // Datasets, root group and HDF5 file are all closed automatically since
    // they are now out of scope. Now we need to determine which process will
    // write the group representing the field description in the HDF5 file. This
    // next block of code performs this by finding all unique hashes and then
    // determining one process that will create (possibly more than one) group
    // for that hash. An alternative would be to communicate the field
    // information to the root processor, but this is a bit convoluted.

    // This set stores the unique hashes.
    std::set<uint64_t> hashToProc;
    // This map takes ranks to hashes this process will write.
    std::map<int, std::vector<uint64_t> > writingProcs;

    // Gather all field hashes to every processor.
    m_comm->AllReduce(all_hashes, LibUtilities::ReduceMax);

    for (int n = 0; n < m_comm->GetSize(); ++n)
    {
        for (int i = 0; i < nMaxFields; ++i)
        {
            uint64_t hash = all_hashes[n*nMaxFields + i];

            // Note hash can be zero if, on this process, nFields < nMaxFields.
            if (hashToProc.find(hash) != hashToProc.end() || hash == 0)
            {
                continue;
            }
            hashToProc.insert(hash);
            writingProcs[n].push_back(hash);
        }
    }

    // Having constructed the map, go ahead and write the attributes out.
    for (auto &sIt : writingProcs)
    {
        int rank = sIt.first;

        // Write out this rank's groups.
        if (m_comm->GetRank() == rank)
        {
            H5::PListSharedPtr serialProps = H5::PList::Default();
            H5::PListSharedPtr writeSR     = H5::PList::Default();

            // Reopen the file
            H5::FileSharedPtr outfile =
                H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
            ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
            H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
            ASSERTL1(root, prfx.str() + "cannot open root group.");

            // Write a HDF5 group for each field
            hashToProc.clear();
            for (int i = 0; i < sIt.second.size(); ++i)
            {
                for (int f = 0; f < nFields; ++f)
                {
                    if (sIt.second[i] !=
                        all_hashes[m_comm->GetRank() * nMaxFields + f] ||
                        hashToProc.find(sIt.second[i]) != hashToProc.end())
                    {
                        continue;
                    }

                    hashToProc.insert(sIt.second[i]);

                    // Just in case we've already written this
                    H5::GroupSharedPtr field_group =
                        root->CreateGroup(fieldNames[f]);
                    ASSERTL1(field_group,
                             prfx.str() + "cannot create field group.");
                    field_group->SetAttribute("FIELDS", fielddefs[f]->m_fields);
                    field_group->SetAttribute("BASIS", fielddefs[f]->m_basis);
                    field_group->SetAttribute("SHAPE", shapeStrings[f]);

                    if (homoLengths[f].size() > 0)
                    {
                        field_group->SetAttribute("HOMOGENEOUSLENGTHS",
                                                  homoLengths[f]);
                    }

                    // If the field has only uniform order, we write the order
                    // into the NUMMODESPERDIR attribute. Otherwise, we'll go
                    // ahead and assume everything is mixed and fix this in the
                    // read later if required.
                    if (!varOrder)
                    {
                        field_group->SetAttribute("NUMMODESPERDIR",
                                                  numModesPerDirUni[f]);
                    }
                    else
                    {
                        std::string numModesPerDir = "MIXORDER";
                        field_group->SetAttribute("NUMMODESPERDIR",
                                                  numModesPerDir);
                    }
                }
            }
        }

        // We block to avoid more than one processor opening the file at a time.
        m_comm->Block();
    }

    // Write the DECOMPOSITION dataset
    if (amRoot)
    {
        H5::PListSharedPtr serialProps = H5::PList::Default();
        H5::PListSharedPtr writeSR     = H5::PList::Default();

        // Reopen the file
        H5::FileSharedPtr outfile =
            H5::File::Open(outFile, H5F_ACC_RDWR, serialProps);
        ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
        H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
        ASSERTL1(root, prfx.str() + "cannot open root group.");

        // Write the DECOMPOSITION dataset
        H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
        ASSERTL1(decomps_dset,
                 prfx.str() + "cannot open DECOMPOSITION dataset.");

        H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
        ASSERTL1(decomps_fspace,
                 prfx.str() + "cannot open DECOMPOSITION filespace.");

        decomps_fspace->SelectRange(0, all_decomps.size());
        decomps_dset->Write(all_decomps, decomps_fspace, writeSR);
    }

    // Initialise the dataset indexes for all MPI processes
    std::vector<uint64_t> idx = m_comm->Scatter(root_rank, all_idxs);
    uint64_t ids_i            = idx[IDS_IDX_IDX];
    uint64_t data_i           = idx[DATA_IDX_IDX];
    uint64_t order_i          = idx[ORDER_IDX_IDX];
    uint64_t homy_i           = idx[HOMY_IDX_IDX];
    uint64_t homz_i           = idx[HOMZ_IDX_IDX];
    uint64_t homs_i           = idx[HOMS_IDX_IDX];

    // Set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr writePL = H5::PList::Default();
    if (m_comm->GetSize() > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // Use collective IO
        writePL = H5::PList::DatasetXfer();
        writePL->SetDxMpioCollective();
    }

    // Reopen the file
    H5::FileSharedPtr outfile =
        H5::File::Open(outFile, H5F_ACC_RDWR, parallelProps);
    ASSERTL1(outfile, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    m_comm->Block();

    // all HDF5 groups have now been created. Open the IDS dataset and
    // associated data space
    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("ELEMENTIDS");
    ASSERTL1(ids_dset, prfx.str() + "cannot open ELEMENTIDS dataset.");
    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, prfx.str() + "cannot open ELEMENTIDS filespace.");

    // Open the DATA dataset and associated data space
    H5::DataSetSharedPtr data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");
    H5::DataSpaceSharedPtr data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");

    // Open the optional datasets and data spaces.
    H5::DataSetSharedPtr order_dset, homy_dset, homz_dset, homs_dset;
    H5::DataSpaceSharedPtr order_fspace , homy_fspace, homz_fspace, homs_fspace;

    if (all_dsetsize[ORDER_CNT_IDX])
    {
        order_dset = root->OpenDataSet("POLYORDERS");
        ASSERTL1(order_dset, prfx.str() + "cannot open POLYORDERS dataset.");
        order_fspace = order_dset->GetSpace();
        ASSERTL1(order_fspace, prfx.str() + "cannot open POLYORDERS filespace.");
    }

    if (all_dsetsize[HOMY_CNT_IDX])
    {
        homy_dset = root->OpenDataSet("HOMOGENEOUSYIDS");
        ASSERTL1(homy_dset, prfx.str() + "cannot open HOMOGENEOUSYIDS dataset.");
        homy_fspace = homy_dset->GetSpace();
        ASSERTL1(homy_fspace, prfx.str() + "cannot open HOMOGENEOUSYIDS filespace.");
    }

    if (all_dsetsize[HOMZ_CNT_IDX])
    {
        homz_dset   = root->OpenDataSet("HOMOGENEOUSZIDS");
        ASSERTL1(homz_dset, prfx.str() + "cannot open HOMOGENEOUSZIDS dataset.");
        homz_fspace = homz_dset->GetSpace();
        ASSERTL1(homz_fspace, prfx.str() + "cannot open HOMOGENEOUSZIDS filespace.");
    }

    if (all_dsetsize[HOMS_CNT_IDX])
    {
        homs_dset   = root->OpenDataSet("HOMOGENEOUSSIDS");
        ASSERTL1(homs_dset, prfx.str() + "cannot open HOMOGENEOUSSIDS dataset.");
        homs_fspace = homs_dset->GetSpace();
        ASSERTL1(homs_fspace, prfx.str() + "cannot open HOMOGENEOUSSIDS filespace.");
    }

    // Write the data
    for (int f = 0; f < nFields; ++f)
    {
        // write the element ids
        std::size_t nFieldElems = fielddefs[f]->m_elementIDs.size();
        ids_fspace->SelectRange(ids_i, nFieldElems);
        ids_dset->Write(fielddefs[f]->m_elementIDs, ids_fspace, writePL);
        ids_i += nFieldElems;

        // write the element values
        std::size_t nFieldVals = fielddata[f].size();
        data_fspace->SelectRange(data_i, nFieldVals);
        data_dset->Write(fielddata[f], data_fspace, writePL);
        data_i += nFieldVals;
    }

    if (order_dset)
    {
        for (int f = 0; f < nFields; ++f)
        {
            std::size_t nOrders = numModesPerDirVar[f].size();
            order_fspace->SelectRange(order_i, nOrders);
            order_dset->Write(numModesPerDirVar[f], order_fspace, writePL);
            order_i += nOrders;
        }
    }

    if (homy_dset)
    {
        for (int f = 0; f < nFields; ++f)
        {
            std::size_t nYIDs = homoYIDs[f].size();
            homy_fspace->SelectRange(homy_i, nYIDs);
            homy_dset->Write(homoYIDs[f], homy_fspace, writePL);
            homy_i += nYIDs;
        }
    }

    if (homz_dset)
    {
        for (int f = 0; f < nFields; ++f)
        {
            std::size_t nZIDs = homoZIDs[f].size();
            homz_fspace->SelectRange(homz_i, nZIDs);
            homz_dset->Write(homoZIDs[f], homz_fspace, writePL);
            homz_i += nZIDs;
        }
    }

    if (homs_dset)
    {
        for (int f = 0; f < nFields; ++f)
        {
            std::size_t nSIDs = homoSIDs[f].size();
            homs_fspace->SelectRange(homs_i, nSIDs);
            homs_dset->Write(homoSIDs[f], homs_fspace, writePL);
            homs_i += nSIDs;
        }
    }

    for (int f = nFields; f < nMaxFields; ++f)
    {
        // this MPI process is handling fewer than nMaxFields fields
        // so, since this is a collective operation
        // just rewrite the element ids and values of the last field
        ids_dset->Write(
            fielddefs[nFields - 1]->m_elementIDs, ids_fspace, writePL);
        data_dset->Write(fielddata[nFields - 1], data_fspace, writePL);

        if (order_dset)
        {
            order_dset->Write(numModesPerDirVar[nFields - 1],
                              order_fspace, writePL);
        }

        if (homy_dset)
        {
            homy_dset->Write(homoYIDs[nFields - 1], homy_fspace, writePL);
        }

        if (homz_dset)
        {
            homz_dset->Write(homoZIDs[nFields - 1], homz_fspace, writePL);
        }

        if (homs_dset)
        {
            homs_dset->Write(homoSIDs[nFields - 1], homs_fspace, writePL);
        }
    }

    m_comm->Block();

    // all data has been written
    if (m_comm->GetRank() == 0)
    {
        tm1 = m_comm->Wtime();
        std::cout << " (" << tm1 - tm0 << "s, HDF5)" << std::endl;
    }
}

/**
 * @brief Import a HDF5 format file.
 *
 * @param finfilename       Input filename
 * @param fielddefs         Field definitions of resulting field
 * @param fielddata         Field data of resulting field
 * @param fieldinfomap      Field metadata of resulting field
 * @param ElementIDs        If specified, contains the list of element IDs on
 *                          this rank. The resulting field definitions will only
 *                          contain data for the element IDs specified in this
 *                          array.
 */
void FieldIOHdf5::v_Import(const std::string &infilename,
                           std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                           std::vector<std::vector<NekDouble> > &fielddata,
                           FieldMetaDataMap &fieldinfomap,
                           const Array<OneD, int> &ElementIDs)
{
    std::stringstream prfx;
    int nRanks = m_comm->GetSize();

    // Set properties for parallel file access (if we're in parallel)
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    H5::PListSharedPtr readPL = H5::PList::Default();
    H5::PListSharedPtr readPLInd = H5::PList::Default();

    if (nRanks > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(m_comm);
        // Use collective IO
        readPL = H5::PList::DatasetXfer();
        readPL->SetDxMpioCollective();
        readPLInd = H5::PList::DatasetXfer();
        readPLInd->SetDxMpioIndependent();
    }

    DataSourceSharedPtr dataSource = H5DataSource::create(
        infilename, parallelProps);

    // Open the root group of the hdf5 file
    H5DataSourceSharedPtr h5 =
        std::static_pointer_cast<H5DataSource>(dataSource);
    ASSERTL1(h5, prfx.str() + "cannot open HDF5 file.");
    H5::GroupSharedPtr root = h5->Get()->OpenGroup("NEKTAR");
    ASSERTL1(root, prfx.str() + "cannot open root group.");

    // Check format version
    unsigned int formatVersion;
    H5::Group::AttrIterator attrIt  = root->attr_begin();
    H5::Group::AttrIterator attrEnd = root->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        if (*attrIt == "FORMAT_VERSION")
        {
            break;
        }
    }

    ASSERTL0(attrIt != attrEnd,
             "Unable to determine Nektar++ HDF5 file version");
    root->GetAttribute("FORMAT_VERSION", formatVersion);

    ASSERTL0(formatVersion <= FORMAT_VERSION,
             "File format if " + infilename + " is higher than supported in "
             "this version of Nektar++");

    // Open the datasets
    H5::DataSetSharedPtr decomps_dset = root->OpenDataSet("DECOMPOSITION");
    ASSERTL1(decomps_dset, prfx.str() + "cannot open DECOMPOSITION dataset.");

    H5::DataSetSharedPtr ids_dset = root->OpenDataSet("ELEMENTIDS");
    ASSERTL1(ids_dset, prfx.str() + "cannot open ELEMENTIDS dataset.");

    H5::DataSetSharedPtr data_dset = root->OpenDataSet("DATA");
    ASSERTL1(data_dset, prfx.str() + "cannot open DATA dataset.");

    // Open the dataset file spaces
    H5::DataSpaceSharedPtr decomps_fspace = decomps_dset->GetSpace();
    ASSERTL1(decomps_fspace,
             prfx.str() + "cannot open DECOMPOSITION filespace.");

    H5::DataSpaceSharedPtr ids_fspace = ids_dset->GetSpace();
    ASSERTL1(ids_fspace, prfx.str() + "cannot open ELEMENTIDS filespace.");

    H5::DataSpaceSharedPtr data_fspace = data_dset->GetSpace();
    ASSERTL1(data_fspace, prfx.str() + "cannot open DATA filespace.");

    // Read entire IDS data set to get list of global IDs.
    std::vector<uint64_t> ids;

    ids_dset->Read(ids, ids_fspace, readPL);

    std::unordered_set<uint64_t> toread;
    if (ElementIDs != NullInt1DArray)
    {
        for (uint64_t i = 0; i < ElementIDs.size(); ++i)
        {
            toread.insert(ElementIDs[i]);
        }
    }

    std::vector<uint64_t> decomps;
    decomps_dset->Read(decomps, decomps_fspace, readPL);

    size_t nDecomps = decomps.size() / MAX_DCMPS;
    size_t cnt = 0, cnt2 = 0;

    // Mapping from each decomposition to offsets in the data array.
    std::vector<OffsetHelper> decompsToOffsets (nDecomps);

    // Mapping from each group's hash to a vector of element IDs. Note this has
    // to be unsigned int, since that's what we use in FieldDefinitions.
    std::map<uint64_t, std::vector<unsigned int> > groupsToElmts;

    // Mapping from each group's hash to each of the decompositions.
    std::map<uint64_t, std::set<uint64_t> > groupsToDecomps;

    // True if we are pulling element IDs from ElementIDs.
    bool selective = toread.size() > 0;

    // Counters for data offsets
    OffsetHelper running;

    for (size_t i = 0; i < nDecomps; ++i, cnt += MAX_DCMPS)
    {
        uint64_t nElmt     = decomps[cnt + ELEM_DCMP_IDX];
        uint64_t groupHash = decomps[cnt + HASH_DCMP_IDX];

        std::vector<uint64_t> tmp;

        if (selective)
        {
            for (size_t j = 0; j < nElmt; ++j)
            {
                uint64_t elmtId = ids[cnt2 + j];
                if (toread.find(elmtId) != toread.end())
                {
                    tmp.push_back(elmtId);
                }
            }
        }
        else
        {
            tmp.insert(
                tmp.begin(), ids.begin() + cnt2, ids.begin() + cnt2 + nElmt);
        }

        std::vector<unsigned int> tmp2(nElmt);
        for (size_t j = 0; j < nElmt; ++j)
        {
            tmp2[j] = ids[cnt2+j];
        }

        cnt2 += nElmt;

        if (tmp.size() > 0)
        {
            groupsToDecomps[groupHash].insert(i);
        }

        groupsToElmts[i] = tmp2;
        decompsToOffsets[i] = running;

        running.data  += decomps[cnt + VAL_DCMP_IDX];
        running.order += decomps[cnt + ORDER_DCMP_IDX];
        running.homy  += decomps[cnt + HOMY_DCMP_IDX];
        running.homz  += decomps[cnt + HOMZ_DCMP_IDX];
        running.homs  += decomps[cnt + HOMS_DCMP_IDX];
    }

    for (auto &gIt : groupsToDecomps)
    {
        // Select region from dataset for this decomposition.
        for (auto &sIt : gIt.second)
        {
            std::stringstream fieldNameStream;
            fieldNameStream << gIt.first;

            FieldDefinitionsSharedPtr fielddef =
                MemoryManager<FieldDefinitions>::AllocateSharedPtr();
            ImportFieldDef(readPLInd, root, decomps, sIt, decompsToOffsets[sIt],
                           fieldNameStream.str(), fielddef);

            fielddef->m_elementIDs = groupsToElmts[sIt];
            fielddefs.push_back(fielddef);

            if (fielddata != NullVectorNekDoubleVector)
            {
                std::vector<NekDouble> decompFieldData;
                ImportFieldData(
                    readPLInd, data_dset, data_fspace,
                    decompsToOffsets[sIt].data, decomps, sIt, fielddef,
                    decompFieldData);
                fielddata.push_back(decompFieldData);
            }
        }
    }

    ImportHDF5FieldMetaData(dataSource, fieldinfomap);
    m_comm->Block();
}

/**
 * @brief Import field definitions from a HDF5 document.
 *
 * @param readPL       Reading parameter list.
 * @param root         Root group containing field definitions.
 * @param group        Group name to process.
 * @param def          On output contains field definitions.
 */
void FieldIOHdf5::ImportFieldDef(
    H5::PListSharedPtr        readPL,
    H5::GroupSharedPtr        root,
    std::vector<uint64_t>    &decomps,
    uint64_t                  decomp,
    OffsetHelper              offset,
    std::string               group,
    FieldDefinitionsSharedPtr def)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldDefsHdf5(): ";

    H5::GroupSharedPtr field = root->OpenGroup(group);
    ASSERTL1(field, prfx.str() + "cannot open field group, " + group + '.');

    def->m_uniOrder = false;

    H5::Group::AttrIterator attrIt  = field->attr_begin();
    H5::Group::AttrIterator attrEnd = field->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        const std::string &attrName = *attrIt;
        if (attrName == "FIELDS")
        {
            field->GetAttribute(attrName, def->m_fields);
        }
        else if (attrName == "SHAPE")
        {
            std::string shapeString;
            field->GetAttribute(attrName, shapeString);

            // check to see if homogeneous expansion and if so
            // strip down the shapeString definition
            size_t loc;
            //---> this finds the first location of 'n'!
            if (shapeString.find("Strips") != std::string::npos)
            {
                def->m_homoStrips = true;
            }

            if ((loc = shapeString.find_first_of("-")) != std::string::npos)
            {
                if (shapeString.find("Exp1D") != std::string::npos)
                {
                    def->m_numHomogeneousDir = 1;
                }
                else // HomogeneousExp1D
                {
                    def->m_numHomogeneousDir = 2;
                }

                shapeString.erase(loc, shapeString.length());
            }

            // get the geometrical shape
            bool valid = false;
            for (unsigned int j = 0; j < SIZE_ShapeType; j++)
            {
                if (ShapeTypeMap[j] == shapeString)
                {
                    def->m_shapeType = (ShapeType)j;
                    valid = true;
                    break;
                }
            }

            ASSERTL0(valid, prfx.str() + std::string(
                         "unable to correctly parse the shape type: ")
                     .append(shapeString).c_str());
        }
        else if (attrName == "BASIS")
        {
            field->GetAttribute(attrName, def->m_basis);
            // check the basis is in range
            for (auto &bIt : def->m_basis)
            {
                BasisType bt = bIt;
                ASSERTL0(bt >= 0 && bt < SIZE_BasisType,
                         prfx.str() +
                         "unable to correctly parse the basis types.");
            }
        }
        else if (attrName == "HOMOGENEOUSLENGTHS")
        {
            field->GetAttribute(attrName, def->m_homogeneousLengths);
        }
        else if (attrName == "NUMMODESPERDIR")
        {
            std::string numModesPerDir;
            field->GetAttribute(attrName, numModesPerDir);

            if (strstr(numModesPerDir.c_str(), "UNIORDER:"))
            {
                def->m_uniOrder = true;
                bool valid = ParseUtils::GenerateVector(
                    numModesPerDir.substr(9), def->m_numModes);
                ASSERTL0(valid,
                         prfx.str() +
                         "unable to correctly parse the number of modes.");
            }
        }
        else if (attrName == "POINTSTYPE")
        {
            std::string pointsString;
            field->GetAttribute(attrName, pointsString);
            def->m_pointsDef = true;

            std::vector<std::string> pointsStrings;
            bool valid = ParseUtils::GenerateVector(
                pointsString, pointsStrings);
            ASSERTL0(valid,
                     prfx.str() +
                     "unable to correctly parse the points types.");
            for (std::vector<std::string>::size_type i = 0;
                 i < pointsStrings.size();
                 i++)
            {
                valid = false;
                for (unsigned int j = 0; j < SIZE_PointsType; j++)
                {
                    if (kPointsTypeStr[j] == pointsStrings[i])
                    {
                        def->m_points.push_back((PointsType)j);
                        valid = true;
                        break;
                    }
                }

                ASSERTL0(
                    valid,
                    prfx.str() +
                    std::string(
                        "unable to correctly parse the points type: ")
                    .append(pointsStrings[i])
                    .c_str());
            }
        }
        else if (attrName == "NUMPOINTSPERDIR")
        {
            std::string numPointsPerDir;
            field->GetAttribute(attrName, numPointsPerDir);
            def->m_numPointsDef = true;

            bool valid = ParseUtils::GenerateVector(
                numPointsPerDir, def->m_numPoints);
            ASSERTL0(valid,
                     prfx.str() +
                     "unable to correctly parse the number of points.");
        }
        else
        {
            std::string errstr("unknown attribute: ");
            errstr += attrName;
            ASSERTL1(false, prfx.str() + errstr.c_str());
        }
    }

    if (def->m_numHomogeneousDir >= 1)
    {
        H5::DataSetSharedPtr homz_dset = root->OpenDataSet("HOMOGENEOUSZIDS");
        H5::DataSpaceSharedPtr homz_fspace = homz_dset->GetSpace();
        uint64_t dsize = decomps[decomp * MAX_DCMPS + HOMZ_DCMP_IDX];
        homz_fspace->SelectRange(offset.homz, dsize);
        homz_dset->Read(def->m_homogeneousZIDs, homz_fspace, readPL);
    }

    if (def->m_numHomogeneousDir >= 2)
    {
        H5::DataSetSharedPtr homy_dset = root->OpenDataSet("HOMOGENEOUSYIDS");
        H5::DataSpaceSharedPtr homy_fspace = homy_dset->GetSpace();
        uint64_t dsize = decomps[decomp * MAX_DCMPS + HOMY_DCMP_IDX];
        homy_fspace->SelectRange(offset.homy, dsize);
        homy_dset->Read(def->m_homogeneousYIDs, homy_fspace, readPL);
    }

    if (def->m_homoStrips)
    {
        H5::DataSetSharedPtr homs_dset = root->OpenDataSet("HOMOGENEOUSSIDS");
        H5::DataSpaceSharedPtr homs_fspace = homs_dset->GetSpace();
        uint64_t dsize = decomps[decomp * MAX_DCMPS + HOMS_DCMP_IDX];
        homs_fspace->SelectRange(offset.homs, dsize);
        homs_dset->Read(def->m_homogeneousSIDs, homs_fspace, readPL);
    }

    if (!def->m_uniOrder)
    {
        H5::DataSetSharedPtr order_dset = root->OpenDataSet("POLYORDERS");
        H5::DataSpaceSharedPtr order_fspace = order_dset->GetSpace();
        uint64_t dsize = decomps[decomp * MAX_DCMPS + ORDER_DCMP_IDX];
        order_fspace->SelectRange(offset.order, dsize);
        order_dset->Read(def->m_numModes, order_fspace, readPL);
    }
}

/**
 * @brief Import the field data from the HDF5 document.
 *
 * @param readPL       Reading parameter list.
 * @param data_dset    Pointer to the `DATA` dataset.
 * @param data_fspace  Pointer to the `DATA` data space.
 * @param data_i       Index in the `DATA` dataset to start reading from.
 * @param decomps      Information from the `DECOMPOSITION` dataset.
 * @param decomp       Index of the decomposition.
 * @param fielddef     Field definitions for this file
 * @param fielddata    On return contains resulting field data.
 */
void FieldIOHdf5::ImportFieldData(
    H5::PListSharedPtr               readPL,
    H5::DataSetSharedPtr             data_dset,
    H5::DataSpaceSharedPtr           data_fspace,
    uint64_t                         data_i,
    std::vector<uint64_t>           &decomps,
    uint64_t                         decomp,
    const FieldDefinitionsSharedPtr  fielddef,
    std::vector<NekDouble>          &fielddata)
{
    std::stringstream prfx;
    prfx << m_comm->GetRank() << ": FieldIOHdf5::ImportFieldData(): ";

    uint64_t nElemVals  = decomps[decomp * MAX_DCMPS + VAL_DCMP_IDX];
    uint64_t nFieldVals = nElemVals;

    data_fspace->SelectRange(data_i, nFieldVals);
    data_dset->Read(fielddata, data_fspace, readPL);
    int datasize = CheckFieldDefinition(fielddef);
    ASSERTL0(
        fielddata.size() == datasize * fielddef->m_fields.size(),
        prfx.str() +
        "input data is not the same length as header information.");
}

/**
 * @brief Import field metadata from @p filename and return the data source
 * which wraps @p filename.
 *
 * @param filename          Input filename.
 * @param fieldmetadatamap  Resulting field metadata from @p dataSource.
 */
DataSourceSharedPtr FieldIOHdf5::v_ImportFieldMetaData(
    const std::string &filename,
    FieldMetaDataMap  &fieldmetadatamap)
{
    H5::PListSharedPtr parallelProps = H5::PList::Default();
    DataSourceSharedPtr ans = H5DataSource::create(filename, parallelProps);
    ImportHDF5FieldMetaData(ans, fieldmetadatamap);
    return ans;
}

/**
 * @brief Import field metadata from @p dataSource.
 *
 * @param dataSource        Input datasource, which should be a H5DataSource.
 * @param fieldmetadatamap  Resulting field metadata from @p dataSource.
 */
void FieldIOHdf5::ImportHDF5FieldMetaData(DataSourceSharedPtr dataSource,
                                          FieldMetaDataMap &fieldmetadatamap)
{
    H5DataSourceSharedPtr hdf =
        std::static_pointer_cast<H5DataSource>(dataSource);

    H5::GroupSharedPtr master = hdf->Get()->OpenGroup("NEKTAR");
    // New metadata format only in HDF
    H5::GroupSharedPtr metadata = master->OpenGroup("Metadata");

    if (metadata)
    {
        H5::Group::AttrIterator param = metadata->attr_begin(),
                                pEnd = metadata->attr_end();
        for (; param != pEnd; ++param)
        {
            std::string paramString = *param;
            if (paramString != "Provenance")
            {
                std::string paramBodyStr;
                metadata->GetAttribute(paramString, paramBodyStr);
                fieldmetadatamap[paramString] = paramBodyStr;
            }
        }
    }
}

}
}
