///////////////////////////////////////////////////////////////////////////////
//
// File FieldIOHdf5.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Field IO to/from HDF5
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOHDF5_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOHDF5_H

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/H5.h>

namespace Nektar
{
namespace LibUtilities
{

namespace H5
{
class Group;
typedef std::shared_ptr<Group> GroupSharedPtr;
}

/**
 * @class Class encapsulating simple HDF5 data source using H5 reader utilities.
 */
class H5DataSource : public DataSource
{
public:
    /// Constructor based on filename.
    H5DataSource(const std::string &fn, H5::PListSharedPtr parallelProps)
        : doc(H5::File::Open(fn, H5F_ACC_RDONLY, parallelProps))
    {
    }

    /// Get H5::FileSharedPtr reference to file.
    H5::FileSharedPtr Get()
    {
        return doc;
    }

    /// Get H5::FileSharedPtr reference to file.
    const H5::FileSharedPtr Get() const
    {
        return doc;
    }

    /// Static constructor for this data source.
    static DataSourceSharedPtr create(
        const std::string &fn, H5::PListSharedPtr parallelProps)
    {
        return DataSourceSharedPtr(new H5DataSource(fn, parallelProps));
    }

private:
    /// HDF5 document.
    H5::FileSharedPtr doc;
};
typedef std::shared_ptr<H5DataSource> H5DataSourceSharedPtr;

/**
 * @class Simple class for writing hierarchical data using HDF5.
 */
class H5TagWriter : public TagWriter
{
public:
    /// Default constructor.
    H5TagWriter(H5::GroupSharedPtr grp) : m_Group(grp) {}

    /// Add a child node.
    TagWriterSharedPtr AddChild(const std::string &name)
    {
        H5::GroupSharedPtr child = m_Group->CreateGroup(name);
        return TagWriterSharedPtr(new H5TagWriter(child));
    }

    /// Set an attribute key/value pair on this tag.
    void SetAttr(const std::string &key, const std::string &val)
    {
        m_Group->SetAttribute(key, val);
    }

private:
    /// HDF5 group for this tag.
    H5::GroupSharedPtr m_Group;
};
typedef std::shared_ptr<H5TagWriter> H5TagWriterSharedPtr;

/**
 * @class Class for operating on HDF5-based FLD files.
 *
 * This class implements a HDF5 reader/writer based on MPI/O that is designed to
 * operate on a single file across all processors of a simulation. The
 * definition follows vaguely similar lines to XML output but is stored somewhat
 * differently to accommodate parallel reading and writing. At a basic level
 * metadata is organised as follows:
 *
 *   - Nektar++ data lies in the root `/NEKTAR` group.
 *   - The contents of a FieldDefinitions object is hashed to construct a unique
 *     identifier for each object, which becomes the name of a group within the
 *     root group. We then use the H5TagWriter to assign the field definitions
 *     to each group.
 *   - In a similar fashion, we create a `Metadata` group to contain field
 *     metadata that is written.
 *
 * We then define five data sets to contain field data:
 *
 *   - The `DATA` dataset contains the double-precision modal coefficient data.
 *   - The `IDS` dataset contains the element IDs of the elements that are
 *     written out.
 *   - The `POLYORDERS` dataset is written if the field data contains variable
 *     polynomial order, and contains the (possibly hetergeneous) mode orders in
 *     each direction for each of the elements.
 *   - The `HOMOGENEOUSZIDS` dataset contains the IDs of z-planes for
 *     homogeneous simulations, if the data are homogeneous.
 *   - The `HOMOGENEOUSYIDS` dataset contains the IDs of y-planes for
 *     homogeneous simulations, if the data are homogeneous.
 *   - The `HOMOGENEOUSSIDS` dataset contains the strip IDs for
 *     homogeneous simulations, if the data are homogeneous and use strips.
 *
 * The ordering is defined according to the `DECOMPOSITION` dataset. A
 * `decomposition' in this class is essentially a single field definition with
 * its accompanying data. Data are written into each dataset by the order of
 * each decomposition. Each decomposition contains the following seven integers
 * that define it per field definition per processor:
 *
 *   - Number of elements in this field definition (index #ELEM_DCMP_IDX).
 *   - Number of entries in the `DATA` array for this field definition
 *     (index #VAL_DCMP_IDX)
 *   - Number of entries in the `POLYORDERS` array for this field definition
 *     (index #ORDER_DCMP_IDX)
 *   - Number of entries in the `HOMOGENEOUSZIDS` array (index #HOMZ_DCMP_IDX).
 *   - Number of entries in the `HOMOGENEOUSYIDS` array (index #HOMY_DCMP_IDX).
 *   - Number of entries in the `HOMOGENEOUSSIDS` array (index #HOMS_DCMP_IDX).
 *   - Hash of the field definition, represented as a 32-bit integer, which
 *     describes the name of the attribute that contains the rest of the field
 *     definition information (e.g. field names, basis type, etc).
 *
 * The number of decompositions is therefore calculated as the field size
 * divided by #MAX_DCMPS which allows us to calculate the offsets of the data
 * for each field definition within the arrays.
 */
class FieldIOHdf5 : public FieldIO
{
public:
    static const unsigned int FORMAT_VERSION;

    static const unsigned int ELEM_DCMP_IDX;
    static const unsigned int VAL_DCMP_IDX;
    static const unsigned int ORDER_DCMP_IDX;
    static const unsigned int HOMY_DCMP_IDX;
    static const unsigned int HOMZ_DCMP_IDX;
    static const unsigned int HOMS_DCMP_IDX;
    static const unsigned int HASH_DCMP_IDX;
    static const unsigned int MAX_DCMPS;

    static const unsigned int ELEM_CNT_IDX;
    static const unsigned int VAL_CNT_IDX;
    static const unsigned int ORDER_CNT_IDX;
    static const unsigned int HOMY_CNT_IDX;
    static const unsigned int HOMZ_CNT_IDX;
    static const unsigned int HOMS_CNT_IDX;
    static const unsigned int MAX_CNTS;

    static const unsigned int IDS_IDX_IDX;
    static const unsigned int DATA_IDX_IDX;
    static const unsigned int ORDER_IDX_IDX;
    static const unsigned int HOMY_IDX_IDX;
    static const unsigned int HOMZ_IDX_IDX;
    static const unsigned int HOMS_IDX_IDX;
    static const unsigned int MAX_IDXS;

    /// Creates an instance of this class
    LIB_UTILITIES_EXPORT static FieldIOSharedPtr create(
        LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    {
        return MemoryManager<FieldIOHdf5>::AllocateSharedPtr(pComm,
                                                             sharedFilesystem);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT FieldIOHdf5(
        LibUtilities::CommSharedPtr pComm,
        bool sharedFilesystem);

    LIB_UTILITIES_EXPORT virtual ~FieldIOHdf5()
    {
    }

    /// Get class name
    inline virtual const std::string &GetClassName() const
    {
        return className;
    }

private:
    struct OffsetHelper {
        OffsetHelper() : data(0), order(0), homy(0), homz(0), homs(0) {}
        OffsetHelper(const OffsetHelper &in) :
            data(in.data), order(in.order), homy(in.homy), homz(in.homz),
            homs(in.homs)
        {
        }
        OffsetHelper& operator=(const OffsetHelper &in) = default;

        uint64_t data, order, homy, homz, homs;
    };

    LIB_UTILITIES_EXPORT virtual void v_Write(
        const std::string &outFile,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const bool backup = false);

    LIB_UTILITIES_EXPORT virtual void v_Import(
        const std::string &infilename,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata =
                                                      NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const Array<OneD, int> &ElementIDs = NullInt1DArray);

    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        const std::string &filename, FieldMetaDataMap &fieldmetadatamap);

    LIB_UTILITIES_EXPORT void ImportHDF5FieldMetaData(
        DataSourceSharedPtr dataSource, FieldMetaDataMap &fieldmetadatamap);

    LIB_UTILITIES_EXPORT void ImportFieldDef(
        H5::PListSharedPtr        readPL,
        H5::GroupSharedPtr        root,
        std::vector<uint64_t>    &decomps,
        uint64_t                  decomp,
        OffsetHelper              offset,
        std::string               group,
        FieldDefinitionsSharedPtr def);

    LIB_UTILITIES_EXPORT void ImportFieldData(
        H5::PListSharedPtr               readPL,
        H5::DataSetSharedPtr             data_dset,
        H5::DataSpaceSharedPtr           data_fspace,
        uint64_t                         data_i,
        std::vector<uint64_t>           &decomps,
        uint64_t                         decomp,
        const FieldDefinitionsSharedPtr  fielddef,
        std::vector<NekDouble>          &fielddata);
};
}
}
#endif
