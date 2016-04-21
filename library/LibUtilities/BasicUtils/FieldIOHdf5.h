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
typedef boost::shared_ptr<Group> GroupSharedPtr;
}

class H5DataSource : public DataSource
{
    H5::FileSharedPtr doc;

public:
    H5DataSource(const std::string &fn)
        : doc(H5::File::Open(fn, H5F_ACC_RDONLY))
    {
    }

    H5::FileSharedPtr Get()
    {
        return doc;
    }
    const H5::FileSharedPtr Get() const
    {
        return doc;
    }

    static DataSourceSharedPtr create(const std::string &fn)
    {
        return DataSourceSharedPtr(new H5DataSource(fn));
    }
};
typedef boost::shared_ptr<H5DataSource> H5DataSourceSharedPtr;

class H5TagWriter : public TagWriter
{
public:
    H5TagWriter(H5::GroupSharedPtr grp) : m_Group(grp) {}

    TagWriterSharedPtr AddChild(const std::string &name)
    {
        H5::GroupSharedPtr child = m_Group->CreateGroup(name);
        return TagWriterSharedPtr(new H5TagWriter(child));
    }

    void SetAttr(const std::string &key, const std::string &val)
    {
        m_Group->SetAttribute(key, val);
    }

private:
    H5::GroupSharedPtr m_Group;
};
typedef boost::shared_ptr<H5TagWriter> H5TagWriterSharedPtr;

/**
 * @class Class for operating on HDF5-based FLD files.
 */
class FieldIOHdf5 : public FieldIO
{
public:
    static const unsigned int ELEM_DCMP_IDX;
    static const unsigned int VAL_DCMP_IDX;
    static const unsigned int HASH_DCMP_IDX;
    static const unsigned int MAX_DCMPS;
    static const unsigned int ELEM_CNT_IDX;
    static const unsigned int VAL_CNT_IDX;
    static const unsigned int MAX_CNTS;
    static const unsigned int IDS_IDX_IDX;
    static const unsigned int DATA_IDX_IDX;
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

    FieldIOHdf5(LibUtilities::CommSharedPtr pComm, bool sharedFilesystem);

    inline virtual const std::string &GetClassName() const
    {
        return className;
    }

private:
    /// Write data in FLD format
    LIB_UTILITIES_EXPORT virtual void v_Write(
        const std::string &outFile,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap);

    LIB_UTILITIES_EXPORT virtual void v_Import(const std::string &infilename,
                          std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                          std::vector<std::vector<NekDouble> > &fielddata =
                              NullVectorNekDoubleVector,
                          FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
                          const Array<OneD, int> ElementiDs = NullInt1DArray);

    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        std::string filename, FieldMetaDataMap &fieldmetadatamap);
    LIB_UTILITIES_EXPORT void v_ImportFieldMetaData(
        DataSourceSharedPtr dataSource, FieldMetaDataMap &fieldmetadatamap);

    /// Imports the field definitions.
    LIB_UTILITIES_EXPORT void ImportFieldDef(H5::PListSharedPtr        readPL,
                                             H5::GroupSharedPtr        root,
                                             std::string               group,
                                             FieldDefinitionsSharedPtr def);

    /// Imports the data fields.
    LIB_UTILITIES_EXPORT void ImportFieldData(
        H5::PListSharedPtr               readPL,
        H5::DataSetSharedPtr             data_dset,
        H5::DataSpaceSharedPtr           data_fspace,
        size_t                           data_i,
        std::vector<std::size_t>        &decomps,
        size_t                           decomp,
        const FieldDefinitionsSharedPtr  fielddef,
        std::vector<NekDouble>          &fielddata);
};
}
}
#endif
