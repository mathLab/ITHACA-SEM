////////////////////////////////////////////////////////////////////////////////
//
//  File: H5.cpp
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
//  Description: Minimal HDF5 wrapper
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/H5.h>
#include <LibUtilities/Foundations/BasisType.h>

#ifdef NEKTAR_USE_MPI
#include <LibUtilities/Communication/CommMpi.h>
#endif

namespace Nektar
{
namespace LibUtilities
{
namespace H5
{

Object::Object() : m_Id(H5I_INVALID_HID)
{
}
Object::Object(hid_t id) : m_Id(id)
{
}
Object::~Object()
{
}

PList::PList() : Object(H5P_DEFAULT)
{
}
PList::PList(hid_t cls) : Object()
{
    H5_CONSTRUCT(m_Id, H5Pcreate, (cls));
}
PList::~PList()
{
    Close();
}
void PList::Close()
{
    H5_CALL(H5Pclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

/// Default options
PListSharedPtr PList::Default()
{
    return PListSharedPtr(new PList());
}

/// Properties for object creation
PListSharedPtr PList::ObjectCreate()
{
    return PListSharedPtr(new PList(H5P_OBJECT_CREATE));
}

/// Properties for file creation
PListSharedPtr PList::FileCreate()
{
    return PListSharedPtr(new PList(H5P_FILE_CREATE));
}

/// Properties for file access
PListSharedPtr PList::FileAccess()
{
    return PListSharedPtr(new PList(H5P_FILE_ACCESS));
}

/// Properties for dataset creation
PListSharedPtr PList::DatasetCreate()
{
    return PListSharedPtr(new PList(H5P_DATASET_CREATE));
}

/// Properties for dataset access
PListSharedPtr PList::DatasetAccess()
{
    return PListSharedPtr(new PList(H5P_DATASET_ACCESS));
}

/// Properties for raw data transfer
PListSharedPtr PList::DatasetXfer()
{
    return PListSharedPtr(new PList(H5P_DATASET_XFER));
}

/// Properties for file mounting
PListSharedPtr PList::FileMount()
{
    return PListSharedPtr(new PList(H5P_FILE_MOUNT));
}

/// Properties for group creation
PListSharedPtr PList::GroupCreate()
{
    return PListSharedPtr(new PList(H5P_GROUP_CREATE));
}

/// Properties for group access
PListSharedPtr PList::GroupAccess()
{
    return PListSharedPtr(new PList(H5P_GROUP_ACCESS));
}

/// Properties for datatype creation
PListSharedPtr PList::DatatypeCreate()
{
    return PListSharedPtr(new PList(H5P_DATATYPE_CREATE));
}

/// Properties for datatype access
PListSharedPtr PList::DatatypeAccess()
{
    return PListSharedPtr(new PList(H5P_DATATYPE_ACCESS));
}

/// Properties for character encoding when encoding strings or object names
PListSharedPtr PList::StringCreate()
{
    return PListSharedPtr(new PList(H5P_STRING_CREATE));
}

/// Properties for attribute creation
PListSharedPtr PList::AttributeCreate()
{
    return PListSharedPtr(new PList(H5P_ATTRIBUTE_CREATE));
}

/// Properties governing the object copying process
PListSharedPtr PList::ObjectCopy()
{
    return PListSharedPtr(new PList(H5P_OBJECT_COPY));
}

/// Properties governing link creation
PListSharedPtr PList::LinkCreate()
{
    return PListSharedPtr(new PList(H5P_LINK_CREATE));
}

/// Properties governing link traversal when accessing objects
PListSharedPtr PList::LinkAccess()
{
    return PListSharedPtr(new PList(H5P_LINK_ACCESS));
}
void PList::SetChunk(const std::vector<hsize_t> &dims)
{
    H5_CALL(H5Pset_chunk, (m_Id, dims.size(), &dims[0]));
}
void PList::SetDeflate(const unsigned level)
{
    H5_CALL(H5Pset_deflate, (m_Id, level));
}
#ifdef NEKTAR_USE_MPI
void PList::SetDxMpioCollective()
{
    H5_CALL(H5Pset_dxpl_mpio, (m_Id, H5FD_MPIO_COLLECTIVE));
}
void PList::SetDxMpioIndependent()
{
    H5_CALL(H5Pset_dxpl_mpio, (m_Id, H5FD_MPIO_INDEPENDENT));
}
void PList::SetMpio(CommSharedPtr comm)
{
    CommMpiSharedPtr mpi_comm = std::dynamic_pointer_cast<CommMpi>(comm);
    ASSERTL0(mpi_comm, "Can't convert communicator to MPI communicator.")
    // TODO: accept hints
    MPI_Info info = MPI_INFO_NULL;
    H5_CALL(H5Pset_fapl_mpio, (m_Id, mpi_comm->GetComm(), info));
}
#else
void PList::SetDxMpioCollective()
{
    ASSERTL0(false, "Trying to use parallel HDF5 without MPI!");
}
void PList::SetDxMpioIndependent()
{
    ASSERTL0(false, "Trying to use parallel HDF5 without MPI!");
}
void PList::SetMpio(CommSharedPtr comm)
{
    ASSERTL0(false, "Trying to use parallel HDF5 without MPI!");
}
#endif

GroupSharedPtr CanHaveGroupsDataSets::CreateGroup(const std::string &name,
                                                  PListSharedPtr createPL,
                                                  PListSharedPtr accessPL)
{
    hid_t grp;
    H5_CONSTRUCT(grp, H5Gcreate2, (m_Id, name.c_str(), H5P_DEFAULT,
                                   createPL->GetId(), accessPL->GetId()));
    GroupSharedPtr ans(new Group(grp));
    return ans;
}

DataSetSharedPtr CanHaveGroupsDataSets::CreateDataSet(const std::string &name,
                                                      DataTypeSharedPtr type,
                                                      DataSpaceSharedPtr space,
                                                      PListSharedPtr createPL,
                                                      PListSharedPtr accessPL)
{
    hid_t ds;
    H5_CONSTRUCT(ds, H5Dcreate2,
                 (m_Id, name.c_str(), type->GetId(), space->GetId(),
                  H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));

    DataSetSharedPtr ans(new DataSet(ds));
    return ans;
}

// Open an existing group.
// The accessPL can be omitted to use the defaults
GroupSharedPtr CanHaveGroupsDataSets::OpenGroup(const std::string &name,
                                                PListSharedPtr accessPL) const
{
    hid_t grp;
    H5_CONSTRUCT(grp, H5Gopen2, (m_Id, name.c_str(), accessPL->GetId()));
    GroupSharedPtr ans(new Group(grp));
    return ans;
}

// Open an existing dataset
// The accessPL can be omitted to use the defaults
DataSetSharedPtr CanHaveGroupsDataSets::OpenDataSet(
    const std::string &name, PListSharedPtr accessPL) const
{
    hid_t ds;
    H5_CONSTRUCT(ds, H5Dopen2, (m_Id, name.c_str(), accessPL->GetId()));
    DataSetSharedPtr ans(new DataSet(ds));
    return ans;
}

bool CanHaveGroupsDataSets::ContainsDataSet(std::string nm)
{
    for(auto it = begin(); it != end(); ++it)
    {
        if(it.GetName() == nm)
        {
            return true;
        }
    }
    return false;
}

CanHaveGroupsDataSets::LinkIterator CanHaveGroupsDataSets::begin()
{
    // Have to use dynamic because of virtual inheritance
    CanHaveGroupsDataSetsSharedPtr thisSh =
        std::dynamic_pointer_cast<CanHaveGroupsDataSets>(shared_from_this());
    return CanHaveGroupsDataSets::LinkIterator(thisSh);
}

CanHaveGroupsDataSets::LinkIterator CanHaveGroupsDataSets::end()
{
    // Have to use dynamic because of virtual inheritance
    CanHaveGroupsDataSetsSharedPtr thisSh =
        std::dynamic_pointer_cast<CanHaveGroupsDataSets>(shared_from_this());
    return CanHaveGroupsDataSets::LinkIterator(thisSh, GetNumElements());
}

CanHaveGroupsDataSets::LinkIterator::LinkIterator(
    CanHaveGroupsDataSetsSharedPtr grp, hsize_t idx)
    : m_grp(grp), m_idx(-1), m_next(idx), m_size(grp->GetNumElements())
{
    ++*this;
}

const std::string &CanHaveGroupsDataSets::LinkIterator::operator*()
{
    return m_currentName;
}
CanHaveGroupsDataSets::LinkIterator
    &CanHaveGroupsDataSets::LinkIterator::operator++()
{
    m_idx = m_next;
    if (m_idx < m_size)
    {
        H5_CALL(H5Literate, (m_grp->GetId(), H5_INDEX_NAME, H5_ITER_NATIVE,
                             &m_next, LinkIterator::helper, this));
    }
    return *this;
}
bool CanHaveGroupsDataSets::LinkIterator::operator==(
    const CanHaveGroupsDataSets::LinkIterator &other) const
{
    return (m_grp == other.m_grp && m_idx == other.m_idx);
}
herr_t CanHaveGroupsDataSets::LinkIterator::helper(hid_t g_id, const char *name,
                                                   const H5L_info_t * info,
                                                   void *op_data)
{
    boost::ignore_unused(g_id, info);
    CanHaveGroupsDataSets::LinkIterator *iter =
        static_cast<CanHaveGroupsDataSets::LinkIterator *>(op_data);
    iter->m_currentName = name;
    return 1;
}

AttributeSharedPtr CanHaveAttributes::CreateAttribute(const std::string &name,
                                                      DataTypeSharedPtr type,
                                                      DataSpaceSharedPtr space)
{
    return Attribute::Create(m_Id, name, type, space);
}

AttributeSharedPtr CanHaveAttributes::OpenAttribute(const std::string &name)
{
    return Attribute::Open(m_Id, name);
}

int CanHaveAttributes::GetNumAttr() const
{
    H5O_info_t info;
    H5_CALL(H5Oget_info, (m_Id, &info));
    return info.num_attrs;
}

CanHaveAttributes::AttrIterator CanHaveAttributes::attr_begin()
{
    // Have to use dynamic because of virtual inheritance
    CanHaveAttributesSharedPtr thisSh =
        std::dynamic_pointer_cast<CanHaveAttributes>(shared_from_this());
    return CanHaveAttributes::AttrIterator(thisSh);
}
CanHaveAttributes::AttrIterator CanHaveAttributes::attr_end()
{
    // Have to use dynamic because of virtual inheritance
    CanHaveAttributesSharedPtr thisSh =
        std::dynamic_pointer_cast<CanHaveAttributes>(shared_from_this());
    return CanHaveAttributes::AttrIterator(thisSh, GetNumAttr());
}

CanHaveAttributes::AttrIterator::AttrIterator(CanHaveAttributesSharedPtr obj,
                                              hsize_t idx)
    : m_obj(obj), m_idx(-1), m_next(idx), m_size(obj->GetNumAttr())
{
    ++*this;
}

const std::string &CanHaveAttributes::AttrIterator::operator*()
{
    return m_currentName;
}
CanHaveAttributes::AttrIterator &CanHaveAttributes::AttrIterator::operator++()
{
    m_idx = m_next;
    if (m_next < m_size)
    {
        H5_CALL(H5Aiterate2, (m_obj->GetId(), H5_INDEX_CRT_ORDER, H5_ITER_INC,
                              &m_next, AttrIterator::helper, this));
    }
    return *this;
}
bool CanHaveAttributes::AttrIterator::operator==(
    const CanHaveAttributes::AttrIterator &other) const
{
    return m_obj == other.m_obj && m_idx == other.m_idx;
}

herr_t CanHaveAttributes::AttrIterator::helper(hid_t g_id, const char *name,
                                               const H5A_info_t * info,
                                               void *op_data)
{
    boost::ignore_unused(g_id, info);
    CanHaveAttributes::AttrIterator *iter =
        static_cast<CanHaveAttributes::AttrIterator *>(op_data);
    iter->m_currentName = name;
    return 1;
}

DataSpaceSharedPtr DataSpace::Null()
{
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_NULL));
    return ans;
}

DataSpaceSharedPtr DataSpace::Scalar()
{
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_SCALAR));
    return ans;
}
DataSpaceSharedPtr DataSpace::OneD(hsize_t size)
{
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate_simple, (1, &size, NULL));
    return ans;
}

DataSpace::DataSpace() : Object()
{
}

DataSpace::DataSpace(hid_t id) : Object(id)
{
}

DataSpace::DataSpace(const std::vector<hsize_t> &dims) : Object()
{
    int rank = dims.size();
    H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], NULL));
}

DataSpace::DataSpace(const hsize_t size, const hsize_t max) : Object()
{
    const hsize_t *max_p = NULL;
    if (max != (H5S_UNLIMITED - 1))
        max_p = &max;
    H5_CONSTRUCT(m_Id, H5Screate_simple, (1, &size, max_p));
}

DataSpace::DataSpace(const std::vector<hsize_t> &dims,
                     const std::vector<hsize_t> &max_dims)
    : Object()
{
    int rank = dims.size();
    H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], &max_dims[0]));
}

DataSpace::~DataSpace()
{
    Close();
}

void DataSpace::Close()
{
    H5_CALL(H5Sclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}
void DataSpace::SelectRange(const hsize_t start, const hsize_t count)
{
    H5_CALL(H5Sselect_hyperslab,
            (m_Id, H5S_SELECT_SET, &start, NULL, &count, NULL));
}
void DataSpace::AppendRange(const hsize_t start, const hsize_t count)
{
    H5_CALL(H5Sselect_hyperslab,
            (m_Id, H5S_SELECT_OR, &start, NULL, &count, NULL));
}
void DataSpace::SelectRange(const std::vector<hsize_t> start,
                            const std::vector<hsize_t> count)
{
    H5_CALL(H5Sselect_hyperslab,
            (m_Id, H5S_SELECT_SET, &start[0], NULL, &count[0], NULL));
}
void DataSpace::AppendRange(const std::vector<hsize_t> start,
                            const std::vector<hsize_t> count)
{
    H5_CALL(H5Sselect_hyperslab,
            (m_Id, H5S_SELECT_OR, &start[0], NULL, &count[0], NULL));
}
void DataSpace::SetSelection(const hsize_t num_elmt,
                             const std::vector<hsize_t> &coords)
{
    if (num_elmt == 0)
    {
        H5_CALL(H5Sselect_none, (m_Id));
    }
    else
    {
        H5_CALL(H5Sselect_elements,
                (m_Id, H5S_SELECT_SET, num_elmt, &coords[0]));
    }
}

void DataSpace::ClearRange()
{
    H5_CALL(H5Sselect_none, (m_Id));
}

hsize_t DataSpace::GetSize()
{
    return H5Sget_simple_extent_npoints(m_Id);
}

std::vector<hsize_t> DataSpace::GetDims()
{
    int ndims = H5Sget_simple_extent_ndims(m_Id);
    std::vector<hsize_t> ret(ndims, 0);
    H5Sget_simple_extent_dims(m_Id, &ret[0], NULL);
    return ret;
}

DataType::DataType(hid_t id) : Object(id)
{
}
DataTypeSharedPtr DataType::String(size_t len)
{
    DataTypeSharedPtr s1  = PredefinedDataType::CS1();
    DataTypeSharedPtr ans = s1->Copy();
    if (len == 0)
    {
        len = H5T_VARIABLE;
    }
    H5_CALL(H5Tset_size, (ans->GetId(), len));
    return ans;
}

void CompoundDataType::Close()
{
    H5_CALL(H5Tclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

void DataType::Close()
{
    H5_CALL(H5Tclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

DataTypeSharedPtr DataType::Copy() const
{
    hid_t ans_id = H5I_INVALID_HID;
    H5_CONSTRUCT(ans_id, H5Tcopy, (m_Id));
    DataTypeSharedPtr ans(new DataType(ans_id));
    return ans;
}

DataTypeSharedPtr PredefinedDataType::CS1()
{
    return DataTypeSharedPtr(new PredefinedDataType(H5T_C_S1));
}

CompoundDataTypeSharedPtr CompoundDataType::Create(size_t sz)
{
    return std::shared_ptr<CompoundDataType>(
        new CompoundDataType(H5Tcreate(H5T_COMPOUND, sz)));
}

CompoundDataType::CompoundDataType(hid_t id) : DataType(id)
{
}

PredefinedDataType::PredefinedDataType(hid_t id) : DataType(id)
{
}

void PredefinedDataType::Close()
{
    // No-op
    m_Id = H5I_INVALID_HID;
}

template <> const hid_t DataTypeTraits<char>::NativeType = H5T_NATIVE_CHAR;

template <> const hid_t DataTypeTraits<int>::NativeType = H5T_NATIVE_INT;

template <>
const hid_t DataTypeTraits<unsigned int>::NativeType = H5T_NATIVE_UINT;

template <>
const hid_t DataTypeTraits<unsigned long>::NativeType = H5T_NATIVE_ULONG;

template <>
const hid_t DataTypeTraits<unsigned long long>::NativeType = H5T_NATIVE_ULLONG;

template <> const hid_t DataTypeTraits<double>::NativeType = H5T_NATIVE_DOUBLE;

template <> const hid_t DataTypeTraits<BasisType>::NativeType = H5T_NATIVE_INT;


AttributeSharedPtr Attribute::Create(hid_t parent, const std::string &name,
                                     DataTypeSharedPtr type,
                                     DataSpaceSharedPtr space)
{
    hid_t id;
    H5_CONSTRUCT(id, H5Acreate2, (parent, name.c_str(), type->GetId(),
                                  space->GetId(), H5P_DEFAULT, H5P_DEFAULT));
    return AttributeSharedPtr(new Attribute(id));
}

AttributeSharedPtr Attribute::Open(hid_t parent, const std::string &name)
{
    hid_t id;
    H5_CONSTRUCT(id, H5Aopen, (parent, name.c_str(), H5P_DEFAULT));
    return AttributeSharedPtr(new Attribute(id));
}

Attribute::~Attribute()
{
    Close();
}
void Attribute::Close()
{
    H5_CALL(H5Aclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

DataSpaceSharedPtr Attribute::GetSpace() const
{
    return DataSpaceSharedPtr(new DataSpace(H5Aget_space(m_Id)));
}

File::File(hid_t id) : Object(id)
{
}
FileSharedPtr File::Create(const std::string &filename, unsigned mode,
                           PListSharedPtr createPL, PListSharedPtr accessPL)
{
    hid_t id;
    H5_CONSTRUCT(id, H5Fcreate, (filename.c_str(), mode, createPL->GetId(),
                                 accessPL->GetId()));
    return FileSharedPtr(new File(id));
}
FileSharedPtr File::Open(const std::string &filename, unsigned mode,
                         PListSharedPtr accessPL)
{
    hid_t id;
    H5_CONSTRUCT(id, H5Fopen, (filename.c_str(), mode, accessPL->GetId()));
    return FileSharedPtr(new File(id));
}

File::~File()
{
    if (m_Id != H5I_INVALID_HID)
    {
        Close();
    }
}

void File::Close()
{
    H5_CALL(H5Fclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}
hsize_t File::GetNumElements()
{
    GroupSharedPtr root = OpenGroup("/");
    return root->GetNumElements();
}
Group::Group(hid_t id) : Object(id)
{
}

Group::~Group()
{
    if (m_Id != H5I_INVALID_HID)
    {
        Close();
    }
}

void Group::Close()
{
    H5_CALL(H5Gclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

hsize_t Group::GetNumElements()
{
    H5G_info_t info;
    H5_CALL(H5Gget_info, (m_Id, &info));
    return info.nlinks;
}

std::vector<std::string> Group::GetElementNames()
{
    std::vector<std::string> ret;
    for(int i = 0; i < GetNumElements(); i++)
    {
        char name[50];
        H5Gget_objname_by_idx(m_Id, (size_t)i, name, 50 );
        ret.push_back(std::string(name));
    }
    return ret;
}

DataSet::DataSet(hid_t id) : Object(id)
{
}

DataSet::~DataSet()
{
    Close();
}

void DataSet::Close()
{
    H5_CALL(H5Dclose, (m_Id));
    m_Id = H5I_INVALID_HID;
}

DataSpaceSharedPtr DataSet::GetSpace() const
{
    return DataSpaceSharedPtr(new DataSpace(H5Dget_space(m_Id)));
}
}
}
}
