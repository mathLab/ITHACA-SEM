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
//  Description: Minimal HDF5 wrapper
//
////////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/H5.h>
#include <boost/make_shared.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {

            Object::Object() :
                    m_Id(H5I_INVALID_HID)
            {

            }
            Object::Object(hid_t id) :
                    m_Id(id)
            {

            }
            Object::~Object()
            {
            }

            PList::PList() :
                    Object(H5P_DEFAULT)
            {
            }
            PList::PList(hid_t cls) :
                    Object()
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
            void PList::SetChunk(const std::vector<hsize_t>& dims)
            {
                H5_CALL(H5Pset_chunk, (m_Id, dims.size(), &dims[0]));
            }
            void PList::SetDeflate(const unsigned level)
            {
                H5_CALL(H5Pset_deflate, (m_Id, level))
            }

            GroupSharedPtr CanHaveGroupsDataSets::CreateGroup(
                    const std::string& name, PListSharedPtr createPL,
                    PListSharedPtr accessPL)
            {
                hid_t grp;
                H5_CONSTRUCT(grp, H5Gcreate,
                        (m_Id, name.c_str(), H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));
                GroupSharedPtr ans(new Group(grp));
                return ans;
            }

            DataSetSharedPtr CanHaveGroupsDataSets::CreateDataSet(
                    const std::string& name, DataTypeSharedPtr type,
                    DataSpaceSharedPtr space, PListSharedPtr createPL,
                    PListSharedPtr accessPL)
            {
                DataSetSharedPtr ans(
                        new DataSet(GetId(), name, type, space, createPL,
                                accessPL));
                return ans;
            }

            AttributeSharedPtr CanHaveAttributes::CreateAttribute(
                    const std::string& name, DataTypeSharedPtr type,
                    DataSpaceSharedPtr space)
            {
                AttributeSharedPtr ans(
                        new Attribute(GetId(), name, type, space));
                return ans;
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

            DataSpace::DataSpace() :
                    Object()
            {
            }

            DataSpace::DataSpace(const std::vector<hsize_t>& dims) :
                    Object()
            {
                int rank = dims.size();
                H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], NULL));
            }

            DataSpace::DataSpace(const hsize_t size, const hsize_t max) :
                    Object()
            {
                const hsize_t* max_p = NULL;
                if (max != (H5S_UNLIMITED - 1))
                    max_p = &max;
                H5_CONSTRUCT(m_Id, H5Screate_simple, (1, &size, max_p));
            }

            DataSpace::DataSpace(const std::vector<hsize_t>& dims,
                    const std::vector<hsize_t>& max_dims) :
                    Object()
            {
                int rank = dims.size();
                H5_CONSTRUCT(m_Id, H5Screate_simple,
                        (rank, &dims[0], &max_dims[0]));
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

            DataType::DataType(hid_t id) :
                    Object(id)
            {
            }
            DataTypeSharedPtr DataType::String(size_t len)
            {
                DataTypeSharedPtr s1 = PredefinedDataType::CS1();
                DataTypeSharedPtr ans = s1->Copy();
                if (len == 0)
                    len = H5T_VARIABLE;
                H5_CALL(H5Tset_size, (ans->GetId(), len));
                return ans;
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

            PredefinedDataType::PredefinedDataType(hid_t id) :
                    DataType(id)
            {
            }

            void PredefinedDataType::Close()
            {
                // No-op
                m_Id = H5I_INVALID_HID;
            }

            template<>
            const hid_t DataTypeTraits<char>::NativeType = H5T_NATIVE_CHAR;

            template<>
            const hid_t DataTypeTraits<int>::NativeType = H5T_NATIVE_INT;

            template<>
            const hid_t DataTypeTraits<unsigned int>::NativeType =
                    H5T_NATIVE_UINT;

            template<>
            const hid_t DataTypeTraits<double>::NativeType = H5T_NATIVE_DOUBLE;

            Attribute::Attribute(hid_t parent, const std::string& name,
                    DataTypeSharedPtr type, DataSpaceSharedPtr space)
            {
                H5_CONSTRUCT(m_Id, H5Acreate,
                        (parent, name.c_str(), type->GetId(), space->GetId(), H5P_DEFAULT, H5P_DEFAULT));
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

            File::File(const std::string& filename, unsigned mode,
                    PListSharedPtr createPL, PListSharedPtr accessPL)
            {
                H5_CONSTRUCT(m_Id, H5Fcreate,
                        (filename.c_str(), mode, createPL->GetId(), accessPL->GetId()));
            }

            File::~File()
            {
                Close();
            }

            void File::Close()
            {
                H5_CALL(H5Fclose, (m_Id));
                m_Id = H5I_INVALID_HID;
            }

            Group::Group(hid_t id) :
                    Object(id)
            {
            }

            Group::~Group()
            {
                Close();
            }

            void Group::Close()
            {
                H5_CALL(H5Gclose, (m_Id));
                m_Id = H5I_INVALID_HID;
            }

            DataSet::DataSet(hid_t parent, const std::string& name,
                    DataTypeSharedPtr type, DataSpaceSharedPtr space,
                    PListSharedPtr createPL, PListSharedPtr accessPL)
            {
                H5_CONSTRUCT(m_Id, H5Dcreate,
                        (parent, name.c_str(), type->GetId(), space->GetId(), H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));
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

        }
    }
}
