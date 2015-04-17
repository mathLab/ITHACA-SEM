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

            GroupSharedPtr CanHaveGroups::CreateGroup(const std::string& name)
            {
                hid_t grp = H5Gcreate(m_Id, name.c_str(), H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
                GroupSharedPtr ans(new Group(grp));
                return ans;
            }

            AttributeSharedPtr CanHaveAttributes::CreateAttribute(
                    const std::string& name, DataTypeSharedPtr type,
                    DataSpaceSharedPtr space)
            {
                AttributeSharedPtr ans(new Attribute(GetId(), name, type, space));
                return ans;
            }

            DataSpaceSharedPtr DataSpace::Null()
            {
                DataSpaceSharedPtr ans(new DataSpace);
                ans->m_Id = H5Screate(H5S_NULL);
                return ans;
            }

            DataSpaceSharedPtr DataSpace::Scalar()
            {
                DataSpaceSharedPtr ans(new DataSpace);
                ans->m_Id = H5Screate(H5S_SCALAR);
                return ans;
            }
            DataSpaceSharedPtr DataSpace::OneD(hsize_t size)
            {
                DataSpaceSharedPtr ans(new DataSpace);
                ans->m_Id = H5Screate_simple(1, &size, NULL);
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
                m_Id = H5Screate_simple(rank, &dims[0], NULL);
            }

            DataSpace::DataSpace(const hsize_t size, const hsize_t max) :
                    Object()
            {
                const hsize_t* max_p = NULL;
                if (max != (H5S_UNLIMITED - 1))
                    max_p = &max;
                m_Id = H5Screate_simple(1, &size, max_p);
            }

            DataSpace::DataSpace(const std::vector<hsize_t>& dims,
                    const std::vector<hsize_t>& max_dims) :
                    Object()
            {
                int rank = dims.size();
                m_Id = H5Screate_simple(rank, &dims[0], &max_dims[0]);
            }

            DataSpace::~DataSpace()
            {
                Close();
            }

            void DataSpace::Close()
            {
                H5Sclose(m_Id);
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
                H5Tset_size(ans->GetId(), len);
                return ans;
            }

            void DataType::Close()
            {
                H5Tclose(m_Id);
                m_Id = H5I_INVALID_HID;
            }

            DataTypeSharedPtr DataType::Copy() const
            {
                hid_t ans_id = H5Tcopy(H5T_C_S1);
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
            hid_t DataTypeTraits<char>::NativeType = H5T_NATIVE_CHAR;

            template<>
            hid_t DataTypeTraits<int>::NativeType = H5T_NATIVE_INT;

            template<>
            hid_t DataTypeTraits<unsigned int>::NativeType = H5T_NATIVE_UINT;

            template<>
            hid_t DataTypeTraits<double>::NativeType = H5T_NATIVE_DOUBLE;

            Attribute::Attribute(hid_t parent,
                    const std::string& name, DataTypeSharedPtr type,
                    DataSpaceSharedPtr space)
            {
                m_Id = H5Acreate(parent, name.c_str(), type->GetId(),
                        space->GetId(), H5P_DEFAULT, H5P_DEFAULT);
            }

            Attribute::~Attribute()
            {
                Close();
            }
            void Attribute::Close()
            {
                H5Aclose(m_Id);
            }


            File::File(const std::string& filename, unsigned mode)
            {
                m_Id = H5Fcreate(filename.c_str(), mode, H5P_DEFAULT,
                        H5P_DEFAULT);
                if (m_Id < 0)
                    throw Error();
            }

            File::~File()
            {
                Close();
            }

            void File::Close()
            {
                H5Fclose(m_Id);
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
                H5Gclose(m_Id);
                m_Id = H5I_INVALID_HID;
            }

        }
    }
}
