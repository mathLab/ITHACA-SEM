///////////////////////////////////////////////////////////////////////////////
//
// File H5.h
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
// Description: Simple OO wrapper around HDF5
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_H5_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_H5_H

#include <string>
#include <exception>
#include <vector>
#include <hdf5.h>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {

            class Error : public std::exception
            {

            };
            /// HDF5 base class
            class Object : public boost::enable_shared_from_this<Object>
            {
                public:
                    virtual void Close() = 0;
                    inline hid_t GetId() const
                    {
                        return m_Id;
                    }
                    // Overload cast to the HDF5 ID type to make objects
                    // transparently usable in the HDF5 API.
                    inline operator hid_t() const
                    {
                        return GetId();
                    }
                protected:
                    Object();
                    Object(hid_t id);
                    virtual ~Object();
                    hid_t m_Id;
            };
            typedef boost::shared_ptr<Object> ObjectSharedPtr;

            // Fwd declare
            class DataType;
            typedef boost::shared_ptr<DataType> DataTypeSharedPtr;
            class DataSpace;
            typedef boost::shared_ptr<DataSpace> DataSpaceSharedPtr;
            class Attribute;
            typedef boost::shared_ptr<Attribute> AttributeSharedPtr;
            class Group;
            typedef boost::shared_ptr<Group> GroupSharedPtr;

            class Location : public virtual Object
            {
            };
            typedef boost::shared_ptr<Location> LocationSharedPtr;

            /// Common functionality of Groups and Files
            class CanHaveGroups : public virtual Object
            {
                public:
                    GroupSharedPtr CreateGroup(const std::string& name);
            };

            /// HDF5 file wrapper
            class File : public Location, public CanHaveGroups
            {
                public:
                    File(const std::string& filename, unsigned mode);
                    ~File();
                    void Close();

            };
            typedef boost::shared_ptr<File> FileSharedPtr;

            /// Base for Groups, DataSets and DataTypes
            class NamedObject : public Location
            {
                public:
                    AttributeSharedPtr CreateAttribute(const std::string& name,
                            DataTypeSharedPtr type, DataSpaceSharedPtr space);

                    template<class T>
                    void SetAttribute(const std::string& name, const T& value);
                    template<class T>
                    void SetAttribute(const std::string& name, const std::vector<T>& value);

                    //void SetAttribute(const std::string& name, const std::string& value);

//                    void SetStringAttribute(const std::string& name,
//                            const std::string& value);
            };
            typedef boost::shared_ptr<NamedObject> NamedObjectSharedPtr;

            /// HDF5 Group wrapper
            class Group : public NamedObject, public CanHaveGroups
            {
                public:
                    ~Group();
                    void Close();
                private:
                    Group(hid_t id);
                    friend class CanHaveGroups;
            };

            class Attribute : public Object
            {
                public:
                    ~Attribute();
                    void Close();
//                    Read();
//                    void Write(const std::string& data);
//                    template<class T>
//                    void Write(const std::vector<T>& data);
                private:
                    Attribute(NamedObjectSharedPtr parent,
                            const std::string& name, DataTypeSharedPtr type,
                            DataSpaceSharedPtr space);
                    friend class NamedObject;
            };
            class DataSet : public NamedObject
            {
                private:
                    DataSet(NamedObjectSharedPtr parent,
                            const std::string& name, DataTypeSharedPtr type,
                            DataSpaceSharedPtr space);
                    friend class NamedObject;
            };
            class DataSpace : public Object
            {
                public:
                    static DataSpaceSharedPtr Null();
                    static DataSpaceSharedPtr Scalar();
                    static DataSpaceSharedPtr OneD(hsize_t size);

                    DataSpace();
                    DataSpace(hsize_t size, hsize_t max = H5S_UNLIMITED - 1);
                    DataSpace(const std::vector<hsize_t>& dims);
                    DataSpace(const std::vector<hsize_t>& dims,
                            const std::vector<hsize_t>& max_dims);
                    ~DataSpace();

                    void Close();
            };

            template<class T>
            struct DataTypeTraits
            {
                    static const void* GetAddress(const T& obj);
                    static DataTypeSharedPtr GetType();
                    static hid_t NativeType;
                    static const bool IsPredefined;
            };

            class DataType : public NamedObject
            {
                public:
                    static DataTypeSharedPtr String(size_t len = 0);
                    virtual void Close();
                    DataTypeSharedPtr Copy() const;
                protected:
                    DataType(hid_t id);
            };

            class PredefinedDataType : public DataType
            {
                public:
                    template<class T>
                    static inline DataTypeSharedPtr Native()
                    {
                        return DataTypeSharedPtr(
                                new PredefinedDataType(
                                        DataTypeTraits<T>::NativeType));
                    }
                    static DataTypeSharedPtr CS1();
                    void Close();
                private:
                    PredefinedDataType (hid_t);

            };
        }
    }
}

#include <LibUtilities/BasicUtils/H5.hpp>
#endif
