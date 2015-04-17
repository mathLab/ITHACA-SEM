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

            // Forward declare
            class Object;
            typedef boost::shared_ptr<Object> ObjectSharedPtr;
            class DataType;
            typedef boost::shared_ptr<DataType> DataTypeSharedPtr;
            class DataSpace;
            typedef boost::shared_ptr<DataSpace> DataSpaceSharedPtr;
            class Attribute;
            typedef boost::shared_ptr<Attribute> AttributeSharedPtr;
            class Group;
            typedef boost::shared_ptr<Group> GroupSharedPtr;
            class File;
            typedef boost::shared_ptr<File> FileSharedPtr;

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

            /// Mixin for objects that contain groups and datasets (Group and File)
            class CanHaveGroups : public virtual Object
            {
                public:
                    GroupSharedPtr CreateGroup(const std::string& name);
            };

            /// Mixin for objects that can have attributes (Group, DataSet, DataType)
            class CanHaveAttributes : public virtual Object
            {
                public:
                    AttributeSharedPtr CreateAttribute(const std::string& name,
                            DataTypeSharedPtr type, DataSpaceSharedPtr space);

                    template<class T>
                    void SetAttribute(const std::string& name, const T& value);
                    template<class T>
                    void SetAttribute(const std::string& name,
                            const std::vector<T>& value);
            };

            /// HDF5 DataSpace wrapper
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

            /// Traits class for HDF5 data types. Specialise this for
            template<class T>
            struct DataTypeTraits
            {
                    /**
                     * Get the address of the start of the data.
                     * Default implementation just uses "&"
                     */
                    static const void* GetAddress(const T& obj);
                    /**
                     * Return a DataType object representing T.
                     * Default implementation just calls PredefinedDataType::Native<T>()
                     */
                    static DataTypeSharedPtr GetType();
                    /***
                     * Define this for a specialision for any HDF5 NATIVE type you want to use.
                     * See http://hdfgroup.org/HDF5/doc/UG/UG_frame11Datatypes.html
                     */
                    static hid_t NativeType;
            };

            /// Wrap and HDF5 data type object
            class DataType : public virtual Object, public CanHaveAttributes
            {
                public:
                    static DataTypeSharedPtr String(size_t len = 0);
                    virtual void Close();
                    DataTypeSharedPtr Copy() const;
                protected:
                    DataType(hid_t id);
            };

            /// Predefined HDF data types that must not be closed when done with.
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

            /// HDF5 Attribute Wrapper
            class Attribute : public Object
            {
                public:
                    ~Attribute();
                    void Close();
                private:
                    Attribute(hid_t parent, const std::string& name,
                            DataTypeSharedPtr type, DataSpaceSharedPtr space);
                    friend class CanHaveAttributes;
            };

            /// HDF5 file wrapper
            class File : public CanHaveGroups
            {
                public:
                    File(const std::string& filename, unsigned mode);
                    ~File();
                    void Close();

            };

            /// HDF5 Group wrapper
            class Group : public CanHaveAttributes, public CanHaveGroups
            {
                public:
                    ~Group();
                    void Close();
                private:
                    Group(hid_t id);
                    friend class CanHaveGroups;
            };

//            class DataSet : public CanHaveAttributes
//            {
//                private:
//                    DataSet(boost::shared_ptr<CanHaveAttributes> parent,
//                            const std::string& name, DataTypeSharedPtr type,
//                            DataSpaceSharedPtr space);
//                    friend class CanHaveAttributes;
//            };


        }
    }
}

#include <LibUtilities/BasicUtils/H5.hpp>
#endif
