namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {
            template<class T>
            DataTypeSharedPtr DataTypeTraits<T>::GetType()
            {
                return PredefinedDataType::Native<T>();
            }

            template<class T>
            const void* DataTypeTraits<T>::GetAddress(const T& obj)
            {
                return &obj;
            }
            template<class T>
            typename DataTypeTraits<T>::ConvertedType DataTypeTraits<T>::Convert(
                    const T& obj)
            {
                return obj;
            }

            template<>
            struct DataTypeTraits<std::string>
            {
                    typedef const char* ConvertedType;
                    static inline const void* GetAddress(ConvertedType& obj)
                    {
                        return &obj;
                    }
                    static inline DataTypeSharedPtr GetType()
                    {
                        return DataType::String();
                    }

                    static ConvertedType Convert(const std::string& obj)
                    {
                        return obj.c_str();
                    }
            };

            template<class T>
            inline DataTypeSharedPtr PredefinedDataType::Native()
            {
                return DataTypeSharedPtr(
                        new PredefinedDataType(DataTypeTraits<T>::NativeType));
            }

            template<typename T>
            void CanHaveAttributes::SetAttribute(const std::string& name,
                    const T& value)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::Scalar();
                AttributeSharedPtr attr = CreateAttribute(name, type, space);

                typename DataTypeTraits<T>::ConvertedType conv = DataTypeTraits<
                        T>::Convert(value);
                H5_CALL(H5Awrite,
                        (attr->GetId(), type->GetId(), DataTypeTraits<T>::GetAddress(
                                conv)));
            }

            template<typename T>
            void CanHaveAttributes::SetAttribute(const std::string& name,
                    const std::vector<T>& value)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::OneD(value.size());
                AttributeSharedPtr attr = CreateAttribute(name, type, space);
                H5_CALL(H5Awrite, (attr->GetId(), type->GetId(), &value[0]));

            }

            template<class T>
            DataSetSharedPtr CanHaveGroupsDataSets::CreateWriteDataSet(
                    const std::string& name, const std::vector<T>& data)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::OneD(data.size());
                DataSetSharedPtr dataset = CreateDataSet(name, type, space);
                dataset->Write(data);
                return dataset;
            }

        }
    }
}
