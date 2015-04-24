namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {
            template<class T>
            typename DataTypeConversionPolicy<T>::ConvertedType DataTypeConversionPolicy<
                    T>::Convert(const T& obj)
            {
                return obj;
            }

            template<class T>
            DataTypeSharedPtr DataTypeTraits<T>::GetType()
            {
                return PredefinedDataType::Native<T>();
            }
            template<class T>
            const void* DataTypeTraits<T>::GetAddress(
                    DataTypeTraits<T>::ConvertedType& obj)
            {
                return &obj;
            }
            template<class T>
            typename DataTypeTraits<T>::ConvertedType DataTypeTraits<T>::Convert(
                    const T& obj)
            {
                return Converter::Convert(obj);
            }

            template<>
            struct DataTypeConversionPolicy<std::string>
            {
                    static const bool MustConvert = true;
                    typedef const char* ConvertedType;
                    typedef const char* ConvertedVectorElemType;
                    inline static ConvertedType Convert(const std::string & obj)
                    {
                        return obj.c_str();
                    }

            };

            template<>
            inline DataTypeSharedPtr DataTypeTraits<std::string>::GetType()
            {
                return DataType::String();
            }

            template<class T>
            DataTypeSharedPtr PredefinedDataType::Native()
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
                typedef std::vector<
                        typename DataTypeConversionPolicy<T>::ConvertedVectorElemType> Vec;
                Vec converted_vals;
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::OneD(value.size());
                AttributeSharedPtr attr = CreateAttribute(name, type, space);

                const void* converted_buf = NULL;
                if (DataTypeConversionPolicy<T>::MustConvert)
                {
                    converted_vals.resize(value.size());
                    for (size_t i = 0; i < value.size(); ++i)
                        converted_vals[i] =
                                DataTypeConversionPolicy<T>::Convert(value[i]);
                    converted_buf = &converted_vals[0];
                }
                else
                {
                    converted_buf = &value[0];
                }

                H5_CALL(H5Awrite,
                        (attr->GetId(), type->GetId(), converted_buf));
            }

            template<class T>
            DataSetSharedPtr CanHaveGroupsDataSets::CreateWriteDataSet(
                    const std::string& name, const std::vector<T>& data, PListSharedPtr createPL, PListSharedPtr accessPL)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::OneD(data.size());
                DataSetSharedPtr dataset = CreateDataSet(name, type, space, createPL, accessPL);
                dataset->Write(data);
                return dataset;
            }

        }
    }
}
