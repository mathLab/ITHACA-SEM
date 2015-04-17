namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {
            template<typename T>
            void CanHaveAttributes::SetAttribute(const std::string& name,
                    const T& value)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::Scalar();
                AttributeSharedPtr attr = CreateAttribute(name, type, space);
                H5Awrite(attr->GetId(), type->GetId(),
                        DataTypeTraits<T>::GetAddress(value));
            }

            template<typename T>
            void CanHaveAttributes::SetAttribute(const std::string& name,
                    const std::vector<T>& value)
            {
                DataTypeSharedPtr type = DataTypeTraits<T>::GetType();
                DataSpaceSharedPtr space = DataSpace::OneD(value.size());
                AttributeSharedPtr attr = CreateAttribute(name, type, space);
                H5Awrite(*attr, *type, &value[0]);

            }

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
            template<>
            inline const void* DataTypeTraits<std::string>::GetAddress(
                    const std::string& obj)
            {
                return &obj[0];
            }

            template<>
            inline DataTypeSharedPtr DataTypeTraits<std::string>::GetType()
            {
                return DataType::String();
            }
        }
    }
}
