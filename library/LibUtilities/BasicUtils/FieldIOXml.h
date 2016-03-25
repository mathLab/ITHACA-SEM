///////////////////////////////////////////////////////////////////////////////
//
// File FieldIOXml.h
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
// Description: Field IO to/from XML
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOXML_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIOXML_H

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
    namespace LibUtilities
    {
        class XmlDataSource : public DataSource
        {
            public:
                XmlDataSource(TiXmlDocument& doc);

                XmlDataSource(const std::string& fn);
                ~XmlDataSource();
                TiXmlDocument& Get();
                const TiXmlDocument& Get() const;

                static DataSourceSharedPtr create(const std::string& fn);
                static DataSourceSharedPtr create(TiXmlDocument& fn);
            private:
                TiXmlDocument* m_doc;
        };

        /// Class for operating on XML FLD files
        class FieldIOXml : public FieldIO
        {
            public:

                /// Creates an instance of this class
                LIB_UTILITIES_EXPORT
                static FieldIOSharedPtr create(
                    LibUtilities::CommSharedPtr pComm,
                    bool sharedFilesystem)
                {
                    return MemoryManager<FieldIOXml>::AllocateSharedPtr(
                        pComm, sharedFilesystem);
                }

                /// Name of class
                LIB_UTILITIES_EXPORT
                static std::string className;

                FieldIOXml(LibUtilities::CommSharedPtr pComm,
                           bool sharedFilesystem);

                LIB_UTILITIES_EXPORT
                std::string SetUpOutput(const std::string outname);

                /// Imports the definition of the fields.
                LIB_UTILITIES_EXPORT
                virtual void ImportFieldDefs(DataSourceSharedPtr dataSource,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        bool expChild);

                /// Imports the data fields.
                LIB_UTILITIES_EXPORT
                void ImportFieldData(DataSourceSharedPtr dataSource,
                        const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata);

                LIB_UTILITIES_EXPORT
                void WriteMultiFldFileIDs(const std::string &outfile,
                        const std::vector<std::string> fileNames,
                        std::vector<std::vector<unsigned int> > &elementList,
                        const FieldMetaDataMap &fieldinfomap =
                                NullFieldMetaDataMap);

                LIB_UTILITIES_EXPORT void SetUpFieldMetaData(
                    const std::string                             outname,
                    const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    const FieldMetaDataMap                       &fieldmetadatamap);

                LIB_UTILITIES_EXPORT
                void ImportMultiFldFileIDs(const std::string &inFile,
                        std::vector<std::string> &fileNames,
                        std::vector<std::vector<unsigned int> > &elementList,
                        FieldMetaDataMap &fieldmetadatamap);

                /// Imports an FLD file.
                LIB_UTILITIES_EXPORT
                void v_Import(const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata =
                                NullVectorNekDoubleVector,
                        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
                        const Array<OneD, int> ElementiDs = NullInt1DArray);

                inline virtual const std::string& GetClassName() const
                {
                    return className;
                }

            private:
                /// Write data in FLD format
                LIB_UTILITIES_EXPORT
                virtual void v_Write(const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata,
                        const FieldMetaDataMap &fieldinfomap =
                                NullFieldMetaDataMap);

                /// Imports the definition of the meta data
                LIB_UTILITIES_EXPORT
                DataSourceSharedPtr v_ImportFieldMetaData(
                        std::string filename,
                        FieldMetaDataMap &fieldmetadatamap);
        };
    }
}
#endif
