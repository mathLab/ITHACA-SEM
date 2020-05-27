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

#include <boost/algorithm/string/predicate.hpp>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace LibUtilities
{

/**
 * @class Class encapsulating simple XML data source using TinyXML.
 */
class XmlDataSource : public DataSource
{
public:
    /// Default constructor.
    XmlDataSource(TiXmlDocument &doc) : m_doc(&doc), m_needsFree(false) { }

    /// Constructor based on filename.
    XmlDataSource(const std::string &fn) : m_needsFree(true)
    {
        m_doc = new TiXmlDocument(fn);
        bool loadOkay = m_doc->LoadFile();
        std::stringstream errstr;
        errstr << "Unable to load file: " << fn << std::endl;
        errstr << "Reason: " << m_doc->ErrorDesc() << std::endl;
        errstr << "Position: Line " << m_doc->ErrorRow() << ", Column "
               << m_doc->ErrorCol() << std::endl;
        ASSERTL0(loadOkay, errstr.str());
    }

    /// Destructor cleans up memory usage.
    ~XmlDataSource()
    {
        if (m_needsFree)
        {
            delete m_doc;
        }
    }

    /// Return the TinyXML document of this source.
    TiXmlDocument &Get()
    {
        return *m_doc;
    }

    /// Return the TinyXML document of this source.
    const TiXmlDocument &Get() const
    {
        return *m_doc;
    }

    /// Create a new XML data source based on the filename.
    static DataSourceSharedPtr create(const std::string &fn)
    {
        return DataSourceSharedPtr(new XmlDataSource(fn));
    }

    /// Create a new XML data source based on a TiXmlDocument.
    static DataSourceSharedPtr create(TiXmlDocument &fn)
    {
        return DataSourceSharedPtr(new XmlDataSource(fn));
    }

private:
    /// Internal TinyXML document storage.
    TiXmlDocument *m_doc;
    /// Boolean dictating whether document needs to be freed or not.
    bool m_needsFree;
};
typedef std::shared_ptr<XmlDataSource> XmlDataSourceSharedPtr;

/**
 * @class Simple class for writing XML hierarchical data using TinyXML.
 */
class XmlTagWriter : public TagWriter
{
public:
    /// Default constructor.
    XmlTagWriter(TiXmlElement *elem) : m_El(elem) {}

    /// Add a child node.
    virtual TagWriterSharedPtr AddChild(const std::string &name)
    {
        TiXmlElement *child = new TiXmlElement(name.c_str());
        m_El->LinkEndChild(child);
        return TagWriterSharedPtr(new XmlTagWriter(child));
    }

    /// Set an attribute key/value pair on this tag.
    virtual void SetAttr(const std::string &key, const std::string &val)
    {
        if (boost::starts_with(key, "XML_"))
        {
            // Auto-expand XML parameters.
            std::string elmtName = key.substr(4);
            TiXmlElement *child = new TiXmlElement(elmtName.c_str());

            // Parse string we're given
            TiXmlDocument doc;
            doc.Parse(val.c_str());

            TiXmlElement *e = doc.FirstChildElement();

            while (e)
            {
                child->LinkEndChild(e->Clone());
                e = e->NextSiblingElement();
            }

            m_El->LinkEndChild(child);
        }
        else
        {
            TiXmlElement *child = new TiXmlElement(key.c_str());
            child->LinkEndChild(new TiXmlText(val.c_str()));
            m_El->LinkEndChild(child);
        }
    }

private:
    /// Internal TinyXML document storage.
    TiXmlElement *m_El;
};
typedef std::shared_ptr<XmlTagWriter> XmlTagWriterSharedPtr;

/**
 * @class Class for operating on XML FLD files.
 *
 * This class is the default for Nektar++ output. It reads and writes one XML
 * file per processor that represents the underlying field data. For serial
 * output, the format of an XML file obeys the following structure:
 *
 * ```
 * <NEKTAR>
 *     <Metadata>
 *         ...
 *     </Metadata>
 *     <ELEMENT FIELDS="..." ...> data1 </ELEMENT>
 *     <ELEMENT FIELDS="..." ...> data2 </ELEMENT>
 *     ...
 * </NEKTAR>
 * ```
 *
 *   - Metadata is converted as key/value pairs in the `<Metadata>` tag.
 *   - There are one or more ELEMENT blocks, whose attributes correspond with
 *     the FieldDefinitions class variables.
 *   - Element data is stored as a base64-encoded zlib-compressed binary
 *     double-precision data using the functions from CompressData.
 *
 * In parallel, each process writes its contributions into an XML file of the
 * form `P0000001.fld` (where 1 is replaced by the rank of the process) inside a
 * directory with the desired output name. These files only include the
 * `ELEMENT` data. Metadata are instead stored in a separate `Info.xml` file,
 * which contains the Metadata and additional tags of the form
 *
 * `<Partition FileName="P0000000.fld"> ID list </Partition>`
 *
 * The ID list enumerates all element IDs on each partition's contribution. For
 * large parallel jobs, this is used to avoid each process reading in every
 * single partition in order to load field data.
 */
class FieldIOXml : public FieldIO
{
public:
    /// Creates an instance of this class
    LIB_UTILITIES_EXPORT static FieldIOSharedPtr create(
        LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    {
        return MemoryManager<FieldIOXml>::AllocateSharedPtr(pComm,
                                                            sharedFilesystem);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT FieldIOXml(
        LibUtilities::CommSharedPtr pComm,
        bool sharedFilesystem);

    LIB_UTILITIES_EXPORT virtual ~FieldIOXml()
    {
    }

    LIB_UTILITIES_EXPORT void ImportFieldDefs(
        DataSourceSharedPtr dataSource,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        bool expChild);

    LIB_UTILITIES_EXPORT void ImportFieldData(
        DataSourceSharedPtr dataSource,
        const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata);

    LIB_UTILITIES_EXPORT void WriteMultiFldFileIDs(
        const std::string &outfile,
        const std::vector<std::string> fileNames,
        std::vector<std::vector<unsigned int> > &elementList,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap);

    LIB_UTILITIES_EXPORT void SetUpFieldMetaData(
        const std::string &outname,
        const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        const FieldMetaDataMap &fieldmetadatamap);

    LIB_UTILITIES_EXPORT void ImportMultiFldFileIDs(
        const std::string &inFile,
        std::vector<std::string> &fileNames,
        std::vector<std::vector<unsigned int> > &elementList,
        FieldMetaDataMap &fieldmetadatamap);

    LIB_UTILITIES_EXPORT void v_Import(
        const std::string &infilename,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata =
            NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const Array<OneD, int> &ElementIDs = NullInt1DArray);

    /// Returns the class name.
    inline virtual const std::string &GetClassName() const
    {
        return className;
    }

private:
    LIB_UTILITIES_EXPORT virtual void v_Write(
        const std::string &outFile,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const bool backup = false);

    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        const std::string &filename, FieldMetaDataMap &fieldmetadatamap);
};

}
}
#endif
