////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIO.cpp
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
//  Description: I/O routines relating to Fields into XML
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <LibUtilities/BasicUtils/CompressData.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <boost/format.hpp>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

namespace berrc = boost::system::errc;

namespace Nektar
{
namespace LibUtilities
{

std::string FieldIOXml::className = GetFieldIOFactory().RegisterCreatorFunction(
    "Xml", FieldIOXml::create, "XML-based output of field data.");

/**
 * @brief Default constructor.
 *
 * @param pComm              Communicator.
 * @param sharedFilesystem   True if the underlying filesystem is shared by the
 *                           compute nodes.
 */
FieldIOXml::FieldIOXml(LibUtilities::CommSharedPtr pComm, bool sharedFilesystem)
    : FieldIO(pComm, sharedFilesystem)
{
}

/**
 * @brief Write an XML file to @p outFile given the field definitions @p
 * fielddefs, field data @p fielddata and metadata @p fieldmetadatamap.
 *
 * The writing strategy is as follows:
 *
 *   - Use FieldIO::SetUpOutput to construct the directory to contain each
 *     partition.
 *   - The root processor writes an `Info.xml` file containing the field
 *     metadata and an index that describes which elements lie in which XML
 *     file.
 *   - Each processor then writes an XML file containing the field definitions
 *     for that processor and output data in base64-encoded zlib-compressed
 *     format.
 *
 * @param outFile           Output filename.
 * @param fielddefs         Input field definitions.
 * @param fielddata         Input field data.
 * @param fieldmetadatamap  Field metadata.
 */
void FieldIOXml::v_Write(const std::string &outFile,
                         std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                         std::vector<std::vector<NekDouble> > &fielddata,
                         const FieldMetaDataMap &fieldmetadatamap,
                         const bool backup)
{
    double tm0 = 0.0, tm1 = 0.0;
    if (m_comm->TreatAsRankZero())
    {
        tm0 = m_comm->Wtime();
    }

    // Check everything seems sensible
    ASSERTL1(fielddefs.size() == fielddata.size(),
             "Length of fielddefs and fielddata incompatible");
    for (int f = 0; f < fielddefs.size(); ++f)
    {
        ASSERTL1(fielddata[f].size() > 0,
                 "Fielddata vector must contain at least one value.");

        ASSERTL1(fielddata[f].size() ==
                     fielddefs[f]->m_fields.size() *
                         CheckFieldDefinition(fielddefs[f]),
                 "Invalid size of fielddata vector.");
    }

    // Prepare to write out data. In parallel, we must create directory and
    // determine the full pathname to the file to write out.  Any existing
    // file/directory which is in the way is removed.
    std::string filename = SetUpOutput(outFile, true, backup);
    SetUpFieldMetaData(outFile, fielddefs, fieldmetadatamap);

    // Create the file (partition)
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    AddInfoTag(XmlTagWriterSharedPtr(new XmlTagWriter(root)), fieldmetadatamap);

    for (int f = 0; f < fielddefs.size(); ++f)
    {
        //---------------------------------------------
        // Write ELEMENTS
        TiXmlElement *elemTag = new TiXmlElement("ELEMENTS");
        root->LinkEndChild(elemTag);

        // Write FIELDS
        std::string fieldsString;
        {
            std::stringstream fieldsStringStream;
            bool first = true;
            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_fields.size();
                 i++)
            {
                if (!first)
                {
                    fieldsStringStream << ",";
                }
                fieldsStringStream << fielddefs[f]->m_fields[i];
                first = false;
            }
            fieldsString = fieldsStringStream.str();
        }
        elemTag->SetAttribute("FIELDS", fieldsString);

        // Write SHAPE
        std::string shapeString;
        {
            std::stringstream shapeStringStream;
            shapeStringStream << ShapeTypeMap[fielddefs[f]->m_shapeType];
            if (fielddefs[f]->m_numHomogeneousDir == 1)
            {
                shapeStringStream << "-HomogenousExp1D";
            }
            else if (fielddefs[f]->m_numHomogeneousDir == 2)
            {
                shapeStringStream << "-HomogenousExp2D";
            }

            if (fielddefs[f]->m_homoStrips)
            {
                shapeStringStream << "-Strips";
            }

            shapeString = shapeStringStream.str();
        }
        elemTag->SetAttribute("SHAPE", shapeString);

        // Write BASIS
        std::string basisString;
        {
            std::stringstream basisStringStream;
            bool first = true;
            for (std::vector<BasisType>::size_type i = 0;
                 i < fielddefs[f]->m_basis.size();
                 i++)
            {
                if (!first)
                {
                    basisStringStream << ",";
                }
                basisStringStream << BasisTypeMap[fielddefs[f]->m_basis[i]];
                first = false;
            }
            basisString = basisStringStream.str();
        }
        elemTag->SetAttribute("BASIS", basisString);

        // Write homogeneuous length details
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            std::string homoLenString;
            {
                std::stringstream homoLenStringStream;
                bool first = true;
                for (int i = 0; i < fielddefs[f]->m_numHomogeneousDir; ++i)
                {
                    if (!first)
                        homoLenStringStream << ",";
                    homoLenStringStream
                        << fielddefs[f]->m_homogeneousLengths[i];
                    first = false;
                }
                homoLenString = homoLenStringStream.str();
            }
            elemTag->SetAttribute("HOMOGENEOUSLENGTHS", homoLenString);
        }

        // Write homogeneuous planes/lines details
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            if (fielddefs[f]->m_homogeneousYIDs.size() > 0)
            {
                std::string homoYIDsString;
                {
                    std::stringstream homoYIDsStringStream;
                    bool first = true;
                    for (int i = 0; i < fielddefs[f]->m_homogeneousYIDs.size();
                         i++)
                    {
                        if (!first)
                        {
                            homoYIDsStringStream << ",";
                        }
                        homoYIDsStringStream
                            << fielddefs[f]->m_homogeneousYIDs[i];
                        first = false;
                    }
                    homoYIDsString = homoYIDsStringStream.str();
                }
                elemTag->SetAttribute("HOMOGENEOUSYIDS", homoYIDsString);
            }

            if (fielddefs[f]->m_homogeneousZIDs.size() > 0)
            {
                std::string homoZIDsString;
                {
                    std::stringstream homoZIDsStringStream;
                    bool first = true;
                    for (int i = 0; i < fielddefs[f]->m_homogeneousZIDs.size();
                         i++)
                    {
                        if (!first)
                        {
                            homoZIDsStringStream << ",";
                        }
                        homoZIDsStringStream
                            << fielddefs[f]->m_homogeneousZIDs[i];
                        first = false;
                    }
                    homoZIDsString = homoZIDsStringStream.str();
                }
                elemTag->SetAttribute("HOMOGENEOUSZIDS", homoZIDsString);
            }

            if (fielddefs[f]->m_homogeneousSIDs.size() > 0)
            {
                std::string homoSIDsString;
                {
                    std::stringstream homoSIDsStringStream;
                    bool first = true;
                    for (int i = 0; i < fielddefs[f]->m_homogeneousSIDs.size();
                         i++)
                    {
                        if (!first)
                        {
                            homoSIDsStringStream << ",";
                        }
                        homoSIDsStringStream
                            << fielddefs[f]->m_homogeneousSIDs[i];
                        first = false;
                    }
                    homoSIDsString = homoSIDsStringStream.str();
                }
                elemTag->SetAttribute("HOMOGENEOUSSIDS", homoSIDsString);
            }
        }

        // Write NUMMODESPERDIR
        std::string numModesString;
        {
            std::stringstream numModesStringStream;

            if (fielddefs[f]->m_uniOrder)
            {
                numModesStringStream << "UNIORDER:";
                // Just dump single definition
                bool first = true;
                for (std::vector<int>::size_type i = 0;
                     i < fielddefs[f]->m_basis.size();
                     i++)
                {
                    if (!first)
                    {
                        numModesStringStream << ",";
                    }
                    numModesStringStream << fielddefs[f]->m_numModes[i];
                    first = false;
                }
            }
            else
            {
                numModesStringStream << "MIXORDER:";
                bool first = true;
                for (std::vector<int>::size_type i = 0;
                     i < fielddefs[f]->m_numModes.size();
                     i++)
                {
                    if (!first)
                    {
                        numModesStringStream << ",";
                    }
                    numModesStringStream << fielddefs[f]->m_numModes[i];
                    first = false;
                }
            }

            numModesString = numModesStringStream.str();
        }
        elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

        // Write ID
        // Should ideally look at ways of compressing this stream
        // if just sequential;
        std::string idString;
        {
            std::stringstream idStringStream;
            idString = ParseUtils::GenerateSeqString(fielddefs[f]->m_elementIDs);
        }
        elemTag->SetAttribute("ID", idString);
        elemTag->SetAttribute("COMPRESSED",
                              LibUtilities::CompressData::GetCompressString());

        // Add this information for future compatibility
        // issues, for exmaple in case we end up using a 128
        // bit machine.
        elemTag->SetAttribute("BITSIZE",
                              LibUtilities::CompressData::GetBitSizeStr());
        std::string base64string;
        ASSERTL0(Z_OK == CompressData::ZlibEncodeToBase64Str(fielddata[f],
                                                             base64string),
                 "Failed to compress field data.");

        elemTag->LinkEndChild(new TiXmlText(base64string));
    }
    doc.SaveFile(filename);

    m_comm->Block();

    // all data has been written
    if (m_comm->TreatAsRankZero())
    {
        tm1 = m_comm->Wtime();
        std::cout << " (" << tm1 - tm0 << "s, XML)" << std::endl;
    }
}

/**
 * @brief Write out a file containing element ID to partition mapping.
 *
 * This function writes out an XML file (usually called `Info.xml`) that
 * contains the element IDs that are contained within each partition, as well as
 * the field metadata map.
 *
 * @param outFile            Output multi-field file name.
 * @param fileNames          List of partition filenames.
 * @param elementList        Vector of element IDs that lie on each process.
 * @param fieldmetadatamap   Field metadata map that is written into @p outFile.
 */
void FieldIOXml::WriteMultiFldFileIDs(
    const std::string                       &outFile,
    const std::vector<std::string>           fileNames,
    std::vector<std::vector<unsigned int> > &elementList,
    const FieldMetaDataMap                  &fieldmetadatamap)
{
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    ASSERTL0(fileNames.size() == elementList.size(),
             "Outfile names and list of elements ids does not match");

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    AddInfoTag(XmlTagWriterSharedPtr(new XmlTagWriter(root)), fieldmetadatamap);

    for (int t = 0; t < fileNames.size(); ++t)
    {
        if (elementList[t].size())
        {
            TiXmlElement *elemIDs = new TiXmlElement("Partition");
            root->LinkEndChild(elemIDs);

            elemIDs->SetAttribute("FileName", fileNames[t]);

            std::string IDstr = ParseUtils::GenerateSeqString(elementList[t]);

            elemIDs->LinkEndChild(new TiXmlText(IDstr));
        }
    }

    doc.SaveFile(outFile);
}

/**
 * @brief Read file containing element ID to partition mapping.
 *
 * This function reads an XML file (usually called `Info.xml`) that contains the
 * element IDs that are contained within each partition, as well as the field
 * metadata map.
 *
 * @param inFile             Input multi-field file name.
 * @param fileNames          List of partition filenames.
 * @param elementList        Vector of element IDs that lie on each process.
 * @param fieldmetadatamap   Field metadata map that is read from @p inFile.
 */
void FieldIOXml::ImportMultiFldFileIDs(
    const std::string                       &inFile,
    std::vector<std::string>                &fileNames,
    std::vector<std::vector<unsigned int> > &elementList,
    FieldMetaDataMap                        &fieldmetadatamap)
{
    boost::ignore_unused(fieldmetadatamap);

    TiXmlDocument doc(inFile);
    bool loadOkay = doc.LoadFile();

    std::stringstream errstr;
    errstr << "Unable to load file: " << inFile << std::endl;
    errstr << "Reason: " << doc.ErrorDesc() << std::endl;
    errstr << "Position: Line " << doc.ErrorRow() << ", Column "
           << doc.ErrorCol() << std::endl;
    ASSERTL0(loadOkay, errstr.str());

    // Handle on XML document
    TiXmlHandle docHandle(&doc);

    // Retrieve main NEKTAR tag - XML specification states one
    // top-level element tag per file.
    TiXmlElement *master = doc.FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");

    // Partition element tag name
    std::string strPartition = "Partition";

    // First attempt to get the first Partition element
    TiXmlElement *fldfileIDs = master->FirstChildElement(strPartition.c_str());
    if (!fldfileIDs)
    {
        // If this files try previous name
        strPartition = "MultipleFldFiles";
        fldfileIDs   = master->FirstChildElement("MultipleFldFiles");
    }
    ASSERTL0(fldfileIDs,
             "Unable to find 'Partition' or 'MultipleFldFiles' tag "
             "within nektar tag.");

    while (fldfileIDs)
    {
        // Read file name of partition file
        const char *attr = fldfileIDs->Attribute("FileName");
        ASSERTL0(attr,
                 "'FileName' not provided as an attribute of '" + strPartition +
                     "' tag.");
        fileNames.push_back(std::string(attr));

        const char *elementIDs = fldfileIDs->GetText();
        ASSERTL0(elementIDs, "Element IDs not specified.");

        std::string elementIDsStr(elementIDs);

        std::vector<unsigned int> idvec;
        ParseUtils::GenerateSeqVector(elementIDsStr, idvec);

        elementList.push_back(idvec);

        fldfileIDs = fldfileIDs->NextSiblingElement(strPartition.c_str());
    }
}

/**
 * @brief Import an XML format file.
 *
 * @param finfilename       Input filename
 * @param fielddefs         Field definitions of resulting field
 * @param fielddata         Field data of resulting field
 * @param fieldinfomap      Field metadata of resulting field
 * @param ElementIDs        If specified, contains the list of element IDs on
 *                          this rank. The resulting field definitions will only
 *                          contain data for the element IDs specified in this
 *                          array.
 */
void FieldIOXml::v_Import(const std::string &infilename,
                          std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                          std::vector<std::vector<NekDouble> > &fielddata,
                          FieldMetaDataMap &fieldinfomap,
                          const Array<OneD, int> &ElementIDs)
{
    std::string infile = infilename;

    fs::path pinfilename(infilename);

    // Check to see whether infile is a directory and therefore read in parallel
    // or serial.
    if (fs::is_directory(pinfilename))
    {
        fs::path infofile("Info.xml");
        fs::path fullpath = pinfilename / infofile;
        infile            = PortablePath(fullpath);

        std::vector<std::string> filenames;
        std::vector<std::vector<unsigned int> > elementIDs_OnPartitions;

        ImportMultiFldFileIDs(
            infile, filenames, elementIDs_OnPartitions, fieldinfomap);

        // Load metadata
        ImportFieldMetaData(infile, fieldinfomap);

        if (ElementIDs == NullInt1DArray) // load all elements
        {
            for (int i = 0; i < filenames.size(); ++i)
            {
                fs::path pfilename(filenames[i]);
                fullpath                       = pinfilename / pfilename;
                std::string fname              = PortablePath(fullpath);
                DataSourceSharedPtr dataSource = XmlDataSource::create(fname);
                ImportFieldDefs(dataSource, fielddefs, false);
                if (fielddata != NullVectorNekDoubleVector)
                {
                    ImportFieldData(dataSource, fielddefs, fielddata);
                }
            }
        }
        else // only load relevant elements from partitions
        {
            int i, j;
            std::map<int, std::vector<int> > FileIDs;
            std::set<int> LoadFile;

            for (i = 0; i < elementIDs_OnPartitions.size(); ++i)
            {
                for (j = 0; j < elementIDs_OnPartitions[i].size(); ++j)
                {
                    FileIDs[elementIDs_OnPartitions[i][j]].push_back(i);
                }
            }

            for (i = 0; i < ElementIDs.size(); ++i)
            {
                auto it = FileIDs.find(ElementIDs[i]);
                if (it != FileIDs.end())
                {
                    for (j = 0; j < it->second.size(); ++j)
                    {
                        LoadFile.insert(it->second[j]);
                    }
                }
            }

            for (auto &iter : LoadFile)
            {
                fs::path pfilename(filenames[iter]);
                fullpath                       = pinfilename / pfilename;
                std::string fname              = PortablePath(fullpath);
                DataSourceSharedPtr dataSource = XmlDataSource::create(fname);
                ImportFieldDefs(dataSource, fielddefs, false);
                if (fielddata != NullVectorNekDoubleVector)
                {
                    ImportFieldData(dataSource, fielddefs, fielddata);
                }
            }
        }
    }
    else
    {
        // serial format case
        DataSourceSharedPtr doc = ImportFieldMetaData(infilename, fieldinfomap);
        ImportFieldDefs(doc, fielddefs, false);
        if (fielddata != NullVectorNekDoubleVector)
        {
            ImportFieldData(doc, fielddefs, fielddata);
        }
    }
}

/**
 * @brief Import field metadata from @p filename and return the data source
 * which wraps @p filename.
 *
 * @param filename          Input filename.
 * @param fieldmetadatamap  Resulting field metadata from @p dataSource.
 */
DataSourceSharedPtr FieldIOXml::v_ImportFieldMetaData(
    const std::string &filename, FieldMetaDataMap &fieldmetadatamap)
{
    DataSourceSharedPtr doc    = XmlDataSource::create(filename);
    XmlDataSourceSharedPtr xml = std::static_pointer_cast<XmlDataSource>(doc);
    TiXmlElement *metadata     = 0;
    TiXmlElement *master       = 0; // Master tag within which all data is
                                    // contained.

    master = xml->Get().FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");
    std::string strLoop = "NEKTAR";

    // Retain original metadata structure for backwards compatibility
    // TODO: Remove old metadata format
    metadata = master->FirstChildElement("FIELDMETADATA");
    if (metadata)
    {
        TiXmlElement *param = metadata->FirstChildElement("P");

        while (param)
        {
            TiXmlAttribute *paramAttr = param->FirstAttribute();
            std::string attrName(paramAttr->Name());
            std::string paramString;

            if (attrName == "PARAM")
            {
                paramString.insert(0, paramAttr->Value());
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "PARAM not provided as an attribute in "
                         "FIELDMETADATA section");
            }

            // Now read body of param
            std::string paramBodyStr;

            TiXmlNode *paramBody = param->FirstChild();

            paramBodyStr += paramBody->ToText()->Value();

            fieldmetadatamap[paramString] = paramBodyStr;
            param                         = param->NextSiblingElement("P");
        }
    }

    // New metadata format
    metadata = master->FirstChildElement("Metadata");
    if (metadata)
    {
        TiXmlElement *param = metadata->FirstChildElement();

        while (param)
        {
            std::string paramString = param->Value();
            if (paramString != "Provenance")
            {
                // Now read body of param
                if (param->NoChildren())
                {
                    fieldmetadatamap[paramString] = "";
                }
                else
                {
                    TiXmlNode *paramBody     = param->FirstChild();
                    std::string paramBodyStr = paramBody->ToText()->Value();
                    fieldmetadatamap[paramString] = paramBodyStr;
                }
            }
            param = param->NextSiblingElement();
        }
    }

    return doc;
}

/**
 * @brief Set up field meta data map.
 *
 * This routine sets up the necessary information for the field metadata map
 * before calling FieldIOXml::WriteMultiFldFileIDs, which involves each process
 * sending its element ID list to the root processor. The root processor writes
 * the `Info.xml` file.
 *
 * @param outname           Output directory.
 * @param fielddefs         Field definitions, needed to grab element IDs.
 * @param fieldmetadatamap  Field metadata map that is also written to the
 *                          `Info.xml` file.
 */
void FieldIOXml::SetUpFieldMetaData(
    const std::string &outname,
    const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    const FieldMetaDataMap &fieldmetadatamap)
{
    ASSERTL0(!outname.empty(), "Empty path given to SetUpFieldMetaData()");

    unsigned int nprocs = m_comm->GetSize();
    unsigned int rank   = m_comm->GetRank();

    fs::path specPath(outname);

    // Compute number of elements on this process and share with other
    // processes. Also construct list of elements on this process from
    // available vector of field definitions.
    std::vector<size_t>       elmtnums(nprocs, 0);
    std::vector<unsigned int> idlist;
    for (size_t i = 0; i < fielddefs.size(); ++i)
    {
        elmtnums[rank] += fielddefs[i]->m_elementIDs.size();
        idlist.insert(idlist.end(),
                      fielddefs[i]->m_elementIDs.begin(),
                      fielddefs[i]->m_elementIDs.end());
    }
    m_comm->AllReduce(elmtnums, LibUtilities::ReduceMax);

    // Collate per-process element lists on root process to generate
    // the info file.
    if (rank == 0)
    {
        std::vector<std::vector<unsigned int> > ElementIDs(nprocs);

        // Populate the list of element ID lists from all processes
        ElementIDs[0] = idlist;
        for (size_t i = 1; i < nprocs; ++i)
        {
            if (elmtnums[i] > 0)
            {
                std::vector<unsigned int> tmp(elmtnums[i]);
                m_comm->Recv(i, tmp);
                ElementIDs[i] = tmp;
            }
        }

        // Set up output names
        std::vector<std::string> filenames;
        for (unsigned int i = 0; i < nprocs; ++i)
        {
            boost::format pad("P%1$07d.%2$s");
            pad % i % GetFileEnding();
            filenames.push_back(pad.str());
        }

        // Write the Info.xml file
        std::string infofile =
            LibUtilities::PortablePath(specPath / fs::path("Info.xml"));

        WriteMultiFldFileIDs(infofile, filenames, ElementIDs, fieldmetadatamap);
    }
    else
    {
        // Send this process's ID list to the root process
        if (elmtnums[rank] > 0)
        {
            m_comm->Send(0, idlist);
        }
    }

}

/**
 * @brief Import field definitions from the target file.
 *
 * @param dataSource  Target XML file
 * @param fielddefs   Output vector that will contain read field definitions.
 * @param expChild    Determines if the field definitions are defined by
 *                    `<EXPANSIONS>` or in `<NEKTAR>`.
 */
void FieldIOXml::ImportFieldDefs(
    DataSourceSharedPtr dataSource,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    bool expChild)
{
    XmlDataSourceSharedPtr xml =
        std::static_pointer_cast<XmlDataSource>(dataSource);
    TiXmlElement *master =
        NULL; // Master tag within which all data is contained.

    master = xml->Get().FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");
    std::string strLoop   = "NEKTAR";
    TiXmlElement *loopXml = master;

    TiXmlElement *expansionTypes;
    if (expChild)
    {
        expansionTypes = master->FirstChildElement("EXPANSIONS");
        ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");
        loopXml = expansionTypes;
        strLoop = "EXPANSIONS";
    }

    // Loop through all nektar tags, finding all of the element tags.
    while (loopXml)
    {
        TiXmlElement *element = loopXml->FirstChildElement("ELEMENTS");

        while (element)
        {
            // Extract the attributes.
            std::string idString;
            std::string shapeString;
            std::string basisString;
            std::string homoLengthsString;
            std::string homoSIDsString;
            std::string homoZIDsString;
            std::string homoYIDsString;
            std::string numModesString;
            std::string numPointsString;
            std::string fieldsString;
            std::string pointsString;
            bool pointDef        = false;
            bool numPointDef     = false;
            TiXmlAttribute *attr = element->FirstAttribute();
            while (attr)
            {
                std::string attrName(attr->Name());
                if (attrName == "FIELDS")
                {
                    fieldsString.insert(0, attr->Value());
                }
                else if (attrName == "SHAPE")
                {
                    shapeString.insert(0, attr->Value());
                }
                else if (attrName == "BASIS")
                {
                    basisString.insert(0, attr->Value());
                }
                else if (attrName == "HOMOGENEOUSLENGTHS")
                {
                    homoLengthsString.insert(0, attr->Value());
                }
                else if (attrName == "HOMOGENEOUSSIDS")
                {
                    homoSIDsString.insert(0, attr->Value());
                }
                else if (attrName == "HOMOGENEOUSZIDS")
                {
                    homoZIDsString.insert(0, attr->Value());
                }
                else if (attrName == "HOMOGENEOUSYIDS")
                {
                    homoYIDsString.insert(0, attr->Value());
                }
                else if (attrName == "NUMMODESPERDIR")
                {
                    numModesString.insert(0, attr->Value());
                }
                else if (attrName == "ID")
                {
                    idString.insert(0, attr->Value());
                }
                else if (attrName == "POINTSTYPE")
                {
                    pointsString.insert(0, attr->Value());
                    pointDef = true;
                }
                else if (attrName == "NUMPOINTSPERDIR")
                {
                    numPointsString.insert(0, attr->Value());
                    numPointDef = true;
                }
                else if (attrName == "COMPRESSED")
                {
                    WARNINGL0(boost::iequals(attr->Value(),
                                             CompressData::GetCompressString()),
                              "Compressed formats do not "
                              "match. Expected: " +
                              CompressData::GetCompressString() +
                              " but got " + std::string(attr->Value()));
                }
                else if (attrName == "BITSIZE")
                {
                    // This information is for future compatibility
                    // issues, for example in case we end up using a 128
                    // bit machine. Currently just do nothing.
                }
                else
                {
                    std::string errstr("Unknown attribute: ");
                    errstr += attrName;
                    NEKERROR(ErrorUtil::ewarning, errstr.c_str());
                }

                // Get the next attribute.
                attr = attr->Next();
            }

            // Check to see if using strips formulation
            bool strips = false;
            if (shapeString.find("Strips") != std::string::npos)
            {
                strips = true;
            }

            // Check to see if homogeneous expansion and if so
            // strip down the shapeString definition
            int numHomoDir = 0;
            size_t loc;
            //---> This finds the first location of  'n'!
            if ((loc = shapeString.find_first_of("-")) != std::string::npos)
            {
                if (shapeString.find("Exp1D") != std::string::npos)
                {
                    numHomoDir = 1;
                }
                else // HomogeneousExp1D
                {
                    numHomoDir = 2;
                }

                shapeString.erase(loc, shapeString.length());
            }

            // Reconstruct the fielddefs.
            std::vector<unsigned int> elementIds;
            {
                bool valid =
                    ParseUtils::GenerateSeqVector(idString, elementIds);
                ASSERTL0(valid, "Unable to correctly parse the element ids.");
            }

            // Get the geometrical shape
            ShapeType shape = (ShapeType)0;
            bool valid = false;
            for (unsigned int j = 0; j < SIZE_ShapeType; j++)
            {
                if (ShapeTypeMap[j] == shapeString)
                {
                    shape = (ShapeType)j;
                    valid = true;
                    break;
                }
            }

            ASSERTL0(valid,
                     std::string("Unable to correctly parse the shape type: ")
                         .append(shapeString)
                         .c_str());

            // Get the basis
            std::vector<std::string> basisStrings;
            std::vector<BasisType> basis;
            valid = ParseUtils::GenerateVector(basisString, basisStrings);
            ASSERTL0(valid, "Unable to correctly parse the basis types.");
            for (std::vector<std::string>::size_type i = 0;
                 i < basisStrings.size();
                 i++)
            {
                valid = false;
                for (unsigned int j = 0; j < SIZE_BasisType; j++)
                {
                    if (BasisTypeMap[j] == basisStrings[i])
                    {
                        basis.push_back((BasisType)j);
                        valid = true;
                        break;
                    }
                }
                ASSERTL0(
                    valid,
                    std::string("Unable to correctly parse the basis type: ")
                        .append(basisStrings[i])
                        .c_str());
            }

            // Get homoLengths
            std::vector<NekDouble> homoLengths;
            if (numHomoDir)
            {
                valid = ParseUtils::GenerateVector(homoLengthsString,
                                                   homoLengths);
                ASSERTL0(valid, "Unable to correctly parse the number of "
                                "homogeneous lengths.");
            }

            // Get Homogeneous strips IDs
            std::vector<unsigned int> homoSIDs;
            if (strips)
            {
                valid = ParseUtils::GenerateVector(homoSIDsString, homoSIDs);
                ASSERTL0(valid,
                         "Unable to correctly parse homogeneous strips IDs.");
            }

            // Get Homogeneous points IDs
            std::vector<unsigned int> homoZIDs;
            std::vector<unsigned int> homoYIDs;

            if (numHomoDir == 1)
            {
                valid = ParseUtils::GenerateSeqVector(homoZIDsString,
                                                      homoZIDs);
                ASSERTL0(valid,
                         "Unable to correctly parse homogeneous planes IDs.");
            }

            if (numHomoDir == 2)
            {
                valid = ParseUtils::GenerateSeqVector(homoZIDsString,
                                                      homoZIDs);
                ASSERTL0(valid, "Unable to correctly parse homogeneous lines "
                                "IDs in z-direction.");
                valid = ParseUtils::GenerateSeqVector(homoYIDsString,
                                                      homoYIDs);
                ASSERTL0(valid, "Unable to correctly parse homogeneous lines "
                                "IDs in y-direction.");
            }

            // Get points type
            std::vector<PointsType> points;

            if (pointDef)
            {
                std::vector<std::string> pointsStrings;
                valid = ParseUtils::GenerateVector(pointsString, pointsStrings);
                ASSERTL0(valid, "Unable to correctly parse the points types.");
                for (std::vector<std::string>::size_type i = 0;
                     i < pointsStrings.size();
                     i++)
                {
                    valid = false;
                    for (unsigned int j = 0; j < SIZE_PointsType; j++)
                    {
                        if (kPointsTypeStr[j] == pointsStrings[i])
                        {
                            points.push_back((PointsType)j);
                            valid = true;
                            break;
                        }
                    }

                    ASSERTL0(valid,
                             std::string(
                                 "Unable to correctly parse the points type: ")
                                 .append(pointsStrings[i])
                                 .c_str());
                }
            }

            // Get numModes
            std::vector<unsigned int> numModes;
            bool UniOrder = false;

            if (strstr(numModesString.c_str(), "UNIORDER:"))
            {
                UniOrder = true;
            }

            valid = ParseUtils::GenerateVector(
                numModesString.substr(9), numModes);
            ASSERTL0(valid, "Unable to correctly parse the number of modes.");

            // Get numPoints
            std::vector<unsigned int> numPoints;
            if (numPointDef)
            {
                valid = ParseUtils::GenerateVector(numPointsString, numPoints);
                ASSERTL0(valid,
                         "Unable to correctly parse the number of points.");
            }

            // Get fields names
            std::vector<std::string> Fields;
            valid = ParseUtils::GenerateVector(fieldsString, Fields);
            ASSERTL0(valid, "Unable to correctly parse the number of fields.");

            FieldDefinitionsSharedPtr fielddef =
                MemoryManager<FieldDefinitions>::AllocateSharedPtr(shape,
                                                                   elementIds,
                                                                   basis,
                                                                   UniOrder,
                                                                   numModes,
                                                                   Fields,
                                                                   numHomoDir,
                                                                   homoLengths,
                                                                   strips,
                                                                   homoSIDs,
                                                                   homoZIDs,
                                                                   homoYIDs,
                                                                   points,
                                                                   pointDef,
                                                                   numPoints,
                                                                   numPointDef);

            fielddefs.push_back(fielddef);

            element = element->NextSiblingElement("ELEMENTS");
        }
        loopXml = loopXml->NextSiblingElement(strLoop);
    }
}

/**
 * @brief Import field data from a target file.
 *
 * @param dataSource  Target XML file
 * @param fielddefs   Field definitions for file
 * @param fielddata   On return, contains field data for each field.
 */
void FieldIOXml::ImportFieldData(
    DataSourceSharedPtr dataSource,
    const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata)
{
    int cntdumps = 0;
    XmlDataSourceSharedPtr xml =
        std::static_pointer_cast<XmlDataSource>(dataSource);

    TiXmlElement *master =
        NULL; // Master tag within which all data is contained.

    master = xml->Get().FirstChildElement("NEKTAR");
    ASSERTL0(master, "Unable to find NEKTAR tag in file.");

    // Loop through all nektar tags, finding all of the element tags.
    while (master)
    {
        TiXmlElement *element = master->FirstChildElement("ELEMENTS");
        ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");
        while (element)
        {
            // Extract the body, which the "data".
            TiXmlNode *elementChild = element->FirstChild();
            ASSERTL0(elementChild,
                     "Unable to extract the data from the element tag.");
            std::string elementStr;
            while (elementChild)
            {
                if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
                {
                    elementStr += elementChild->ToText()->ValueStr();
                }
                elementChild = elementChild->NextSibling();
            }

            std::vector<NekDouble> elementFieldData;

            // Convert from base64 to binary.
            const char *CompressStr = element->Attribute("COMPRESSED");
            if (CompressStr)
            {
                WARNINGL0(boost::iequals(CompressStr,
                                         CompressData::GetCompressString()),
                          "Compressed formats do not match. "
                          "Expected: " +
                          CompressData::GetCompressString() +
                          " but got " + std::string(CompressStr));
            }

            ASSERTL0(Z_OK == CompressData::ZlibDecodeFromBase64Str(
                                 elementStr, elementFieldData),
                     "Failed to decompress field data.");
            fielddata.push_back(elementFieldData);

            int datasize = CheckFieldDefinition(fielddefs[cntdumps]);
            ASSERTL0(
                fielddata[cntdumps].size() ==
                    datasize * fielddefs[cntdumps]->m_fields.size(),
                "Input data is not the same length as header infoarmation");

            cntdumps++;

            element = element->NextSiblingElement("ELEMENTS");
        }
        master = master->NextSiblingElement("NEKTAR");
    }
}
}
}
