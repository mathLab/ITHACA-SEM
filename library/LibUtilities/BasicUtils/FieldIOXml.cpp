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
//  Description: I/O routines relating to Fields into XML
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <LibUtilities/BasicUtils/CompressData.h>

#include "zlib.h"

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

// Buffer size for zlib compression/decompression
#define CHUNK 16384

namespace Nektar
{
    namespace LibUtilities
    {

        XmlDataSource::XmlDataSource(TiXmlDocument& doc) :
                m_doc(&doc)
        {
        }

        XmlDataSource::XmlDataSource(const std::string& fn) :
                m_doc(new TiXmlDocument(fn))
        {
            bool loadOkay = m_doc->LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << fn << std::endl;
            errstr << "Reason: " << m_doc->ErrorDesc() << std::endl;
            errstr << "Position: Line " << m_doc->ErrorRow() << ", Column "
                    << m_doc->ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());
        }
        XmlDataSource::~XmlDataSource()
        {
            delete m_doc;
        }

        TiXmlDocument& XmlDataSource::Get()
        {
            return *m_doc;
        }
        const TiXmlDocument& XmlDataSource::Get() const
        {
            return *m_doc;
        }

        DataSourceSharedPtr XmlDataSource::create(const std::string& fn)
        {
            return DataSourceSharedPtr(new XmlDataSource(fn));
        }
        DataSourceSharedPtr XmlDataSource::create(TiXmlDocument& doc)
        {
            return DataSourceSharedPtr(new XmlDataSource(doc));
        }

        typedef boost::shared_ptr<XmlDataSource> XmlDataSourceSharedPtr;

        std::string FieldIOXml::className =
                GetFieldIOFactory().RegisterCreatorFunction("Xml",
                        FieldIOXml::create, "XML-based output of field data.");

        FieldIOXml::FieldIOXml(LibUtilities::CommSharedPtr pComm,
                               bool sharedFilesystem) :
            FieldIO(pComm, sharedFilesystem)
        {
        }

        /**
         *
         */
        void FieldIOXml::v_Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            std::stringstream prfx;
            prfx << m_comm->GetRank() << ": FieldIOXml::v_Write(): ";
	    double tm0 = 0.0, tm1 = 0.0;
            if (0 == m_comm->GetRank())
	    {
	        cout << prfx.str() << "entering..." << endl;
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

            // Prepare to write out data. In parallel, we must create directory
            // and determine the full pathname to the file to write out.
            // Any existing file/directory which is in the way is removed.
            std::string filename = SetUpOutput(outFile);
            SetUpFieldMetaData(outFile, fielddefs, fieldmetadatamap);

            // Create the file (partition)
            TiXmlDocument doc;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            TiXmlElement * root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            AddInfoTag(root, fieldmetadatamap);

            for (int f = 0; f < fielddefs.size(); ++f)
            {
                //---------------------------------------------
                // Write ELEMENTS
                TiXmlElement * elemTag = new TiXmlElement("ELEMENTS");
                root->LinkEndChild(elemTag);

                // Write FIELDS
                std::string fieldsString;
                {
                    std::stringstream fieldsStringStream;
                    bool first = true;
                    for (std::vector<int>::size_type i = 0;
                            i < fielddefs[f]->m_fields.size(); i++)
                    {
                        if (!first)
                            fieldsStringStream << ",";
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
                    shapeStringStream
                            << ShapeTypeMap[fielddefs[f]->m_shapeType];
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
                            i < fielddefs[f]->m_basis.size(); i++)
                    {
                        if (!first)
                            basisStringStream << ",";
                        basisStringStream
                                << BasisTypeMap[fielddefs[f]->m_basis[i]];
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
                        for (int i = 0; i < fielddefs[f]->m_numHomogeneousDir;
                                ++i)
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
                            for (int i = 0;
                                    i < fielddefs[f]->m_homogeneousYIDs.size();
                                    i++)
                            {
                                if (!first)
                                    homoYIDsStringStream << ",";
                                homoYIDsStringStream
                                        << fielddefs[f]->m_homogeneousYIDs[i];
                                first = false;
                            }
                            homoYIDsString = homoYIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSYIDS",
                                homoYIDsString);
                    }

                    if (fielddefs[f]->m_homogeneousZIDs.size() > 0)
                    {
                        std::string homoZIDsString;
                        {
                            std::stringstream homoZIDsStringStream;
                            bool first = true;
                            for (int i = 0;
                                    i < fielddefs[f]->m_homogeneousZIDs.size();
                                    i++)
                            {
                                if (!first)
                                    homoZIDsStringStream << ",";
                                homoZIDsStringStream
                                        << fielddefs[f]->m_homogeneousZIDs[i];
                                first = false;
                            }
                            homoZIDsString = homoZIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSZIDS",
                                homoZIDsString);
                    }

                    if(fielddefs[f]->m_homogeneousSIDs.size() > 0)
                    {
                        std::string homoSIDsString;
                        {
                            std::stringstream homoSIDsStringStream;
                            bool first = true;
                            for(int i = 0; i < fielddefs[f]->m_homogeneousSIDs.size(); i++)
                            {
                                if (!first)
                                {
                                    homoSIDsStringStream << ",";
                                }
                                homoSIDsStringStream << fielddefs[f]->m_homogeneousSIDs[i];
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
                                i < fielddefs[f]->m_basis.size(); i++)
                        {
                            if (!first)
                                numModesStringStream << ",";
                            numModesStringStream << fielddefs[f]->m_numModes[i];
                            first = false;
                        }
                    }
                    else
                    {
                        numModesStringStream << "MIXORDER:";
                        bool first = true;
                        for (std::vector<int>::size_type i = 0;
                                i < fielddefs[f]->m_numModes.size(); i++)
                        {
                            if (!first)
                                numModesStringStream << ",";
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
                    GenerateSeqString(fielddefs[f]->m_elementIDs, idString);
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
                ASSERTL0(Z_OK == CompressData::ZlibEncodeToBase64Str(
                                                fielddata[f], base64string),
                         "Failed to compress field data.");

                elemTag->LinkEndChild(new TiXmlText(base64string));

            }
            doc.SaveFile(filename);

	    m_comm->Block();
	    // all data has been written
	    
            if (0 == m_comm->GetRank())
	    {
	        tm1 = m_comm->Wtime();
	        cout << prfx.str() << "leaving after " << tm1-tm0 << " secs..." << endl;
            }
        }

        /**
         *
         */
        void FieldIOXml::v_ImportFile(const std::string& fname,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                DataSourceSharedPtr dataSource)
        {
            if (!dataSource)
                dataSource = XmlDataSource::create(fname);

            ImportFieldDefs(dataSource, fielddefs, false);
            if (fielddata != NullVectorNekDoubleVector)
            {
                ImportFieldData(dataSource, fielddefs, fielddata);
            }
        }

        DataSourceSharedPtr FieldIOXml::v_ImportFieldMetaData(
                std::string filename, FieldMetaDataMap &fieldmetadatamap)
        {
            DataSourceSharedPtr ans = XmlDataSource::create(filename);
            v_ImportFieldMetaData(ans, fieldmetadatamap);
            return ans;
        }

        void FieldIOXml::v_ImportFieldMetaData(DataSourceSharedPtr dataSource,
                FieldMetaDataMap &fieldmetadatamap)
        {
            XmlDataSourceSharedPtr xml = boost::static_pointer_cast
                    < XmlDataSource > (dataSource);
            TiXmlElement* master = 0; // Master tag within which all data is contained.
            TiXmlElement* metadata = 0;

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
                        ASSERTL0(false,
                                "PARAM not provided as an attribute in FIELDMETADATA section");
                    }

                    // Now read body of param
                    std::string paramBodyStr;

                    TiXmlNode *paramBody = param->FirstChild();

                    paramBodyStr += paramBody->ToText()->Value();

                    fieldmetadatamap[paramString] = paramBodyStr;
                    param = param->NextSiblingElement("P");
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
                        TiXmlNode *paramBody = param->FirstChild();
                        std::string paramBodyStr = paramBody->ToText()->Value();

                        fieldmetadatamap[paramString] = paramBodyStr;
                    }
                    param = param->NextSiblingElement();
                }
            }

        }

        /**
         * The bool decides if the FieldDefs are in <EXPANSIONS> or in <NEKTAR>.
         */
        void FieldIOXml::ImportFieldDefs(DataSourceSharedPtr dataSource,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                bool expChild)
        {
            XmlDataSourceSharedPtr xml = boost::static_pointer_cast
                    < XmlDataSource > (dataSource);
            TiXmlElement* master = NULL; // Master tag within which all data is contained.

            master = xml->Get().FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");
            std::string strLoop = "NEKTAR";
            TiXmlElement* loopXml = master;

            TiXmlElement *expansionTypes;
            if (expChild)
            {
                expansionTypes = master->FirstChildElement("EXPANSIONS");
                ASSERTL0(expansionTypes,
                        "Unable to find EXPANSIONS tag in file.");
                loopXml = expansionTypes;
                strLoop = "EXPANSIONS";
            }

            // Loop through all nektar tags, finding all of the element tags.
            while (loopXml)
            {
                TiXmlElement* element = loopXml->FirstChildElement("ELEMENTS");
                ASSERTL0(element,
                        "Unable to find ELEMENTS tag within nektar tag.");

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
                    bool pointDef = false;
                    bool numPointDef = false;
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
                            homoSIDsString.insert(0,attr->Value());
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
                            if(!boost::iequals(attr->Value(),
                                               CompressData::GetCompressString()))
                            {
                                WARNINGL0(false, "Compressed formats do not "
                                          "match. Expected: "
                                          + CompressData::GetCompressString()
                                          + " but got "+ string(attr->Value()));
                            }
                        }
                        else if (attrName =="BITSIZE")
                        {
                            // This information is for future compatibility
                            // issues, for example in case we end up using a 128
                            // bit machine. Currently just do nothing.
                        }
                        else
                        {
                            std::string errstr("Unknown attribute: ");
                            errstr += attrName;
                            ASSERTL1(false, errstr.c_str());
                        }

                        // Get the next attribute.
                        attr = attr->Next();
                    }

                    // Check to see if using strips formulation
                    bool strips = false;
                    if(shapeString.find("Strips")!=string::npos)
                    {
                        strips = true;
                    }

                    // Check to see if homogeneous expansion and if so
                    // strip down the shapeString definition
                    int numHomoDir = 0;
                    size_t loc;
                    //---> This finds the first location of  'n'!
                    if ((loc = shapeString.find_first_of("-")) != string::npos)
                    {
                        if (shapeString.find("Exp1D") != string::npos)
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
                        bool valid = ParseUtils::GenerateSeqVector(
                                idString.c_str(), elementIds);
                        ASSERTL0(valid,
                                "Unable to correctly parse the element ids.");
                    }

                    // Get the geometrical shape
                    ShapeType shape;
                    bool valid = false;
                    for (unsigned int j = 0; j < SIZE_ShapeType; j++)
                    {
                        if (ShapeTypeMap[j] == shapeString)
                        {
                            shape = (ShapeType) j;
                            valid = true;
                            break;
                        }
                    }

                    ASSERTL0(valid,
                            std::string(
                                    "Unable to correctly parse the shape type: ").append(
                                    shapeString).c_str());

                    // Get the basis
                    std::vector < std::string > basisStrings;
                    std::vector<BasisType> basis;
                    valid = ParseUtils::GenerateOrderedStringVector(
                            basisString.c_str(), basisStrings);
                    ASSERTL0(valid,
                            "Unable to correctly parse the basis types.");
                    for (std::vector<std::string>::size_type i = 0;
                            i < basisStrings.size(); i++)
                    {
                        valid = false;
                        for (unsigned int j = 0; j < SIZE_BasisType; j++)
                        {
                            if (BasisTypeMap[j] == basisStrings[i])
                            {
                                basis.push_back((BasisType) j);
                                valid = true;
                                break;
                            }
                        }
                        ASSERTL0(valid,
                                std::string(
                                        "Unable to correctly parse the basis type: ").append(
                                        basisStrings[i]).c_str());
                    }

                    // Get homoLengths
                    std::vector<NekDouble> homoLengths;
                    if (numHomoDir)
                    {
                        valid = ParseUtils::GenerateUnOrderedVector(
                                homoLengthsString.c_str(), homoLengths);
                        ASSERTL0(valid,
                                "Unable to correctly parse the number of homogeneous lengths.");
                    }

                    // Get Homogeneous strips IDs
                    std::vector<unsigned int> homoSIDs;
                    if(strips)
                    {
                        valid = ParseUtils::GenerateSeqVector(homoSIDsString.c_str(), homoSIDs);
                        ASSERTL0(valid, "Unable to correctly parse homogeneous strips IDs.");
                    }

                    // Get Homogeneous points IDs
                    std::vector<unsigned int> homoZIDs;
                    std::vector<unsigned int> homoYIDs;

                    if (numHomoDir == 1)
                    {
                        valid = ParseUtils::GenerateSeqVector(
                                homoZIDsString.c_str(), homoZIDs);
                        ASSERTL0(valid,
                                "Unable to correctly parse homogeneous planes IDs.");
                    }

                    if (numHomoDir == 2)
                    {
                        valid = ParseUtils::GenerateSeqVector(
                                homoZIDsString.c_str(), homoZIDs);
                        ASSERTL0(valid,
                                "Unable to correctly parse homogeneous lines IDs in z-direction.");
                        valid = ParseUtils::GenerateSeqVector(
                                homoYIDsString.c_str(), homoYIDs);
                        ASSERTL0(valid,
                                "Unable to correctly parse homogeneous lines IDs in y-direction.");
                    }

                    // Get points type
                    std::vector<PointsType> points;

                    if (pointDef)
                    {
                        std::vector < std::string > pointsStrings;
                        valid = ParseUtils::GenerateOrderedStringVector(
                                pointsString.c_str(), pointsStrings);
                        ASSERTL0(valid,
                                "Unable to correctly parse the points types.");
                        for (std::vector<std::string>::size_type i = 0;
                                i < pointsStrings.size(); i++)
                        {
                            valid = false;
                            for (unsigned int j = 0; j < SIZE_PointsType; j++)
                            {
                                if (kPointsTypeStr[j] == pointsStrings[i])
                                {
                                    points.push_back((PointsType) j);
                                    valid = true;
                                    break;
                                }
                            }

                            ASSERTL0(valid,
                                    std::string(
                                            "Unable to correctly parse the points type: ").append(
                                            pointsStrings[i]).c_str());
                        }
                    }

                    // Get numModes
                    std::vector<unsigned int> numModes;
                    bool UniOrder = false;

                    if (strstr(numModesString.c_str(), "UNIORDER:"))
                    {
                        UniOrder = true;
                    }

                    valid = ParseUtils::GenerateOrderedVector(
                            numModesString.c_str() + 9, numModes);
                    ASSERTL0(valid,
                            "Unable to correctly parse the number of modes.");

                    // Get numPoints
                    std::vector<unsigned int> numPoints;
                    if (numPointDef)
                    {
                        valid = ParseUtils::GenerateOrderedVector(
                                numPointsString.c_str(), numPoints);
                        ASSERTL0(valid,
                                "Unable to correctly parse the number of points.");
                    }

                    // Get fields names
                    std::vector < std::string > Fields;
                    valid = ParseUtils::GenerateOrderedStringVector(
                            fieldsString.c_str(), Fields);
                    ASSERTL0(valid,
                            "Unable to correctly parse the number of fields.");

                    FieldDefinitionsSharedPtr fielddef = MemoryManager<
                        FieldDefinitions>::AllocateSharedPtr(
                            shape, elementIds, basis, UniOrder, numModes,
                            Fields, numHomoDir, homoLengths, strips, homoSIDs,
                            homoZIDs, homoYIDs, points, pointDef, numPoints,
                            numPointDef);

                    fielddefs.push_back(fielddef);

                    element = element->NextSiblingElement("ELEMENTS");
                }
                loopXml = loopXml->NextSiblingElement(strLoop);
            }
        }

        /**
         *
         */
        void FieldIOXml::ImportFieldData(DataSourceSharedPtr dataSource,
                const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata)
        {
            int cntdumps = 0;
            XmlDataSourceSharedPtr xml = boost::static_pointer_cast
                    < XmlDataSource > (dataSource);

            TiXmlElement* master = NULL; // Master tag within which all data is contained.

            master = xml->Get().FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");

            // Loop through all nektar tags, finding all of the element tags.
            while (master)
            {
                TiXmlElement* element = master->FirstChildElement("ELEMENTS");
                ASSERTL0(element,
                        "Unable to find ELEMENTS tag within nektar tag.");
                while (element)
                {
                    // Extract the body, which the "data".
                    TiXmlNode* elementChild = element->FirstChild();
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
                    if(CompressStr)
                    {
                        if(!boost::iequals(CompressStr,
                                           CompressData::GetCompressString()))
                        {
                            WARNINGL0(false, "Compressed formats do not match. "
                                      "Expected: "
                                      + CompressData::GetCompressString()
                                      + " but got "+ string(CompressStr));
                        }
                    }

                    ASSERTL0(Z_OK == CompressData::ZlibDecodeFromBase64Str(
                                                        elementStr,
                                                        elementFieldData),
                             "Failed to decompress field data.");
                    fielddata.push_back(elementFieldData);

                    int datasize = CheckFieldDefinition(fielddefs[cntdumps]);
                    ASSERTL0(
                            fielddata[cntdumps].size()
                                    == datasize
                                            * fielddefs[cntdumps]->m_fields.size(),
                            "Input data is not the same length as header infoarmation");

                    cntdumps++;

                    element = element->NextSiblingElement("ELEMENTS");
                }
                master = master->NextSiblingElement("NEKTAR");
            }
        }

    }
}
