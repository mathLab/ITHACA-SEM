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
//  Description: I/O routines relating to Fields
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIO.h>

namespace Nektar
{
    namespace LibUtilities
    {
        /**
         *
         */
        void Write(const std::string &outFile,
                   std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                   std::vector<std::vector<NekDouble> > &fielddata, 
                   FieldMetaDataMap &fieldmetadatamap)
        {
           ASSERTL1(fielddefs.size() == fielddata.size(),
                     "Length of fielddefs and fielddata incompatible");

            TiXmlDocument doc;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            cout << "Writing outfile: " << outFile << endl;

            TiXmlElement * root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            //---------------------------------------------
            // write field info section 
            if(fieldmetadatamap != NullFieldMetaDataMap)
            {
                TiXmlElement * infoTag = new TiXmlElement("FIELDMETADATA");
                root->LinkEndChild(infoTag);
                
                FieldMetaDataMap::iterator infoit;
                
                for(infoit = fieldmetadatamap.begin(); infoit != fieldmetadatamap.end(); ++infoit)
                {
                    NekDouble val = infoit->second;
                    stringstream s;
                    s << val; 
                    TiXmlElement * v = new TiXmlElement( "P" );
                    v->SetAttribute("PARAM",(infoit->first).c_str());
                    v->LinkEndChild(new TiXmlText(s.str()));
                    infoTag->LinkEndChild(v);
                }

            }

            for (int f = 0; f < fielddefs.size(); ++f)
            {

                ASSERTL1(fielddata[f].size() > 0,
                        "Fielddata vector must contain at least one value.");
                
                ASSERTL1(fielddata[f].size() ==
                             fielddefs[f]->m_fields.size() *
                             CheckFieldDefinition(fielddefs[f]),
                         "Invalid size of fielddata vector.");

                //---------------------------------------------
                // Write ELEMENTS
                TiXmlElement * elemTag = new TiXmlElement("ELEMENTS");
                root->LinkEndChild(elemTag);

                // Write FIELDS
                std::string fieldsString;
                {
                    std::stringstream fieldsStringStream;
                    bool first = true;
                    for (std::vector<int>::size_type i = 0; i
                    < fielddefs[f]->m_fields.size(); i++)
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
                    shapeStringStream << ShapeTypeMap[fielddefs[f]->m_shapeType];
                    if(fielddefs[f]->m_numHomogeneousDir == 1)
                    {
                        shapeStringStream << "-HomogenousExp1D";
                    }
                    else if (fielddefs[f]->m_numHomogeneousDir == 2)
                    {
                        shapeStringStream << "-HomogenousExp2D";
                    }

                    shapeString = shapeStringStream.str();
                }
                elemTag->SetAttribute("SHAPE", shapeString);

                // Write BASIS
                std::string basisString;
                {
                    std::stringstream basisStringStream;
                    bool first = true;
                    for (std::vector<BasisType>::size_type i = 0; i < fielddefs[f]->m_basis.size(); i++)
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
                if(fielddefs[f]->m_numHomogeneousDir)
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
                if(fielddefs[f]->m_numHomogeneousDir)
                {
                    if(fielddefs[f]->m_homogeneousYIDs.size() > 0)
                    {
                        std::string homoYIDsString;
                        {
                            std::stringstream homoYIDsStringStream;
                            bool first = true;
                            for(int i = 0; i < fielddefs[f]->m_homogeneousYIDs.size(); i++)
                            {
                                if (!first)
                                    homoYIDsStringStream << ",";
                                homoYIDsStringStream << fielddefs[f]->m_homogeneousYIDs[i];
                                first = false;
                            }
                            homoYIDsString = homoYIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSYIDS", homoYIDsString);
                    }
                    
                    if(fielddefs[f]->m_homogeneousZIDs.size() > 0)
                    {
                        std::string homoZIDsString;
                        {
                            std::stringstream homoZIDsStringStream;
                            bool first = true;
                            for(int i = 0; i < fielddefs[f]->m_homogeneousZIDs.size(); i++)
                            {
                                if (!first)
                                    homoZIDsStringStream << ",";
                                homoZIDsStringStream << fielddefs[f]->m_homogeneousZIDs[i];
                                first = false;
                            }
                            homoZIDsString = homoZIDsStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSZIDS", homoZIDsString);
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
                        for (std::vector<int>::size_type i = 0; i
                                 < fielddefs[f]->m_basis.size(); i++)
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
                        for (std::vector<int>::size_type i = 0; i
                                 < fielddefs[f]->m_numModes.size(); i++)
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
                    bool first = true;
                    for (std::vector<unsigned int>::size_type i = 0; 
                         i < fielddefs[f]->m_elementIDs.size(); i++)
                    {
                        if (!first)
                            idStringStream << ",";
                        idStringStream << fielddefs[f]->m_elementIDs[i];
                        first = false;
                    }
                    idString = idStringStream.str();
                }
                elemTag->SetAttribute("ID", idString);

                // Write binary data
                std::string compressedDataString;
                {
                    // Serialize the fielddata vector to the stringstream.
                    std::stringstream archiveStringStream(std::string((char*) &fielddata[f][0], sizeof(fielddata[f][0]) / sizeof(char) * fielddata[f].size()));

                    // Compress the serialized data.
                    std::stringstream compressedData;
                    {
                        boost::iostreams::filtering_streambuf<boost::iostreams::input>   out;
                        out.push(boost::iostreams::zlib_compressor());
                        out.push(archiveStringStream);
                        boost::iostreams::copy(out, compressedData);

                        // try decompressing to see if compression accurate
                        try
                        {
                            std::stringstream elementDecompressedData;
                            boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
                            in.push(boost::iostreams::zlib_decompressor());
                            in.push(compressedData);
                            
                            boost::iostreams::copy(in, elementDecompressedData);
                        }
                        catch (boost::iostreams::zlib_error e)
                        {
                            if (e.error() == boost::iostreams::zlib::stream_end)
                            {
                                cout <<  "Stream end zlib error";
                            }
                            else if (e.error() == boost::iostreams::zlib::stream_error)
                            {
                                cout <<   "Stream zlib error";
                            }
                            else if (e.error() == boost::iostreams::zlib::version_error)
                            {
                                cout << "Version zlib error";
                            }
                            else if (e.error() == boost::iostreams::zlib::data_error)
                            {
                                cout <<  "Data zlib error";
                            }
                            else if (e.error() == boost::iostreams::zlib::mem_error)
                            {
                                cout << "Memory zlib error";
                            }
                            else if (e.error() == boost::iostreams::zlib::buf_error)
                            {
                                cout <<"Buffer zlib error";
                            }
                            else
                            {
                                cout <<  "Unknown zlib error";
                            }

                            cout << endl << "Retrying compression will not check again" << endl;
                            compressedData.clear();
                            compressedData.str("");
                            boost::iostreams::copy(out, compressedData);
                        }
                    }

                    // If the string length is not divisible by 3,
                    // pad it. There is a bug in transform_width
                    // that will make it reference past the end
                    // and crash.
                    switch (compressedData.str().length() % 3)
                    {
                    case 1:
                        compressedData << '\0';
                    case 2:
                        compressedData << '\0';
                        break;
                    }
                    compressedDataString = compressedData.str();
                }

                // Convert from binary to base64.
                typedef boost::archive::iterators::base64_from_binary<
                        boost::archive::iterators::transform_width<
                        std::string::const_iterator, 6, 8> > base64_t;
                std::string base64string(base64_t(compressedDataString.begin()),
                        base64_t(compressedDataString.end()));
                elemTag->LinkEndChild(new TiXmlText(base64string));

            }
            doc.SaveFile(outFile);
        }


        /**
         *
         */
        void Import(const std::string& infilename,
                    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    std::vector<std::vector<NekDouble> > &fielddata,
                    FieldMetaDataMap &fieldmetadatamap)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());


            //---------------------------------------------
            // read field meta data  section 
            ImportFieldMetaData(doc,fieldmetadatamap);
            ImportFieldDefs(doc, fielddefs, false);
            ImportFieldData(doc, fielddefs, fielddata);
        }


        void ImportFieldMetaData(std::string filename,
                                 FieldMetaDataMap &fieldmetadatamap)
        {
            TiXmlDocument doc(filename);
            bool loadOkay = doc.LoadFile();
            
            std::stringstream errstr;
            errstr << "Unable to load file: " << filename << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());
                    
            ImportFieldMetaData(doc,fieldmetadatamap);
        }
        

        void ImportFieldMetaData(TiXmlDocument &doc,
                                FieldMetaDataMap &fieldmetadatamap)
        {
            
            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.
            
            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");
            std::string strLoop = "NEKTAR";


            TiXmlElement* metadata = master->FirstChildElement("FIELDMETADATA");
            
            if(!metadata) // section not available so just exit
            {
                return;
            }
            else
            {
                
                TiXmlElement *param = metadata->FirstChildElement("P");
                    
                while (param)
                {
                    TiXmlAttribute *paramAttr = param->FirstAttribute();
                    std::string attrName(paramAttr->Name());
                    std::string paramString;
                    
                    if(attrName == "PARAM")
                    {
                        paramString.insert(0,paramAttr->Value());
                    }
                    else
                    {
                        ASSERTL0(false,"PARAM not provided as an attribute in FIELDMETADATA section");
                    }

                    // Now read body of param
                    std::string paramBodyStr;
                    
                    TiXmlNode *paramBody = param->FirstChild();
                    
                    paramBodyStr += paramBody->ToText()->Value();
                    
                    NekDouble value;
                    std::istringstream paramStrm(paramBodyStr.c_str());
                    
                    try
                    {
                        while(!paramStrm.fail())
                        {
                            paramStrm >> value;
                        }
                    }
                    catch(...)
                    {
                        ASSERTL0(false,"Failied to read PARAM data");
                    }

                    fieldmetadatamap[paramString] = value; 
                    param = param->NextSiblingElement("P");
                }
            }
        }
        

        /**
         * The bool decides if the FieldDefs are in <EXPANSIONS> or in <NEKTAR>.
         */
        void ImportFieldDefs(TiXmlDocument &doc, std::vector<FieldDefinitionsSharedPtr> &fielddefs, 
                             bool expChild)
        {
            ASSERTL1(fielddefs.size() == 0, "Expected an empty fielddefs vector.");

            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");
            std::string strLoop = "NEKTAR";
            TiXmlElement* loopXml = master;

            TiXmlElement *expansionTypes;
            if(expChild)
            {
                expansionTypes = master->FirstChildElement("EXPANSIONS");
                ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");
                loopXml = expansionTypes;
                strLoop = "EXPANSIONS";
            }

            // Loop through all nektar tags, finding all of the element tags.
            while (loopXml)
            {
                TiXmlElement* element = loopXml->FirstChildElement("ELEMENTS");
                ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");

                while (element)
                {
                    // Extract the attributes.
                    std::string idString;
                    std::string shapeString;
                    std::string basisString;
                    std::string homoLengthsString;
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
                            homoLengthsString.insert(0,attr->Value());
                        }
                        else if (attrName == "HOMOGENEOUSZIDS")
                        {
                            homoZIDsString.insert(0,attr->Value());
                        }
                        else if (attrName == "HOMOGENEOUSYIDS")
                        {
                            homoYIDsString.insert(0,attr->Value());
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
                        else
                        {
                            std::string errstr("Unknown attribute: ");
                            errstr += attrName;
                            ASSERTL1(false, errstr.c_str());
                        }

                        // Get the next attribute.
                        attr = attr->Next();
                    }

                    // Check to see if homogeneous expansion and if so
                    // strip down the shapeString definition
                    int numHomoDir = 0;
                    size_t loc;
                    //---> This finds the first location of  'n'!
                    if((loc = shapeString.find_first_of("-"))!=string::npos)
                    {
                        if(shapeString.find("Exp1D")!=string::npos)
                        {
                            numHomoDir = 1;
                        }
                        else // HomogeneousExp1D
                        {
                            numHomoDir = 2;
                        }

                        shapeString.erase(loc,shapeString.length());
                    }

                    // Reconstruct the fielddefs.
                    std::vector<unsigned int> elementIds;
                    {
                        bool valid = ParseUtils::GenerateSeqVector(idString.c_str(), elementIds);
                        ASSERTL0(valid, "Unable to correctly parse the element ids.");
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

                    ASSERTL0(valid, std::string("Unable to correctly parse the shape type: ").append(shapeString).c_str());

                    // Get the basis
                    std::vector<std::string> basisStrings;
                    std::vector<BasisType> basis;
                    valid = ParseUtils::GenerateOrderedStringVector(basisString.c_str(), basisStrings);
                    ASSERTL0(valid, "Unable to correctly parse the basis types.");
                    for (std::vector<std::string>::size_type i = 0; i < basisStrings.size(); i++)
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
                        ASSERTL0(valid, std::string("Unable to correctly parse the basis type: ").append(basisStrings[i]).c_str());
                    }

                    // Get homoLengths
                    std::vector<NekDouble> homoLengths;
                    if(numHomoDir)
                    {
                        valid = ParseUtils::GenerateUnOrderedVector(homoLengthsString.c_str(), homoLengths);
                        ASSERTL0(valid, "Unable to correctly parse the number of homogeneous lengths.");
                    }
					
					// Get Homogeneous points IDs
					std::vector<unsigned int> homoZIDs;
					std::vector<unsigned int> homoYIDs;
					
					if(numHomoDir == 1)
                    {
                        valid = ParseUtils::GenerateSeqVector(homoZIDsString.c_str(), homoZIDs);
                        ASSERTL0(valid, "Unable to correctly parse homogeneous planes IDs.");
                    }
					
					if(numHomoDir == 2)
					{
						valid = ParseUtils::GenerateSeqVector(homoZIDsString.c_str(), homoZIDs);
                        ASSERTL0(valid, "Unable to correctly parse homogeneous lines IDs in z-direction.");
						valid = ParseUtils::GenerateSeqVector(homoYIDsString.c_str(), homoYIDs);
                        ASSERTL0(valid, "Unable to correctly parse homogeneous lines IDs in y-direction.");
					}
					

                    // Get points type
                    std::vector<PointsType> points;

                    if(pointDef)
                    {
                        std::vector<std::string> pointsStrings;
                        valid = ParseUtils::GenerateOrderedStringVector(pointsString.c_str(), pointsStrings);
                        ASSERTL0(valid, "Unable to correctly parse the points types.");
                        for (std::vector<std::string>::size_type i = 0; i < pointsStrings.size(); i++)
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

                            ASSERTL0(valid, std::string("Unable to correctly parse the points type: ").append(pointsStrings[i]).c_str());
                        }
                    }

                    // Get numModes
                    std::vector<unsigned int> numModes;
                    bool UniOrder = false;

                    if(strstr(numModesString.c_str(),"UNIORDER:"))
                    {
                        UniOrder  = true;
                    }

                    valid = ParseUtils::GenerateOrderedVector(numModesString.c_str()+9, numModes);
                    ASSERTL0(valid, "Unable to correctly parse the number of modes.");

                    // Get numPoints
                    std::vector<unsigned int> numPoints;
                    if(numPointDef)
                    {
                        valid = ParseUtils::GenerateOrderedVector(numPointsString.c_str(), numPoints);
                        ASSERTL0(valid, "Unable to correctly parse the number of points.");
                    }

                    // Get fields names
                    std::vector<std::string> Fields;
                    valid = ParseUtils::GenerateOrderedStringVector(fieldsString.c_str(), Fields);
                    ASSERTL0(valid, "Unable to correctly parse the number of fields.");

                    FieldDefinitionsSharedPtr fielddef  = MemoryManager<FieldDefinitions>::AllocateSharedPtr(shape, elementIds, basis, UniOrder, numModes, Fields, numHomoDir, homoLengths, homoZIDs, homoYIDs, points, pointDef, numPoints, numPointDef);
                    
                    fielddefs.push_back(fielddef);

                    element = element->NextSiblingElement("ELEMENTS");
                }
                loopXml = loopXml->NextSiblingElement(strLoop);
            }
        }


        /**
         *
         */
        void ImportFieldData(TiXmlDocument &doc, const std::vector<FieldDefinitionsSharedPtr> &fielddefs, std::vector<std::vector<NekDouble> > &fielddata)
        {
            int cntdumps = 0;
            ASSERTL1(fielddata.size() == 0, "Expected an empty fielddata vector.");

            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");

            // Loop through all nektar tags, finding all of the element tags.
            while (master)
            {
                TiXmlElement* element = master->FirstChildElement("ELEMENTS");
                ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");
                while (element)
                {
                    // Extract the body, which the "data".
                    TiXmlNode* elementChild = element->FirstChild();
                    ASSERTL0(elementChild, "Unable to extract the data from the element tag.");
                    std::string elementStr;
                    while(elementChild)
                    {
                        if (elementChild->Type() == TiXmlNode::TEXT)
                        {
                            elementStr += elementChild->ToText()->ValueStr();
                        }
                        elementChild = elementChild->NextSibling();
                    }

                    // Convert from base64 to binary.
                    typedef boost::archive::iterators::transform_width<
                            boost::archive::iterators::binary_from_base64<
                            std::string::const_iterator>, 8, 6 > binary_t;
                    std::stringstream elementCompressedData(std::string(binary_t(elementStr.begin()), binary_t(elementStr.end())));

                    // Decompress the binary data.
                    std::stringstream elementDecompressedData;
                    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
                    in.push(boost::iostreams::zlib_decompressor());
                    in.push(elementCompressedData);
                    try
                    {
                        boost::iostreams::copy(in, elementDecompressedData);
                    }
                    catch (boost::iostreams::zlib_error e)
                    {
                        if (e.error() == boost::iostreams::zlib::stream_end)
                        {
                            ASSERTL0(false, "Stream end zlib error");
                        }
                        else if (e.error() == boost::iostreams::zlib::stream_error)
                        {
                            ASSERTL0(false, "Stream zlib error");
                        }
                        else if (e.error() == boost::iostreams::zlib::version_error)
                        {
                            ASSERTL0(false, "Version zlib error");
                        }
                        else if (e.error() == boost::iostreams::zlib::data_error)
                        {
                            ASSERTL0(false, "Data zlib error");
                        }
                        else if (e.error() == boost::iostreams::zlib::mem_error)
                        {
                            ASSERTL0(false, "Memory zlib error");
                        }
                        else if (e.error() == boost::iostreams::zlib::buf_error)
                        {
                            ASSERTL0(false, "Buffer zlib error");
                        }
                        else
                        {
                            ASSERTL0(false, "Unknown zlib error");
                        }
                    }

                    // Deserialize the array.
                    std::string vData = elementDecompressedData.str();
                    NekDouble* readFieldData = (NekDouble*) vData.c_str();
                    std::vector<NekDouble> elementFieldData(readFieldData, readFieldData + elementDecompressedData.str().length() * sizeof(*elementDecompressedData.str().c_str()) / sizeof(NekDouble));
                    fielddata.push_back(elementFieldData);

                    int datasize = CheckFieldDefinition(fielddefs[cntdumps]);
                    ASSERTL0(fielddata[cntdumps].size() == datasize*fielddefs[cntdumps]->m_fields.size(),"Input data is not the same length as header infoarmation");

                    cntdumps++;

                    element = element->NextSiblingElement("ELEMENTS");
                }
                master = master->NextSiblingElement("NEKTAR");
            }
        }


        /**
         *
         */
        int CheckFieldDefinition(const FieldDefinitionsSharedPtr &fielddefs)
        {
            int i;
            ASSERTL0(fielddefs->m_elementIDs.size() > 0, "Fielddefs vector must contain at least one element of data .");

            unsigned int numbasis = 0;

            // Determine nummodes vector lists are correct length
            switch(fielddefs->m_shapeType)
            {
            case eSegment:
                numbasis = 1;
                if(fielddefs->m_numHomogeneousDir)
                {
                    numbasis += fielddefs->m_numHomogeneousDir;
                }
                
                break;
            case eTriangle:  
            case eQuadrilateral:
                if(fielddefs->m_numHomogeneousDir)
                {
                    numbasis = 3;
                }
                else
                {
                    numbasis = 2;
                }
                break;
            case eTetrahedron:
            case ePyramid:
            case ePrism:
            case eHexahedron:
                numbasis = 3;
                break;
            default:
                ASSERTL0(false, "Unsupported shape type.");
                break;
            }
            
            unsigned int datasize = 0;
            
            ASSERTL0(fielddefs->m_basis.size() == numbasis, "Length of basis vector is incorrect");

            if(fielddefs->m_uniOrder == true)
            {
                unsigned int cnt = 0;
                // calculate datasize
                switch(fielddefs->m_shapeType)
                {
                case eSegment:
                {
                    int l = fielddefs->m_numModes[cnt++];
                    if(fielddefs->m_numHomogeneousDir == 1)
                    {
                        datasize += l*fielddefs->m_numModes[cnt++];
                    }
                    else if(fielddefs->m_numHomogeneousDir == 2)
                    {
                        int m = fielddefs->m_numModes[cnt++];
                        datasize += l*m*fielddefs->m_numModes[cnt++];
                    }
                    else
                    {
                        datasize += l;
                    }
                }
                break;
                case eTriangle:
                {
                    int l = fielddefs->m_numModes[cnt++];
                    int m = fielddefs->m_numModes[cnt++];
                    
                    if(fielddefs->m_numHomogeneousDir == 1)
                    {
                        datasize += StdTriData::getNumberOfCoefficients(l,m)*
                                    fielddefs->m_homogeneousZIDs.size();
                    }
                    else
                    {
                        datasize += StdTriData::getNumberOfCoefficients(l,m);
                    }
                }
                break;
                case eQuadrilateral:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        if(fielddefs->m_numHomogeneousDir == 1)
                        {
                            datasize += l*m*fielddefs->m_homogeneousZIDs.size();
                        }
                        else
                        {
                            datasize += l*m;
                        }
                    }
                    break;
                case eTetrahedron:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdTetData::getNumberOfCoefficients(l,m,n);
                    }
                    break;
                case ePyramid:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdPyrData::getNumberOfCoefficients(l,m,n);
                    }
                    break;
                case  ePrism:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += StdPrismData::getNumberOfCoefficients(l,m,n);
                    }
                    break;
                case eHexahedron:
                    {
                        int l = fielddefs->m_numModes[cnt++];
                        int m = fielddefs->m_numModes[cnt++];
                        int n = fielddefs->m_numModes[cnt++];
                        datasize += l*m*n;
                    }
                    break;
                default:
                    ASSERTL0(false, "Unsupported shape type.");
                    break;
                }
                
                datasize *= fielddefs->m_elementIDs.size();
            }
            else
            {
                unsigned int cnt = 0;
                // calculate data length
                for(i = 0; i < fielddefs->m_elementIDs.size(); ++i)
                {
                    switch(fielddefs->m_shapeType)
                    {
                    case eSegment:
                        datasize += fielddefs->m_numModes[cnt++];
                        break;
                    case eTriangle:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            datasize += StdTriData::getNumberOfCoefficients(l,m);
                        }
                        break;
                    case eQuadrilateral:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            datasize += l*m;
                        }
                        break;
                    case eTetrahedron:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdTetData::getNumberOfCoefficients(l,m,n);
                        }
                        break;
                    case ePyramid:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdPyrData::getNumberOfCoefficients(l,m,n);
                        }
                        break;
                    case  ePrism:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += StdPrismData::getNumberOfCoefficients(l,m,n);
                        }
                        break;
                    case eHexahedron:
                        {
                            int l = fielddefs->m_numModes[cnt++];
                            int m = fielddefs->m_numModes[cnt++];
                            int n = fielddefs->m_numModes[cnt++];
                            datasize += l*m*n;
                        }
                        break;
                    default:
                        ASSERTL0(false, "Unsupported shape type.");
                        break;
                    }
                }
            }
            
            return datasize;
        }
        

    }
}
