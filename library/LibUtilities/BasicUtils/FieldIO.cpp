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

#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/asio/ip/host_name.hpp>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicConst/GitRevision.h>

#include "zlib.h"
#include <set>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

// Buffer size for zlib compression/decompression
#define CHUNK 16384

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace ptime = boost::posix_time;
namespace ip = boost::asio::ip;
namespace berrc = boost::system::errc;

namespace Nektar
{
    namespace LibUtilities
    {
        /**
         * This function allows for data to be written to an FLD file when a
         * session and/or communicator is not instantiated. Typically used in
         * utilities which do not take XML input and operate in serial only.
         */
        void Write(     const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> >   &fielddata,
                        const FieldMetaDataMap &fieldinfomap)
        {
#ifdef NEKTAR_USE_MPI
            int size;
            int init;
            MPI_Initialized(&init);

            // If MPI has been initialised we can check the number of processes
            // and, if > 1, tell the user he should not be running this
            // function in parallel. If it is not initialised, we do not
            // initialise it here, and assume the user knows what they are
            // doing.
            if (init)
            {
                MPI_Comm_size( MPI_COMM_WORLD, &size );
                ASSERTL0(size == 1,
                     "This static function is not available in parallel. Please"
                     "instantiate a FieldIO object for parallel use.");
            }
#endif
            CommSharedPtr c = GetCommFactory().CreateInstance("Serial", 0, 0);
            FieldIO f(c);
            f.Write(outFile, fielddefs, fielddata, fieldinfomap);
        }


        /**
         * This function allows for data to be imported from an FLD file when
         * a session and/or communicator is not instantiated. Typically used in
         * utilities which only operate in serial.
         */
        LIB_UTILITIES_EXPORT void Import(
                        const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata,
                        FieldMetaDataMap &fieldinfomap,
                        const Array<OneD, int> ElementiDs)
        {
#ifdef NEKTAR_USE_MPI
            int size;
            int init;
            MPI_Initialized(&init);

            // If MPI has been initialised we can check the number of processes
            // and, if > 1, tell the user he should not be running this
            // function in parallel. If it is not initialised, we do not
            // initialise it here, and assume the user knows what they are
            // doing.
            if (init)
            {
                MPI_Comm_size( MPI_COMM_WORLD, &size );
                ASSERTL0(size == 1,
                     "This static function is not available in parallel. Please"
                     "instantiate a FieldIO object for parallel use.");
            }
#endif
            CommSharedPtr c = GetCommFactory().CreateInstance("Serial", 0, 0);
            FieldIO f(c);
            f.Import(infilename, fielddefs, fielddata, fieldinfomap, ElementiDs);
        }


        /**
         *
         */
        FieldIO::FieldIO(
                LibUtilities::CommSharedPtr pComm)
            : m_comm(pComm)
        {
        }


        /**
         *
         */
        void FieldIO::Write(const std::string &outFile,
                   std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                   std::vector<std::vector<NekDouble> > &fielddata, 
                   const FieldMetaDataMap &fieldmetadatamap)
        {
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
            std::string filename = SetUpOutput(outFile, fielddefs, fieldmetadatamap);

            // Create the file (partition)
            TiXmlDocument doc;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            TiXmlElement * root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            AddInfoTag(root,fieldmetadatamap);

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
                    GenerateSeqString(fielddefs[f]->m_elementIDs,idString);
                }
                elemTag->SetAttribute("ID", idString);

                std::string compressedDataString;
                ASSERTL0(Z_OK == Deflate(fielddata[f], compressedDataString),
                        "Failed to compress field data.");

                // If the string length is not divisible by 3,
                // pad it. There is a bug in transform_width
                // that will make it reference past the end
                // and crash.
                switch (compressedDataString.length() % 3)
                {
                case 1:
                    compressedDataString += '\0';
                case 2:
                    compressedDataString += '\0';
                    break;
                }

                // Convert from binary to base64.
                typedef boost::archive::iterators::base64_from_binary<
                        boost::archive::iterators::transform_width<
                        std::string::const_iterator, 6, 8> > base64_t;
                std::string base64string(base64_t(compressedDataString.begin()),
                        base64_t(compressedDataString.end()));
                elemTag->LinkEndChild(new TiXmlText(base64string));

            }
            doc.SaveFile(filename);
        }


        /**
         *
         */
        void FieldIO::Import(const std::string& infilename,
                    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    std::vector<std::vector<NekDouble> > &fielddata,
                    FieldMetaDataMap &fieldmetadatamap,
                    const Array<OneD, int> ElementIDs)
        {

            std::string infile = infilename;

            fs::path pinfilename(infilename);            
            
            if(fs::is_directory(pinfilename)) // check to see that infile is a directory
            {
                fs::path infofile("Info.xml");
                fs::path fullpath = pinfilename / infofile; 
                infile = PortablePath(fullpath);

                std::vector<std::string> filenames; 
                std::vector<std::vector<unsigned int> > elementIDs_OnPartitions;
                
            
                ImportMultiFldFileIDs(infile,filenames, elementIDs_OnPartitions,
                                      fieldmetadatamap);
                
                if(ElementIDs == NullInt1DArray) //load all fields
                {
                    for(int i = 0; i < filenames.size(); ++i)
                    {
                        fs::path pfilename(filenames[i]);
                        fullpath = pinfilename / pfilename; 
                        string fname = PortablePath(fullpath); 

                        TiXmlDocument doc1(fname);
                        bool loadOkay1 = doc1.LoadFile();
                        
                        std::stringstream errstr;
                        errstr << "Unable to load file: " << fname << std::endl;
                        errstr << "Reason: " << doc1.ErrorDesc() << std::endl;
                        errstr << "Position: Line " << doc1.ErrorRow() << ", Column " << doc1.ErrorCol() << std::endl;
                        ASSERTL0(loadOkay1, errstr.str());
                        
                        ImportFieldDefs(doc1, fielddefs, false);
                        if(fielddata != NullVectorNekDoubleVector)
                        {
                            ImportFieldData(doc1, fielddefs, fielddata);
                        }
                    }
                    
                }
                else // only load relevant partitions
                {
                    int i,j;
                    map<int,vector<int> > FileIDs;
                    map<int,vector<int> >::iterator it;
                    set<int> LoadFile;

                    for(i = 0; i < elementIDs_OnPartitions.size(); ++i)
                    {
                        for(j = 0; j < elementIDs_OnPartitions[i].size(); ++j)
                        {
                            FileIDs[elementIDs_OnPartitions[i][j]].push_back(i);
                        }
                    }
                    
                    for(i = 0; i < ElementIDs.num_elements(); ++i)
                    {
                        it = FileIDs.find(ElementIDs[i]);
                        if (it != FileIDs.end())
                        {
                            for (j = 0; j < it->second.size(); ++j)
                            {
                                LoadFile.insert(it->second[j]);
                            }
                        }
                    }
                    
                    set<int>::iterator iter; 
                    for(iter = LoadFile.begin(); iter != LoadFile.end(); ++iter)
                    {
                        fs::path pfilename(filenames[*iter]);
                        fullpath = pinfilename / pfilename; 
                        string fname = PortablePath(fullpath); 
                        TiXmlDocument doc1(fname);
                        bool loadOkay1 = doc1.LoadFile();
                        
                        std::stringstream errstr;
                        errstr << "Unable to load file: " << fname << std::endl;
                        errstr << "Reason: " << doc1.ErrorDesc() << std::endl;
                        errstr << "Position: Line " << doc1.ErrorRow() << ", Column " << doc1.ErrorCol() << std::endl;
                        ASSERTL0(loadOkay1, errstr.str());
                        
                        ImportFieldDefs(doc1, fielddefs, false);
                        if(fielddata != NullVectorNekDoubleVector)
                        {
                            ImportFieldData(doc1, fielddefs, fielddata);
                        }
                    }
                }
            }
            else // serial format case 
            {
                
                TiXmlDocument doc(infile);
                bool loadOkay = doc.LoadFile();
                
                std::stringstream errstr;
                errstr << "Unable to load file: " << infile << std::endl;
                errstr << "Reason: " << doc.ErrorDesc() << std::endl;
                errstr << "Position: Line " << doc.ErrorRow() << ", Column " << 
                    doc.ErrorCol() << std::endl;
                ASSERTL0(loadOkay, errstr.str());
                
                ImportFieldMetaData(doc,fieldmetadatamap);
                ImportFieldDefs(doc, fielddefs, false);
                if(fielddata != NullVectorNekDoubleVector)
                {
                    ImportFieldData(doc, fielddefs, fielddata);
                }
            }
        }


        /**
         *
         */
        void FieldIO::WriteMultiFldFileIDs(const std::string &outFile,
                                  const std::vector<std::string> fileNames,
                                  std::vector<std::vector<unsigned int> > &elementList,
                                  const FieldMetaDataMap &fieldmetadatamap)
        {
            TiXmlDocument doc;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            ASSERTL0(fileNames.size() == elementList.size(),"Outfile names and list of elements ids does not match");

            TiXmlElement * root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            AddInfoTag(root,fieldmetadatamap);

            for (int t = 0; t < fileNames.size(); ++t)
            {
                if(elementList[t].size())
                {
                    TiXmlElement * elemIDs = new TiXmlElement("Partition");
                    root->LinkEndChild(elemIDs);
                    
                    elemIDs->SetAttribute("FileName",fileNames[t]);
                    
                    string IDstring;
                    
                    GenerateSeqString(elementList[t],IDstring);

                    elemIDs->LinkEndChild(new TiXmlText(IDstring));
                }
            }

            doc.SaveFile(outFile);
        }


        /** 
         *
         */
         void FieldIO::ImportMultiFldFileIDs(const std::string &inFile,
                                    std::vector<std::string> &fileNames,
                                    std::vector<std::vector<unsigned int> > &elementList,
                                    FieldMetaDataMap &fieldmetadatamap)
         {        
             TiXmlDocument doc(inFile);
             bool loadOkay = doc.LoadFile();
             
             
             std::stringstream errstr;
             errstr << "Unable to load file: " << inFile<< std::endl;
             errstr << "Reason: " << doc.ErrorDesc() << std::endl;
             errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
             ASSERTL0(loadOkay, errstr.str());
             
             // Handle on XML document
             TiXmlHandle docHandle(&doc);
             
             // Retrieve main NEKTAR tag - XML specification states one
             // top-level element tag per file.
             TiXmlElement* master = doc.FirstChildElement("NEKTAR");
             ASSERTL0(master, "Unable to find NEKTAR tag in file.");

             // Partition element tag name
             std::string strPartition = "Partition";

             // First attempt to get the first Partition element
             TiXmlElement* fldfileIDs = master->FirstChildElement(strPartition.c_str());
             if (!fldfileIDs)
             {
                 // If this files try previous name
                 strPartition = "MultipleFldFiles";
                 fldfileIDs = master->FirstChildElement("MultipleFldFiles");
             }
             ASSERTL0(fldfileIDs,
                      "Unable to find 'Partition' or 'MultipleFldFiles' tag "
                      "within nektar tag.");

             while (fldfileIDs)
             {
                // Read file name of partition file
                const char *attr = fldfileIDs->Attribute("FileName");
                ASSERTL0(attr, "'FileName' not provided as an attribute of '"
                                + strPartition + "' tag.");
                fileNames.push_back(std::string(attr));

                const char* elementIDs = fldfileIDs->GetText();
                ASSERTL0(elementIDs, "Element IDs not specified.");

                std::string elementIDsStr(elementIDs);

                std::vector<unsigned int> idvec;
                ParseUtils::GenerateSeqVector(elementIDsStr.c_str(),idvec);

                elementList.push_back(idvec);

                fldfileIDs = fldfileIDs->NextSiblingElement(strPartition.c_str());
             }
         }

        void FieldIO::ImportFieldMetaData(std::string filename,
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
        

        void FieldIO::ImportFieldMetaData(TiXmlDocument &doc,
                                FieldMetaDataMap &fieldmetadatamap)
        {
            
            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = 0;    // Master tag within which all data is contained.
            TiXmlElement* metadata = 0;
            
            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");
            std::string strLoop = "NEKTAR";

            // Retain original metadata structure for backwards compatibility
            // TODO: Remove old metadata format
            metadata = master->FirstChildElement("FIELDMETADATA");
            if(metadata)
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

                    fieldmetadatamap[paramString] = paramBodyStr;
                    param = param->NextSiblingElement("P");
                }
            }

            // New metadata format
            metadata = master->FirstChildElement("Metadata");
            if(metadata)
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
        void FieldIO::ImportFieldDefs(TiXmlDocument &doc, std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                             bool expChild)
        {
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
        void FieldIO::ImportFieldData(TiXmlDocument &doc, const std::vector<FieldDefinitionsSharedPtr> &fielddefs, std::vector<std::vector<NekDouble> > &fielddata)
        {
            int cntdumps = 0;

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
                    std::string vCompressed(binary_t(elementStr.begin()),
                                            binary_t(elementStr.end()));

                    std::vector<NekDouble> elementFieldData;
                    ASSERTL0(Z_OK == Inflate(vCompressed, elementFieldData),
                            "Failed to decompress field data.");

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
         * \brief add information about provenance and fieldmetadata
         */
        void FieldIO::AddInfoTag(TiXmlElement * root,
                             const FieldMetaDataMap &fieldmetadatamap)
        {
            FieldMetaDataMap ProvenanceMap;

            // Nektar++ release version from VERSION file
            ProvenanceMap["NektarVersion"] = string(NEKTAR_VERSION);

            // Date/time stamp
            ptime::time_facet *facet = new ptime::time_facet("%d-%b-%Y %H:%M:%S");
            std::stringstream wss;
            wss.imbue(locale(wss.getloc(), facet));
            wss << ptime::second_clock::local_time();
            ProvenanceMap["Timestamp"] = wss.str();

            // Hostname
            boost::system::error_code ec;
            ProvenanceMap["Hostname"] = ip::host_name(ec);

            // Git information
            // If built from a distributed package, do not include this
            if (NekConstants::kGitSha1 != "GITDIR-NOTFOUND")
            {
                ProvenanceMap["GitSHA1"]   = NekConstants::kGitSha1;
                ProvenanceMap["GitBranch"] = NekConstants::kGitBranch;
            }

            TiXmlElement * infoTag = new TiXmlElement("Metadata");
            root->LinkEndChild(infoTag);

            TiXmlElement * v;
            FieldMetaDataMap::const_iterator infoit;

            TiXmlElement * provTag = new TiXmlElement("Provenance");
            infoTag->LinkEndChild(provTag);
            for (infoit = ProvenanceMap.begin(); infoit != ProvenanceMap.end(); ++infoit)
            {
                v = new TiXmlElement( (infoit->first).c_str() );
                v->LinkEndChild(new TiXmlText((infoit->second).c_str()));
                provTag->LinkEndChild(v);
            }

            //---------------------------------------------
            // write field info section
            if(fieldmetadatamap != NullFieldMetaDataMap)
            {
                for(infoit = fieldmetadatamap.begin(); infoit != fieldmetadatamap.end(); ++infoit)
                {
                    v = new TiXmlElement( (infoit->first).c_str() );
                    v->LinkEndChild(new TiXmlText((infoit->second).c_str()));
                    infoTag->LinkEndChild(v);
                }
            }
        }


        /**
         *
         */
        void FieldIO::GenerateSeqString(const std::vector<unsigned int> &elmtids,
                                      std::string &idString)
        {
            std::stringstream idStringStream;
            bool setdash = true;
            unsigned int endval;

            idStringStream << elmtids[0];
            for (int i = 1; i < elmtids.size(); ++i)
            {
                if(elmtids[i] == elmtids[i-1]+1)
                {
                    if(setdash)
                    {
                        idStringStream << "-";
                        setdash = false;
                    }

                    if(i == elmtids.size()-1) // last element
                    {
                        idStringStream << elmtids[i];
                    }
                    else
                    {
                        endval = elmtids[i];
                    }
                }
                else
                {
                    if(setdash == false) // finish off previous dash sequence
                    {
                        idStringStream << endval;
                        setdash = true;
                    }


                    idStringStream << "," << elmtids[i];
                }
            }
            idString = idStringStream.str();
        }


        /**
         *
         */
        std::string FieldIO::SetUpOutput(const std::string outname,
                const std::vector<FieldDefinitionsSharedPtr>& fielddefs,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            ASSERTL0(!outname.empty(), "Empty path given to SetUpOutput()");

            int nprocs = m_comm->GetSize();
            int rank   = m_comm->GetRank();

            // Directory name if in parallel, regular filename if in serial
            fs::path specPath (outname);

            // Remove any existing file which is in the way
            if(m_comm->RemoveExistingFiles())
            {
                try
                {
                    fs::remove_all(specPath);
                }
                catch (fs::filesystem_error& e)
                {
                    ASSERTL0(e.code().value() == berrc::no_such_file_or_directory,
                             "Filesystem error: " + string(e.what()));
                }
            }

            // serial processing just add ending.
            if(nprocs == 1)
            {
                cout << "Writing: " << specPath << endl;
                return LibUtilities::PortablePath(specPath);
            }

            // Compute number of elements on this process and share with other
            // processes. Also construct list of elements on this process from
            // available vector of field definitions.
            std::vector<unsigned int> elmtnums(nprocs,0);
            std::vector<unsigned int> idlist;
            int i;
            for (i = 0; i < fielddefs.size(); ++i)
            {
                elmtnums[rank] += fielddefs[i]->m_elementIDs.size();
                idlist.insert(idlist.end(), fielddefs[i]->m_elementIDs.begin(),
                                            fielddefs[i]->m_elementIDs.end());
            }
            m_comm->AllReduce(elmtnums,LibUtilities::ReduceMax);

            // Create the destination directory
            try
            {
                fs::create_directory(specPath);
            }
            catch (fs::filesystem_error& e)
            {
                ASSERTL0(false, "Filesystem error: " + string(e.what()));
            }

            // Collate per-process element lists on root process to generate
            // the info file.
            if (rank == 0)
            {
                std::vector<std::vector<unsigned int> > ElementIDs(nprocs);

                // Populate the list of element ID lists from all processes
                ElementIDs[0] = idlist;
                for (i = 1; i < nprocs; ++i)
                {
                    std::vector<unsigned int> tmp(elmtnums[i]);
                    m_comm->Recv(i, tmp);
                    ElementIDs[i] = tmp;
                }

                // Set up output names
                std::vector<std::string> filenames;
                for(int i = 0; i < nprocs; ++i)
                {
                    boost::format pad("P%1$07d.fld");
                    pad % i;
                    filenames.push_back(pad.str());
                }

                // Write the Info.xml file
                string infofile = LibUtilities::PortablePath(
                                            specPath / fs::path("Info.xml"));

                cout << "Writing: " << specPath << endl;
                WriteMultiFldFileIDs(infofile, filenames, ElementIDs,
                                     fieldmetadatamap);
            }
            else
            {
                // Send this process's ID list to the root process
                m_comm->Send(0, idlist);
            }

            // Pad rank to 8char filenames, e.g. P0000000.fld
            boost::format pad("P%1$07d.fld");
            pad % m_comm->GetRank();

            // Generate full path name
            fs::path poutfile(pad.str());
            fs::path fulloutname = specPath / poutfile;

            // Return the full path to the partition for this process
            return LibUtilities::PortablePath(fulloutname);
        }


        /**
         * Compress a vector of NekDouble values into a string using zlib.
         */
        int FieldIO::Deflate(std::vector<NekDouble>& in,
                        string& out)
        {
            int ret;
            unsigned have;
            z_stream strm;
            unsigned char* input = (unsigned char*)(&in[0]);
            string buffer;
            buffer.resize(CHUNK);

            /* allocate deflate state */
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);

            ASSERTL0(ret == Z_OK, "Error initializing Zlib.");

            strm.avail_in = in.size() * sizeof(NekDouble) / sizeof(char);
            strm.next_in = input;

            // Deflate input until output buffer is no longer full.
            do {
                strm.avail_out = CHUNK;
                strm.next_out = (unsigned char*)(&buffer[0]);

                ret = deflate(&strm, Z_FINISH);

                // Deflate can return Z_OK, Z_STREAM_ERROR, Z_BUF_ERROR or
                // Z_STREAM_END. All, except Z_STREAM_ERROR are ok.
                ASSERTL0(ret != Z_STREAM_ERROR, "Zlib stream error");

                have = CHUNK - strm.avail_out;
                out += buffer.substr(0, have);

            } while (strm.avail_out == 0);

            // Check all input was processed.
            ASSERTL0(strm.avail_in == 0, "Not all input was used.");

            // Check stream is complete.
            ASSERTL0(ret == Z_STREAM_END, "Stream not finished");

            // Clean-up and return
            (void)deflateEnd(&strm);
            return Z_OK;
        }


        /**
         * Decompress a zlib-compressed string into a vector of NekDouble
         * values.
         */
        int FieldIO::Inflate(std::string& in,
                    std::vector<NekDouble>& out)
        {
            int ret;
            unsigned have;
            z_stream strm;
            string buffer;
            buffer.resize(CHUNK);
            string output;

            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.avail_in = 0;
            strm.next_in = Z_NULL;
            ret = inflateInit(&strm);
            ASSERTL0(ret == Z_OK, "Error initializing zlib decompression.");

            strm.avail_in = in.size();
            strm.next_in = (unsigned char*)(&in[0]);

            do {
                strm.avail_out = CHUNK;
                strm.next_out = (unsigned char*)(&buffer[0]);

                ret = inflate(&strm, Z_NO_FLUSH);

                ASSERTL0(ret != Z_STREAM_ERROR, "Stream error occured.");

                switch (ret) {
                    case Z_NEED_DICT:
                        ret = Z_DATA_ERROR;     /* and fall through */
                    case Z_DATA_ERROR:
                    case Z_MEM_ERROR:
                        (void)inflateEnd(&strm);
                        return ret;
                }

                have = CHUNK - strm.avail_out;
                output += buffer.substr(0, have);

            } while (strm.avail_out == 0);

            (void)inflateEnd(&strm);

            if (ret == Z_STREAM_END)
            {
                NekDouble* readFieldData = (NekDouble*) output.c_str();
                unsigned int len = output.size() * sizeof(*output.c_str())
                                                 / sizeof(NekDouble);
                out.assign( readFieldData, readFieldData + len);
                return Z_OK;
            }
            else
            {
                return Z_DATA_ERROR;
            }
        }

        /**
         *
         */
        int FieldIO::CheckFieldDefinition(const FieldDefinitionsSharedPtr &fielddefs)
        {
            int i;

            if(fielddefs->m_elementIDs.size() == 0) // empty partition
            {
                return 0;
            }
            //ASSERTL0(fielddefs->m_elementIDs.size() > 0, "Fielddefs vector must contain at least one element of data .");

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
