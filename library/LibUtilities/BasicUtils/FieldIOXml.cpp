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
        std::string FieldIOXml::className
            = GetFieldIOFactory().RegisterCreatorFunction(
                    "Xml",
                    FieldIOXml::create,
                    "XML-based output of field data.");

        FieldIOXml::FieldIOXml(LibUtilities::CommSharedPtr pComm) : FieldIO(pComm)
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
         * Compress a vector of NekDouble values into a string using zlib.
         */
        int FieldIOXml::Deflate(std::vector<NekDouble>& in,
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

    }
}
