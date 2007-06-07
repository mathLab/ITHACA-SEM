////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source$
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/pchSpatialDomains.h>

#include <SpatialDomains/BoundaryConditions.h>
#include <SpatialDomains/ParseUtils.hpp>

namespace Nektar
{
    namespace SpatialDomains
    {

        BoundaryConditions::BoundaryConditions(const MeshGraph *meshGraph):
            m_MeshGraph(meshGraph)
        {
        }

        BoundaryConditions::~BoundaryConditions()
        {
        }

        void BoundaryConditions::Read(std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") + 
                infilename).c_str());

            Read(doc);
        }

        // \brief Read segments (and general MeshGraph) given TiXmlDocument.
        void BoundaryConditions::Read(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            TiXmlNode* node = NULL;
            TiXmlElement* conditions = NULL;

            /// Look for all data in CONDITIONS block.
            conditions = docHandle.FirstChildElement("NEKTAR").FirstChildElement("CONDITIONS").Element();

            ASSERTL0(conditions, "Unable to find BOUNDARYCONDITIONS tag in file.");

            // Now read all the different tagged sections
            ReadParameters(conditions);
            ReadVariables(conditions);
            ReadBoundaryRegions(conditions);
            ReadExpansionTypes(conditions);
            ReadBoundaryConditions(conditions);
            ReadForcingFunctions(conditions);
            ReadInitialConditions(conditions);
        }

        void BoundaryConditions::ReadParameters(TiXmlElement *conditions)
        {
            TiXmlElement *parametersElement = conditions->FirstChildElement("PARAMETERS");

            // See if we have parameters defined.  They are optional so we go on if not.
            if (parametersElement)
            {
                TiXmlElement *parameter = parametersElement->FirstChildElement("P");

                // Multiple nodes will only occur if there is a comment in between
                // definitions.
                while (parameter)
                {
                    TiXmlNode *node = parameter->FirstChild();

                    while (node && node->Type() != TiXmlNode::TEXT)
                    {
                        node = node->NextSibling();
                    }

                    if (node)
                    {
                        // Format is "paramName = value"
                        std::string line = node->ToText()->Value();
                        std::string symbol;
                        double value;

                        if (ParseUtils::ParseRealAssignment(line.c_str(), symbol, value))
                        {
                            m_Parameters[symbol] = value;
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }
            }
        }

        void BoundaryConditions::ReadVariables(TiXmlElement *conditions)
        {
            TiXmlElement *variablesElement = conditions->FirstChildElement("VARIABLES");

            int varIndex = 0;   // Current index, should be zero-based.

            // See if we have parameters defined.  They are optional so we go on if not.
            if (variablesElement)
            {
                TiXmlNode *node = variablesElement->FirstChild();

                // Multiple nodes will only occur if there is a comment in between
                // definitions.
                while (node)
                {
                    int nodeType = node->Type();
                    if (nodeType == TiXmlNode::TEXT)
                    {
                        // Format is "number variable"
                        // Due to the way XML is processed, all variable
                        // definitions will be on the same line unless seperated
                        // by a comment or other node type (not allowed).
                        std::string line = node->ToText()->Value();

                        try
                        {
                            // Skip leading blanks.
                            line = line.substr(line.find_first_not_of(" "));

                            // Repeat reading "number variable" until string is depleted.
                            while (!line.empty())
                            {
                                // Find the end of the number.
                                std::string::size_type endPos = line.find_first_of(" ");

                                // Extract the number and remove it from the line
                                // including the trailing spaces.
                                int indx = atoi(line.substr(0, endPos).c_str());
                                ASSERTL0(indx == varIndex, "Index error.  Variable indices must be zero-based and sequential.");
                                ++varIndex;

                                line = line.substr(endPos);
                                line = line.substr(line.find_first_not_of(" "));

                                // Find the end of the variable, which may be the end
                                // of the line as well (find_first_of will return
                                // string::npos in this case).
                                endPos = line.find_first_of(" ");
                                m_Variables.push_back(line.substr(0,endPos));

                                // Need to go to the next definition, or if at the end
                                // clear the line so our loop will terminate.
                                if (endPos != std::string::npos)
                                {
                                    line = line.substr(endPos);
                                    line = line.substr(line.find_first_not_of(" "));
                                }
                                else
                                {
                                    line.clear();
                                }
                            }
                        }
                        catch(...)
                        {
                            NEKERROR(ErrorUtil::efatal, "Error processing variable definitions.");
                        }
                    }
                    else if (nodeType != TiXmlNode::COMMENT)
                    {
                        NEKERROR(ErrorUtil::ewarning, "Unknown node type in VARIABLES block.")
                    }

                    node = node->NextSibling();
                }

                ASSERTL0(varIndex > 0, "Number of variables must be greater than zero.")
            }
        }

        void BoundaryConditions::ReadBoundaryRegions(TiXmlElement *conditions)
        {
            TiXmlElement *boundaryRegionsElement = conditions->FirstChildElement("BOUNDARYREGIONS");
            ASSERTL0(boundaryRegionsElement, "Unable to find BOUNDARYREGIONS block.");

            int regionIndx = 0;

            // See if we have boundary regions defined.
            TiXmlNode *node = boundaryRegionsElement->FirstChild();

            // Multiple nodes will only occur if there is a comment in between
            // definitions.
            while (node)
            {
                int nodeType = node->Type();
                if (nodeType == TiXmlNode::TEXT)
                {
                    // Format is "number <B composites /B>"
                    // Due to the way XML is processed, all variable
                    // definitions will be on the same line unless seperated
                    // by a comment or other node type (not allowed).
                    std::string line = node->ToText()->Value();

                    try
                    {
                        // Skip leading blanks.
                        line = line.substr(line.find_first_not_of(" "));

                        int indx = atoi(line.c_str());
                        ASSERTL0(indx == regionIndx,
                            "Index error.  Boundary region indices must be zero-based and sequential.");
                        ++regionIndx;

                        // Now should read another element.
                        node = node->NextSibling();

                        // Skip any comments.
                        while (node->Type() == TiXmlNode::COMMENT)
                        {
                            node = node->NextSibling();
                        }

                        ASSERTL0(node->Type() == TiXmlNode::ELEMENT,
                            "Error processing boundary region definitions. Boundary region ill-formed.");

                        ASSERTL0(node->ValueStr() == "B", "Only 'B' tag is recognized in boundary region definition.");

                        TiXmlNode *child = node->FirstChild();

                        ASSERTL0(child, "No boundary region definition was found between tags.");

                        while (child)
                        {
                            // At this time process composites only.
                            int type = child->Type();
                            if (type == TiXmlNode::TEXT)
                            {
                                //...and finally we arrive at the composites defining the boundary region.
                                // All composites will be together unless another node is between them.
                                std::string line = child->ToText()->Value();

                                std::string::size_type indxBeg = line.find_first_of('[') + 1;
                                std::string::size_type indxEnd = line.find_last_of(']') - 1;

                                ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading boundary region definition:") + line).c_str());

                                std::string indxStr = line.substr(indxBeg, indxEnd - indxBeg + 1);
                                typedef vector<unsigned int> SeqVectorType;
                                SeqVectorType seqVector;

                                if (!indxStr.empty())
                                {
                                    vector<Composite> compVector;

                                    if (ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector))
                                    {
                                        for (SeqVectorType::iterator iter = seqVector.begin(); iter != seqVector.end(); ++iter)
                                        {
                                            Composite composite = m_MeshGraph->GetComposite(*iter);
                                            if (composite)
                                            {
                                                compVector.push_back(composite);
                                            }
                                            else
                                            {
                                                char str[64];
                                                ::sprintf(str, "%d", *iter);
                                                NEKERROR(ErrorUtil::ewarning, (std::string("Undefined composite: ") + str).c_str());

                                            }
                                        }

                                        m_BoundaryRegions.push_back(compVector);
                                    }
                                    else
                                    {
                                        NEKERROR(ErrorUtil::efatal, (std::string("Cannot read composite definition: ") + line).c_str());
                                    }
                                }
                            }
                            else if (type != TiXmlNode::COMMENT)
                            {
                                NEKERROR(ErrorUtil::efatal, "Error processing boundary region definitions.");
                            }

                            child = child->NextSibling();
                        }
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal, "Error processing boundary region definitions.");
                    }
                }
                else if (nodeType != TiXmlNode::COMMENT)
                {
                    NEKERROR(ErrorUtil::ewarning, "Unknown node type in BOUNDARYREGION block.")
                }

                node = node->NextSibling();
            }
        }

        void BoundaryConditions::ReadExpansionTypes(TiXmlElement *conditions)
        {
            NEKERROR(ErrorUtil::ewarning, "ExpansionTypes not currently implemented.");
        }

        void BoundaryConditions::ReadBoundaryConditions(TiXmlElement *conditions)
        {
        }

        void BoundaryConditions::ReadForcingFunctions(TiXmlElement *conditions)
        {
        }

        void BoundaryConditions::ReadInitialConditions(TiXmlElement *conditions)
        {
        }
    }
}
