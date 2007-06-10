////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/BoundaryConditions.cpp,v $
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
                TiXmlElement *variableElement = variablesElement->FirstChildElement("V");

                // Sequential counter for the composite numbers.
                int nextVariableNumber = -1;

                while (variableElement)
                {
                    /// All elements are of the form: "<V ID="#"> name = value </V>", with
                    /// ? being the element type.

                    nextVariableNumber++;

                    TiXmlAttribute *attr = variableElement->FirstAttribute();
                    std::string attrName(attr->Name());

                    ASSERTL0(attrName == "ID", (std::string("Unknown attribute: ") + attrName).c_str());
                    int indx;
                    int err = attr->QueryIntValue(&indx);
                    ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
                    ASSERTL0(indx == nextVariableNumber, "Composite IDs must begin with zero and be sequential.");

                    TiXmlNode* variableChild = variableElement->FirstChild();
                    // This is primarily to skip comments that may be present.
                    // Comments appear as nodes just like elements.
                    // We are specifically looking for text in the body
                    // of the definition.
                    while(variableChild && variableChild->Type() != TiXmlNode::TEXT)
                    {
                        variableChild = variableChild->NextSibling();
                    }

                    ASSERTL0(variableChild, "Unable to read variable definition body.");
                    std::string variableName = variableChild->ToText()->ValueStr();

                    std::istringstream variableStrm(variableName);

                    variableStrm >> variableName;
                    m_Variables.push_back(variableName);

                    variableElement = variableElement->NextSiblingElement("V");
                }

                ASSERTL0(nextVariableNumber > 0, "Number of variables must be greater than zero.")
            }
        }

        void BoundaryConditions::ReadBoundaryRegions(TiXmlElement *conditions)
        {
            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");
            ASSERTL0(boundaryRegions, "Unable to find BOUNDARYREGIONS block.");

            int regionIndx = 0;

            // See if we have boundary regions defined.
            TiXmlElement *boundaryRegionsElement = boundaryRegions->FirstChildElement("B");

            // Sequential counter for the composite numbers.
            int nextBoundaryRegionNumber = -1;

            while (boundaryRegionsElement)
            {
                /// All elements are of the form: "<B ID="#"> ... </B>", with
                /// ? being the element type.

                nextBoundaryRegionNumber++;

                TiXmlAttribute *attr = boundaryRegionsElement->FirstAttribute();
                std::string attrName(attr->Name());

                ASSERTL0(attrName == "ID", (std::string("Unknown attribute: ") + attrName).c_str());
                int indx;
                int err = attr->QueryIntValue(&indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
                ASSERTL0(indx == nextBoundaryRegionNumber, "Boundary region IDs must begin with zero and be sequential.");

                TiXmlNode* boundaryRegionChild = boundaryRegionsElement->FirstChild();
                // This is primarily to skip comments that may be present.
                // Comments appear as nodes just like elements.
                // We are specifically looking for text in the body
                // of the definition.
                while(boundaryRegionChild && boundaryRegionChild->Type() != TiXmlNode::TEXT)
                {
                    boundaryRegionChild = boundaryRegionChild->NextSibling();
                }

                ASSERTL0(boundaryRegionChild, "Unable to read variable definition body.");
                std::string boundaryRegionStr = boundaryRegionChild->ToText()->ValueStr();

                std::string::size_type indxBeg = boundaryRegionStr.find_first_of('[') + 1;
                std::string::size_type indxEnd = boundaryRegionStr.find_last_of(']') - 1;

                ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading boundary region definition:") + boundaryRegionStr).c_str());

                std::string indxStr = boundaryRegionStr.substr(indxBeg, indxEnd - indxBeg + 1);
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
                        NEKERROR(ErrorUtil::efatal, (std::string("Cannot read composite definition: ") + boundaryRegionStr).c_str());
                    }
                }

                boundaryRegionsElement = boundaryRegionsElement->NextSiblingElement("B");
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
