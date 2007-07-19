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
#include <SpatialDomains/Equation.hpp>

#include <string>

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

                ASSERTL0(nextVariableNumber > -1, "Number of variables must be greater than zero.")
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
                typedef vector<unsigned int> SeqVector;
                SeqVector seqVector;

                if (!indxStr.empty())
                {
                    BoundaryRegionShPtr boundaryRegion(MemoryManager<BoundaryRegion>::AllocateSharedPtr());

                    if (ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector))
                    {
                        for (SeqVector::iterator iter = seqVector.begin(); iter != seqVector.end(); ++iter)
                        {
                            Composite composite = m_MeshGraph->GetComposite(*iter);
                            if (composite)
                            {
                                boundaryRegion->push_back(composite);
                            }
                            else
                            {
                                char str[64];
                                ::sprintf(str, "%d", *iter);
                                NEKERROR(ErrorUtil::ewarning, (std::string("Undefined composite: ") + str).c_str());

                            }
                        }

                        m_BoundaryRegions.push_back(boundaryRegion);
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
            // Read REGION tags
            TiXmlElement *boundaryConditionsElement = conditions->FirstChildElement("BOUNDARYCONDITIONS");
            ASSERTL0(boundaryConditionsElement, "Boundary conditions must be specified.");

            TiXmlElement *regionElement = boundaryConditionsElement->FirstChildElement("REGION");

            ASSERTL0(regionElement, "One or more boundary conditions must be specified.");

            // Read R (Robin), D (Dirichlet), N (Neumann) [What about Cauchy?] tags

            while (regionElement)
            {
                BoundaryConditionMapShPtr boundaryConditions = MemoryManager<BoundaryConditionMap>::AllocateSharedPtr();

                TiXmlAttribute *attr = regionElement->FirstAttribute();
                ASSERTL0(attr,
                    "The REF attribute must be present to specify the boundary region to which the condition applies.");
                std::string attrName(attr->Name());
                ASSERTL0(attrName == "REF",
                    "The REF attribute must be present to specify the boundary region to which the condition applies.");

                int boundaryRegionID;
                int err = attr->QueryIntValue(&boundaryRegionID);
                ASSERTL0(err == TIXML_SUCCESS, "Error reading boundary region reference.");

                // Find the boundary region corresponding to this ID.
                std::string boundaryRegionIDStr;
                std::ostringstream boundaryRegionIDStrm(boundaryRegionIDStr);
                boundaryRegionIDStrm << boundaryRegionID;
                ASSERTL0(boundaryRegionID < m_BoundaryRegions.size(),
                    (std::string("Boundary region ID not found: ") + boundaryRegionIDStr).c_str());

                //// Need to also make sure that we only specify a region ID once.  Since they
                //// must be specified in order we can just check to see if that index exists.
                //// It should be the next one in the container.
                //ASSERTL0(boundaryRegionID == m_BoundaryConditions.size(),
                //    (std::string("Next boundary condition must be ID: ") + boundaryRegionIDStr).c_str());

                // Here is the boundary region.
                // m_BoundaryRegions[boundaryRegionID];

                TiXmlElement *conditionElement = regionElement->FirstChildElement();

                while (conditionElement)
                {
                    // Check type.
                    std::string conditionType = conditionElement->Value();
                    std::string attrData;

                    // All have var specified, or else all variables are zero.
                    attr = conditionElement->FirstAttribute();

                    Variable::iterator iter;

                    if (attr)
                    {
                        attrName = attr->Name();

                        ASSERTL0(attrName == "VAR", "First attribute must be VAR.");

                        attrData = attr->Value();

                        if (!attrData.empty())
                        {
                            iter = std::find(m_Variables.begin(), m_Variables.end(), attrData);

                            ASSERTL0(iter != m_Variables.end(), (std::string("Cannot find variable: ") + attrData).c_str());
                        }
                    }

                    if (conditionType == "N")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Neumann and are set to zero.
                            for (Variable::iterator varIter = m_Variables.begin();
                                varIter != m_Variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr("00.0"));
                                (*boundaryConditions)[*varIter]  = neumannCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            // Read the VALUE attribute.  It is the next and only other attribute.
                            attr = attr->Next();

                            if (attr)
                            {
                                attrName = attr->Name();

                                ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                attrData = attr->Value();

                                ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(attrData));
                                (*boundaryConditions)[*iter]  = neumannCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr("0"));
                                (*boundaryConditions)[*iter]  = neumannCondition;
                            }
                        }
                    }
                    else if (conditionType == "D")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Dirichlet and are set to zero.
                            for (Variable::iterator varIter = m_Variables.begin();
                                varIter != m_Variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr("0"));
                                (*boundaryConditions)[*varIter] = dirichletCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            // Read the VALUE attribute.  It is the next and only other attribute.
                            attr = attr->Next();

                            if (attr)
                            {
                                attrName = attr->Name();

                                ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                attrData = attr->Value();

                                ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(attrData));
                                (*boundaryConditions)[*iter]  = dirichletCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr("0"));
                                (*boundaryConditions)[*iter]  = dirichletCondition;
                            }
                        }
                    }
                    else if (conditionType == "R")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Robin and are set to zero.
                            for (Variable::iterator varIter = m_Variables.begin();
                                varIter != m_Variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr("0", "0"));
                                (*boundaryConditions)[*varIter] = robinCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            // Read the A and B attributes.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string attrName1, attrName2;
                                std::string attrData1, attrData2;

                                attrName1 = attr->Name();

                                ASSERTL0(attrName1 == "A", (std::string("Unknown attribute: ") + attrName1).c_str());

                                attrData1 = attr->Value();

                                ASSERTL0(!attrData1.empty(), "A attributes must have associated values.");

                                attr = attr->Next();
                                ASSERTL0(attr, "Unable to read B attribute.");

                                attrName2 = attr->Name();
                                ASSERTL0(attrName2 == "B", (std::string("Unknown attribute: ") + attrName2).c_str());

                                attrData2 = attr->Value();

                                ASSERTL0(!attrData2.empty(), "B attributes must have associated values.");

                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(attrData1, attrData2));
                                (*boundaryConditions)[*iter]  = robinCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr("0", "0"));
                                (*boundaryConditions)[*iter]  = robinCondition;
                            }
                        }
                    }
                    else if (conditionType == "C")
                    {
                        NEKERROR(ErrorUtil::ewarning, "Cauchy type boundary conditions not implemented.");
                    }

                    conditionElement = conditionElement->NextSiblingElement();
                }

                m_BoundaryConditions[boundaryRegionID] = boundaryConditions;
                regionElement = regionElement->NextSiblingElement("REGION");
            }
       }

        void BoundaryConditions::ReadForcingFunctions(TiXmlElement *conditions)
        {
            TiXmlElement *forcingFunctionsElement = conditions->FirstChildElement("FORCING");

            if (forcingFunctionsElement)
            {
                TiXmlElement *forcingFunction = forcingFunctionsElement->FirstChildElement("F");

                // All forcing functions are initialized to "0" so they only have to be
                // partially specified.  That is, not all variables have to have functions
                // specified.  For those that are missing it is assumed they are "0".
                for (Variable::iterator varIter = m_Variables.begin();
                    varIter != m_Variables.end(); ++varIter)
                {
                    ForcingFunctionShPtr forcingFunctionShPtr(MemoryManager<ForcingFunction>::AllocateSharedPtr("0"));
                    m_ForcingFunctions[*varIter] = forcingFunctionShPtr;
                }

                while (forcingFunction)
                {
                    TiXmlAttribute *variableAttr = forcingFunction->FirstAttribute();

                    ASSERTL0(variableAttr, "The variable must be specified for the forcing function.");
                    std::string variableAttrName = variableAttr->Name();
                    ASSERTL0(variableAttrName == "VAR", (std::string("Error in forcing function attribute name: ") + variableAttrName).c_str());

                    std::string variableStr = variableAttr->Value();

                    TiXmlAttribute *functionAttr = variableAttr->Next();
                    if (functionAttr)
                    {
                        ForcingFunctionsMap::iterator forcingFcnsIter = m_ForcingFunctions.find(variableStr);

                        if (forcingFcnsIter != m_ForcingFunctions.end())
                        {
                            m_ForcingFunctions[variableStr]->SetEquation(functionAttr->Value());
                        }
                        else
                        {
                            NEKERROR(ErrorUtil::efatal, (std::string("Error setting forcing function for variable: ") + variableStr).c_str());
                        }
                    }

                    forcingFunction = forcingFunction->NextSiblingElement("F");
                }
            }
        }

        void BoundaryConditions::ReadInitialConditions(TiXmlElement *conditions)
        {
            TiXmlElement *initialConditionsElement = conditions->FirstChildElement("INITIALCONDITIONS");
            ASSERTL0(initialConditionsElement, "Initial conditions must be specified.");

            TiXmlElement *initialCondition = initialConditionsElement->FirstChildElement("I");

            // All initial conditions are initialized to "0" so they only have to be
            // partially specified.  That is, not all variables have to have functions
            // specified.  For those that are missing it is assumed they are "0".
            for (Variable::iterator varIter = m_Variables.begin();
                varIter != m_Variables.end(); ++varIter)
            {
                InitialConditionShPtr initialConditionShPtr(MemoryManager<InitialCondition>::AllocateSharedPtr("00.0"));
                m_InitialConditions[*varIter] = initialConditionShPtr;
            }

            while (initialCondition)
            {
                TiXmlAttribute *initialConditionAttr = initialCondition->FirstAttribute();

                ASSERTL0(initialConditionAttr, "The variable must be specified for the forcing function.");
                std::string intialConditionAttrName = initialConditionAttr->Name();
                ASSERTL0(intialConditionAttrName == "VAR", (std::string("Error in initial condition attribute name: ") + intialConditionAttrName).c_str());

                std::string initialConditionStr = initialConditionAttr->Value();

                TiXmlAttribute *functionAttr = initialConditionAttr->Next();
                if (functionAttr)
                {
                    InitialConditionsMap::iterator initialConditionFcnsIter =
                        m_InitialConditions.find(initialConditionStr);

                    if (initialConditionFcnsIter != m_InitialConditions.end())
                    {
                        m_InitialConditions[initialConditionStr]->SetEquation(functionAttr->Value());
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Error setting initial condition for variable: ") + initialConditionStr).c_str());
                    }
                }

                initialCondition = initialCondition->NextSiblingElement("I");
            }
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(int indx)
        {
            ConstForcingFunctionShPtr returnval;

            if (indx >= m_Variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index: ") + errStr).c_str());
            }

            return GetForcingFunction(m_Variables[indx]);
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(const string &var)
        {
            ConstForcingFunctionShPtr returnval;
            
            // Check that var is defined in forcing function list.
            ForcingFunctionsMap::iterator ffIter = m_ForcingFunctions.find(var);

            bool found = false;
            if( ffIter != m_ForcingFunctions.end() )
            {
                returnval = ffIter->second;
                found = true;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable used in obtaining forcing function: ") + var).c_str());
            }

            if (found)
            {
                if (returnval->GetEquation() == "00.0")
                {
                    NEKERROR(ErrorUtil::efatal,
                        (std::string("Default forcing function used for variable: ") + var).c_str());
                }
            }

            return returnval;
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(const char *var)
        {
            return GetForcingFunction(std::string(var));
        }

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(int indx)
        {
            ConstInitialConditionShPtr returnval;

            if (indx >= m_Variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index: ") + errStr).c_str());
            }

            return GetInitialCondition(m_Variables[indx]);
        }

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(const string &var)
        {
            ConstInitialConditionShPtr returnval;
            
            // Check that var is defined in forcing function list.
            InitialConditionsMap::iterator ffIter = m_InitialConditions.find(var);

            bool found = false;
            if( ffIter != m_InitialConditions.end() )
            {
                returnval = ffIter->second;
                found = true;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable used in obtaining initial condition: ") + var).c_str());
            }

            if (found)
            {
                if (returnval->GetEquation() == "00.0")
                {
                    NEKERROR(ErrorUtil::efatal,
                        (std::string("Default initial condition used for variable: ") + var).c_str());
                }
            }

            return returnval;
        }

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(const char *var)
        {
            return GetInitialCondition(std::string(var));
        }
    }
}
