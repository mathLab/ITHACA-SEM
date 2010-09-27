////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditions.cpp
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

#include <SpatialDomains/Equation.h>
#include <SpatialDomains/Conditions.h>
#include <SpatialDomains/ParseUtils.hpp>
#include <cctype>
#include <algorithm>
#include <string>

namespace Nektar
{
    namespace SpatialDomains
    {
        ParamMap BoundaryConditions::m_parameters;

        BoundaryConditions::BoundaryConditions(const MeshGraph *meshGraph):
            m_meshGraph(meshGraph)
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

            ASSERTL0(conditions, "Unable to find CONDITIONS tag in file.");

            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");

            // Now read all the different tagged sections
            ReadSolverInfo(conditions);

            ReadParameters(conditions);

            ReadFunctions(conditions);

            ReadVariables(conditions);

            if(boundaryRegions)
            {
                ReadBoundaryRegions(conditions);

                ReadBoundaryConditions(conditions);
            }

            ReadForcingFunctions(conditions);

            ReadInitialConditions(conditions);

            ReadExactSolution(conditions);

            ReadUserDefinedEqn(conditions);
        }

        void BoundaryConditions::ReadParameters(TiXmlElement *conditions)
        {
            TiXmlElement *parametersElement = conditions->FirstChildElement("PARAMETERS");

            // See if we have parameters defined.  They are optional so we go on if not.
            if (parametersElement)
            {
                TiXmlElement *parameter = parametersElement->FirstChildElement("P");
                LibUtilities::ExpressionEvaluator expEvaluator;

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

                        /// Pull out lhs and rhs and eliminate any spaces.
                        int beg=line.find_first_not_of(" ");
                        int end=line.find_first_of("=");
                        std::string lhs = line.substr(line.find_first_not_of(" "), end-beg-1);
                        lhs = lhs.substr(0, lhs.find_last_not_of(" ")+1);

                        std::string rhs = line.substr(line.find_last_of("=")+1);
                        rhs = rhs.substr(rhs.find_first_not_of(" "));
                        rhs = rhs.substr(0, rhs.find_last_not_of(" ")+1);

                        /// We want the list of parameters to have their RHS evaluated,
                        /// so we use the expression evaluator to do the dirty work.
                        if (!lhs.empty() && !rhs.empty())
                        {
                            NekDouble value=0.0;
                            expEvaluator.DefineFunction("", rhs);
                            value =  expEvaluator.Evaluate();
                            m_parameters[lhs] = value;
                            expEvaluator.SetParameter(lhs, value);
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }

                // Set ourselves up for evaluation later.
                Equation::SetConstParameters(m_parameters);
            }
        }

        void BoundaryConditions::ReadFunctions(TiXmlElement *conditions)
        {
            TiXmlElement *functionElement = conditions->FirstChildElement("FUNCTIONS");

            if (functionElement)
            {
                TiXmlElement *function = functionElement->FirstChildElement("F");

                while (function)
                {
                    std::string lhsString = function->Attribute("LHS");
                    ASSERTL0(!lhsString.empty(), "Unable to find LHS value.");

                    std::string fcnString = function->Attribute("VALUE");
                    ASSERTL0(!fcnString.empty(), "Unable to find function value.");

                    FunctionMap::iterator fcnIter = m_functions.find(lhsString);
                    ASSERTL0(fcnIter == m_functions.end(),
                        (std::string("Function value: ") + lhsString + std::string("already specified.")).c_str());

                    m_functions[lhsString] = fcnString;
                    function = function->NextSiblingElement("F");
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

                    int indx;
                    int err = variableElement->QueryIntAttribute("ID", &indx);
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
                    m_variables.push_back(variableName);

                    variableElement = variableElement->NextSiblingElement("V");
                }

                ASSERTL0(nextVariableNumber > -1, "Number of variables must be greater than zero.");
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

                int indx;
                int err = boundaryRegionsElement->QueryIntAttribute("ID", &indx);
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

                if (!indxStr.empty())
                {
                    // Extract the composites from the string and return them in a list.
                    BoundaryRegionShPtr boundaryRegion(MemoryManager<BoundaryRegion>::AllocateSharedPtr());
                    m_meshGraph->GetCompositeList(indxStr, *boundaryRegion);

                    m_boundaryRegions.push_back(boundaryRegion);
                }

                boundaryRegionsElement = boundaryRegionsElement->NextSiblingElement("B");
            }
        }

        void BoundaryConditions::ReadBoundaryConditions(TiXmlElement *conditions)
        {
            // Read REGION tags
            TiXmlElement *boundaryConditionsElement = conditions->FirstChildElement("BOUNDARYCONDITIONS");
            ASSERTL0(boundaryConditionsElement, "Boundary conditions must be specified.");

            TiXmlElement *regionElement = boundaryConditionsElement->FirstChildElement("REGION");

            ASSERTL0(regionElement, "One or more boundary conditions must be specified.");

            // Read R (Robin), D (Dirichlet), N (Neumann), P (Periodic) [What about Cauchy?] tags

            while (regionElement)
            {
                BoundaryConditionMapShPtr boundaryConditions = MemoryManager<BoundaryConditionMap>::AllocateSharedPtr();

                int boundaryRegionID;
                int err = regionElement->QueryIntAttribute("REF", &boundaryRegionID);
                ASSERTL0(err == TIXML_SUCCESS, "Error reading boundary region reference.");

                // Find the boundary region corresponding to this ID.
                std::string boundaryRegionIDStr;
                std::ostringstream boundaryRegionIDStrm(boundaryRegionIDStr);
                boundaryRegionIDStrm << boundaryRegionID;
                ASSERTL0(boundaryRegionID < m_boundaryRegions.size(),
                (std::string("Boundary region ID not found: ") + boundaryRegionIDStr).c_str());

                //// Need to also make sure that we only specify a region ID once.  Since they
                //// must be specified in order we can just check to see if that index exists.
                //// It should be the next one in the container.
                //ASSERTL0(boundaryRegionID == m_boundaryConditions.size(),
                //    (std::string("Next boundary condition must be ID: ") + boundaryRegionIDStr).c_str());

                // Here is the boundary region.
                // m_boundaryRegions[boundaryRegionID];

                TiXmlElement *conditionElement = regionElement->FirstChildElement();

                while (conditionElement)
                {
                    // Check type.
                    std::string conditionType = conditionElement->Value();
                    std::string attrData;

                    // All have var specified, or else all variables are zero.
                    TiXmlAttribute *attr = conditionElement->FirstAttribute();

                    Variable::iterator iter;
                    std::string attrName;

                    attrData = conditionElement->Attribute("VAR");

                    if (!attrData.empty())
                    {
                        iter = std::find(m_variables.begin(), m_variables.end(), attrData);
                        ASSERTL0(iter != m_variables.end(), (std::string("Cannot find variable: ") + attrData).c_str());
                    }

                    if (conditionType == "N")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Neumann and are set to zero.
                            for (Variable::iterator varIter = m_variables.begin();
                                varIter != m_variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr("00.0"));
                                (*boundaryConditions)[*varIter]  = neumannCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string equation, userDefined;

                                while(attr) {

                                    attrName = attr->Name();

                                    if (attrName=="USERDEFINEDTYPE") {

                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                        // Suppose to go here?
                                        SubstituteFunction(attrData);

                                        userDefined = attrData;
                                     }
                                     else if(attrName=="VALUE")
                                     {
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must be specified.");

                                        SubstituteFunction(attrData);

                                        equation = attrData;
                                      }
                                      attr = attr->Next();
                                    }

                                    BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(equation, userDefined));
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
                            for (Variable::iterator varIter = m_variables.begin();
                                varIter != m_variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr("0"));
                                (*boundaryConditions)[*varIter] = dirichletCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string equation, userDefined;

                                while(attr) {

                                   attrName = attr->Name();

                                    if (attrName=="USERDEFINEDTYPE") {

                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                        SubstituteFunction(attrData);

                                        userDefined = attrData;
                                    }
                                    else if(attrName=="VALUE")
                                    {
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                        SubstituteFunction(attrData);

                                        equation = attrData;
                                    }
                                   attr = attr->Next();
                                }

                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(equation, userDefined));
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
                    else if (conditionType == "R") // Read du/dn +  PRIMCOEFF u = VALUE
                    {
                        if (attrData.empty())
                        {
                            // All variables are Robin and are set to zero.
                            for (Variable::iterator varIter = m_variables.begin();
                                varIter != m_variables.end(); ++varIter)
                            {
                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr("0", "0"));
                                (*boundaryConditions)[*varIter] = robinCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must
                            // point to the variable.  Read the A and
                            // B attributes.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string attrName1, attrName2;
                                std::string attrData1, attrData2;
                                std::string equation1, equation2, userDefined;

                                while(attr){

                                attrName1 = attr->Name();

                                if (attrName1=="USERDEFINEDTYPE") {

                                    // Do stuff for the user defined attribute
                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                    SubstituteFunction(attrData1);

                                    userDefined = attrData1;

                                 } else if(attrName1 == "VALUE"){

                                    ASSERTL0(attrName1 == "VALUE", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "VALUE attributes must have associated values.");

                                    SubstituteFunction(attrData1);

                                    equation1 = attrData1;

                                    attr = attr->Next();
                                    ASSERTL0(attr, "Unable to read PRIMCOEFF attribute.");

                                    attrName1= attr->Name();
                                    ASSERTL0(attrName1 == "PRIMCOEFF", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "PRIMCOEFF attributes must have associated values.");

                                    SubstituteFunction(attrData1);

                                    equation2 = attrData1;

                                 }
                                 attr = attr->Next();

                                }

                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(equation1, equation2, userDefined));
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
                    else if (conditionType == "P")
                    {
                        if (attrData.empty())
                        {
                            attr = attr->Next();

                            if (attr)
                            {
                                attrName = attr->Name();

                                ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                int beg = attrData.find_first_of("[");
                                int end = attrData.find_first_of("]");
                                std::string periodicBndRegionIndexStr = attrData.substr(beg+1,end-beg-1);
                                ASSERTL0(beg < end, (std::string("Error reading periodic boundary region definition for boundary region: ")
                                                      + boundaryRegionIDStrm.str()).c_str());

                                vector<unsigned int> periodicBndRegionIndex;
                                bool parseGood = ParseUtils::GenerateSeqVector(periodicBndRegionIndexStr.c_str(), periodicBndRegionIndex);

                                ASSERTL0(parseGood && (periodicBndRegionIndex.size()==1), (std::string("Unable to read periodic boundary condition for boundary region: ")
                                                                              + boundaryRegionIDStrm.str()).c_str());

                                BoundaryConditionShPtr periodicCondition(MemoryManager<PeriodicBoundaryCondition>::AllocateSharedPtr(periodicBndRegionIndex[0]));

                                for (Variable::iterator varIter = m_variables.begin();
                                     varIter != m_variables.end(); ++varIter)
                                {
                                    (*boundaryConditions)[*varIter] = periodicCondition;
                                }
                            }
                            else
                            {
                                ASSERTL0(false, "Periodic boundary conditions should be explicitely defined");
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

                                int beg = attrData.find_first_of("[");
                                int end = attrData.find_first_of("]");
                                std::string periodicBndRegionIndexStr = attrData.substr(beg+1,end-beg-1);
                                ASSERTL0(beg < end, (std::string("Error reading periodic boundary region definition for boundary region: ")
                                                     + boundaryRegionIDStrm.str()).c_str());

                                vector<unsigned int> periodicBndRegionIndex;
                                bool parseGood = ParseUtils::GenerateSeqVector(periodicBndRegionIndexStr.c_str(), periodicBndRegionIndex);

                                ASSERTL0(parseGood && (periodicBndRegionIndex.size()==1), (std::string("Unable to read periodic boundary condition for boundary region: ") + boundaryRegionIDStrm.str()).c_str());

                                BoundaryConditionShPtr periodicCondition(MemoryManager<PeriodicBoundaryCondition>::AllocateSharedPtr(periodicBndRegionIndex[0]));
                                (*boundaryConditions)[*iter]  = periodicCondition;
                            }
                            else
                            {
                                ASSERTL0(false, "Periodic boundary conditions should be explicitely defined");
                            }
                        }
                    }
                    else if (conditionType == "C")
                    {
                        NEKERROR(ErrorUtil::ewarning, "Cauchy type boundary conditions not implemented.");
                    }

                    conditionElement = conditionElement->NextSiblingElement();
                }

                m_boundaryConditions[boundaryRegionID] = boundaryConditions;
                regionElement = regionElement->NextSiblingElement("REGION");
            }
       }


        void BoundaryConditions::ReadForcingFunctions(TiXmlElement *conditions)
        {
            TiXmlElement *forcingFunctionsElement = conditions->FirstChildElement("FORCING");

            if (forcingFunctionsElement)
            {
                TiXmlElement *forcingFunction = forcingFunctionsElement->FirstChildElement("F");

                // All forcing functions are initialized to "0" so
                // they only have to be partially specified.  That is,
                // not all variables have to have functions specified.
                // For those that are missing it is assumed they are
                // "0".
                for (Variable::iterator varIter = m_variables.begin();
                    varIter != m_variables.end(); ++varIter)
                {
                    ForcingFunctionShPtr forcingFunctionShPtr(MemoryManager<ForcingFunction>::AllocateSharedPtr("0"));
                    m_forcingFunctions[*varIter] = forcingFunctionShPtr;
                }

                while (forcingFunction)
                {
                    std::string variableStr = forcingFunction->Attribute("VAR");
                    ASSERTL0(!variableStr.empty(), "The variable must be specified for the forcing function.");

                    std::string fcnStr = forcingFunction->Attribute("VALUE");
                    ASSERTL0(!fcnStr.empty(),
                        (std::string("Forcing function for var: ") + variableStr + std::string(" must be specified.")).c_str());

                    /// Check the RHS against the functions defined in
                    /// m_functions.  If the name on the RHS is found
                    /// in the function map, then use the function
                    /// contained in the function map in place of the
                    /// RHS.
                    SubstituteFunction(fcnStr);

                    ForcingFunctionsMap::iterator forcingFcnsIter = m_forcingFunctions.find(variableStr);

                    if (forcingFcnsIter != m_forcingFunctions.end())
                    {
                        m_forcingFunctions[variableStr]->SetEquation(fcnStr);
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Error setting forcing function for variable: ") + variableStr).c_str());
                    }

                    forcingFunction = forcingFunction->NextSiblingElement("F");
                }
            }
        }

        void BoundaryConditions::ReadExactSolution(TiXmlElement *solution)
        {
            TiXmlElement *exactSolutionElement = solution->FirstChildElement("EXACTSOLUTION");

            // All exact solution functions are initialized to "0".
            // This means that even if the section is not specified
            // they will have a default value of "0". Alternatively
            // they only have to be partially specified.  That is, not
            // all variables have to have functions specified.  For
            // those that are missing it is assumed they are "0".
            for (Variable::iterator varIter = m_variables.begin();
                 varIter != m_variables.end(); ++varIter)
            {
                ExactSolutionShPtr exactSolutionShPtr(MemoryManager<ExactSolution>::AllocateSharedPtr("0"));
                m_exactSolution[*varIter] = exactSolutionShPtr;
            }

            if (exactSolutionElement)
            {
                TiXmlElement *exactSolution = exactSolutionElement->FirstChildElement("F");


                while (exactSolution)
                {
                    std::string variableStr = exactSolution->Attribute("VAR");
                    ASSERTL0(!variableStr.empty(), "The variable must be specified for the exact solution.");

                    std::string fcnStr = exactSolution->Attribute("VALUE");
                    ASSERTL0(!fcnStr.empty(), (std::string("The exact solution function must be specified for variable: ") + variableStr).c_str());

                    /// Check the RHS against the functions defined in
                    /// m_functions.  If the name on the RHS is found
                    /// in the function map, then use the function
                    /// contained in the function map in place of the
                    /// RHS.
                    SubstituteFunction(fcnStr);

                    ExactSolutionMap::iterator exactSolutionIter = m_exactSolution.find(variableStr);

                    if (exactSolutionIter != m_exactSolution.end())
                    {
                        m_exactSolution[variableStr]->SetEquation(fcnStr);
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Error setting exact solution for variable: ") + variableStr).c_str());
                    }

                    exactSolution = exactSolution->NextSiblingElement("F");
                }
            }
        }

        void BoundaryConditions::ReadSolverInfo(TiXmlElement *function)
        {
            TiXmlElement *solverInfoElement = function->FirstChildElement("SOLVERINFO");

            if (solverInfoElement)
            {
                TiXmlElement *solverInfo = solverInfoElement->FirstChildElement("I");

                while (solverInfo)
                {
                    std::string solverProperty = solverInfo->Attribute("PROPERTY");
                    // make sure that solver property is capitalised
                    transform(solverProperty.begin(), solverProperty.end(), solverProperty.begin(), (int(*)(int))std::toupper);
                    ASSERTL0(!solverProperty.empty(), "Unable to find PROPERTY value.");

                    std::string solverValue    = solverInfo->Attribute("VALUE");
                    ASSERTL0(!solverValue.empty(),"Unable to find VALUE string");

                    SolverInfoMap::iterator solverInfoIter = m_solverInfo.find(solverProperty);

                    ASSERTL0(solverInfoIter == m_solverInfo.end(),
                             (std::string("SolverInfo value: ") + solverProperty
                              + std::string(" already specified.")).c_str());

                    // Set Variable
                    m_solverInfo[solverProperty] = solverValue;
                    solverInfo = solverInfo->NextSiblingElement("I");
                }
            }
        }


        void BoundaryConditions::ReadUserDefinedEqn(TiXmlElement *function)
        {
            TiXmlElement *userDefinedEqnElement = function->FirstChildElement("USERDEFINEDEQNS");

            if (userDefinedEqnElement)
            {
                TiXmlElement *userDefinedEqn = userDefinedEqnElement->FirstChildElement("F");

                while (userDefinedEqn)
                {
                    std::string lhsString = userDefinedEqn->Attribute("LHS");
                    ASSERTL0(!lhsString.empty(), "Unable to find LHS value.");

                    std::string fcnString = userDefinedEqn->Attribute("VALUE");
                    ASSERTL0(!fcnString.empty(),"Unable to find eqn value");

                    /// Check the RHS against the functions defined in
                    /// m_functions.  If the name on the RHS is found
                    /// in the function map, then use the function
                    /// contained in the function map in place of the
                    /// RHS.
                    SubstituteFunction(fcnString);

                    UserDefinedEqnMap::iterator userDefinedEqnIter = m_userDefinedEqn.find(lhsString);

                    ASSERTL0(userDefinedEqnIter == m_userDefinedEqn.end(),
                             (std::string("UserDefinedEqn value: ") + lhsString
                              + std::string(" already specified.")).c_str());

                    // Set Variable
                    UserDefinedEqnShPtr userDefinedEqnShPtr(MemoryManager<UserDefinedEqn>::AllocateSharedPtr(fcnString));
                    m_userDefinedEqn[lhsString] = userDefinedEqnShPtr;
                    userDefinedEqn = userDefinedEqn->NextSiblingElement("F");
                }
            }
        }


        void BoundaryConditions::ReadInitialConditions(TiXmlElement *conditions)
        {
            TiXmlElement *initialConditionsElement = conditions->FirstChildElement("INITIALCONDITIONS");

            if (initialConditionsElement)
            {
                TiXmlElement *initialCondition = initialConditionsElement->FirstChildElement();

                // All initial conditions are initialized to "0" so
                // they only have to be partially specified.  That is,
                // not all variables have to have functions specified.
                // For those that are missing it is assumed they are
                // "0".
                for (Variable::iterator varIter = m_variables.begin();
                    varIter != m_variables.end(); ++varIter)
                {
                    InitialConditionShPtr initialConditionShPtr(MemoryManager<InitialCondition>::AllocateSharedPtr("00.0"));
                    m_initialConditions[*varIter] = initialConditionShPtr;
                }

                while (initialCondition)
                {
                    std::string conditionType = initialCondition->Value();

                    if(conditionType == "F")
                    {
                        std::string variableStr = initialCondition->Attribute("VAR");
                        ASSERTL0(!variableStr.empty(), "The variable must be specified for the initial solution.");

                        std::string fcnStr = initialCondition->Attribute("VALUE");
                        ASSERTL0(!fcnStr.empty(),
                                 (std::string("The initial condition function must be specified for variable: ") + variableStr).c_str());

                        /// Check the RHS against the functions defined in
                        /// m_functions.  If the name on the RHS is found
                        /// in the function map, then use the function
                        /// contained in the function map in place of the
                        /// RHS.
                        SubstituteFunction(fcnStr);

                        InitialConditionsMap::iterator initialConditionFcnsIter = m_initialConditions.find(variableStr);

                        if (initialConditionFcnsIter != m_initialConditions.end())
                        {
                            m_initialConditions[variableStr]->SetEquation(fcnStr);
                        }
                        else
                        {
                            NEKERROR(ErrorUtil::efatal,
                                     (std::string("Error setting initial condition for variable: ") + variableStr).c_str());
                        }
                    }
                    else if(conditionType == "R")
                    {
                        std::string restartStr = "RESTART";
                        std::string fileStr = initialCondition->Attribute("FILE");
                        ASSERTL0(!fileStr.empty(), "A file  must be specified for the restart initial solution.");

                        InitialConditionShPtr initialConditionShPtr(MemoryManager<InitialCondition>::AllocateSharedPtr("00.0"));

                        m_initialConditions[restartStr] = initialConditionShPtr;
                        m_initialConditions[restartStr]->SetEquation(fileStr);
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::ewarning,(std::string("Identifier ")+conditionType+std::string(" in Initial Conditions not recognised")).c_str());
                    }

                    initialCondition = initialCondition->NextSiblingElement();
                }
            }
        }

        bool BoundaryConditions::CheckForParameter(const std::string &parmName)
        {
            bool returnval;

            ParamMap::iterator paramMapIter = m_parameters.find(parmName);

            if(paramMapIter == m_parameters.end())
            {
                returnval = false;
            }
            else
            {
                returnval = true;
            }

            return returnval;
        }

        NekDouble BoundaryConditions::GetParameter(const std::string &parmName)
        {
            ParamMap::iterator paramMapIter = m_parameters.find(parmName);

            ASSERTL0(paramMapIter != m_parameters.end(),
                (std::string("Unable to find requested parameter: ") + parmName).c_str());

            return paramMapIter->second;
        }


        const std::string &BoundaryConditions::GetFunction(const std::string &lhs)
        {
            FunctionMap::iterator fcnIter = m_functions.find(lhs);

            ASSERTL1(fcnIter != m_functions.end(),
                (std::string("Unable to find requested function: ") + lhs).c_str());

            return fcnIter->second;
        }

        const std::string &BoundaryConditions::GetSolverInfo(const std::string &property)
        {
            SolverInfoMap::iterator slvIter = m_solverInfo.find(property);

            ASSERTL1(slvIter != m_solverInfo.end(),
                (std::string("Unable to find requested property: ") + property).c_str());

            return slvIter->second;
        }


        bool BoundaryConditions::SolverInfoExists(const std::string &property)
        {

            bool returnval = true;
            SolverInfoMap::iterator slvIter = m_solverInfo.find(property);

            if(slvIter == m_solverInfo.end())
            {
                returnval = false;
            }

            return returnval;
        }


        Equation BoundaryConditions::GetFunctionAsEquation(const std::string &lhs)
        {
            FunctionMap::iterator fcnIter = m_functions.find(lhs);

            ASSERTL1(fcnIter != m_functions.end(),
                (std::string("Unable to find requested function: ") + lhs).c_str());

            return Equation(fcnIter->second);
        }


        bool BoundaryConditions::SubstituteFunction(std::string &str)
        {
            FunctionMap::iterator fcnIter = m_functions.find(str);
            bool returnval = false;

            if (fcnIter != m_functions.end())
            {
                str = fcnIter->second;
                returnval = true;
            }

            return returnval;
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(int indx) const
        {
            ConstForcingFunctionShPtr returnval;

            if (indx >= m_variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index (ForcingFunction): ") + errStr).c_str());
            }

            return GetForcingFunction(m_variables[indx]);
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(const std::string &var) const
        {
            ConstForcingFunctionShPtr returnval;

            // Check that var is defined in forcing function list.
            ForcingFunctionsMap::const_iterator ffIter = m_forcingFunctions.find(var);

            bool found = false;
            if( ffIter != m_forcingFunctions.end() )
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

        bool BoundaryConditions::ExactSolutionExists(int indx) const
        {
            bool returnval = false;

            // Check that var is defined in forcing function list.
            ExactSolutionMap::const_iterator exSolnIter = m_exactSolution.find(m_variables[indx]);

            if(exSolnIter != m_exactSolution.end())
            {
                returnval = true;
            }

             return returnval;
        }

        ConstExactSolutionShPtr BoundaryConditions::GetExactSolution(int indx) const
        {
            ConstExactSolutionShPtr returnval;

            if (indx >= m_variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index (ExactSolution): ") + errStr).c_str());
            }

            return GetExactSolution(m_variables[indx]);
        }

        ConstExactSolutionShPtr BoundaryConditions::GetExactSolution(const std::string &var) const
        {
            ConstExactSolutionShPtr returnval;

            // Check that var is defined in forcing function list.
            ExactSolutionMap::const_iterator exSolnIter = m_exactSolution.find(var);

            if(exSolnIter != m_exactSolution.end())
            {
                returnval = exSolnIter->second;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable used in obtaining exact solution: ") + var).c_str());
            }

            return returnval;
        }


        bool BoundaryConditions::UserDefinedEqnExists(const std::string &var) const
        {
            bool returnval = false;

            // Check that var is defined in forcing function list.
            UserDefinedEqnMap::const_iterator userDefIter = m_userDefinedEqn.find(var);

            if(userDefIter != m_userDefinedEqn.end())
            {
                returnval = true;
            }

             return returnval;
        }


        ConstUserDefinedEqnShPtr BoundaryConditions::GetUserDefinedEqn(int indx) const
        {
            ConstUserDefinedEqnShPtr returnval;

            if (indx >= m_variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index (UserDefineEqn): ") + errStr).c_str());
            }

            return GetUserDefinedEqn(m_variables[indx]);
        }

        ConstUserDefinedEqnShPtr BoundaryConditions::GetUserDefinedEqn(const std::string &var) const
        {
            ConstUserDefinedEqnShPtr returnval;

            // Check that var is defined in forcing function list.
            UserDefinedEqnMap::const_iterator userDefIter = m_userDefinedEqn.find(var);

            if(userDefIter != m_userDefinedEqn.end())
            {
                returnval = userDefIter->second;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unable to find user defined equation LHS variable: ") + var).c_str());
            }

            return returnval;
        }


        bool BoundaryConditions::InitialConditionExists(int indx) const
        {
            bool returnval = false;

            // Check that var is defined in forcing function list.
            InitialConditionsMap::const_iterator ffIter = m_initialConditions.find(m_variables[indx]);

            bool found = false;
            if( ffIter != m_initialConditions.end() )
            {
                returnval = true;
            }

             return returnval;
        }

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(int indx) const
        {
            ConstInitialConditionShPtr returnval;

            if (indx >= m_variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index (InitialCondition): ") + errStr).c_str());
            }

            return GetInitialCondition(m_variables[indx]);
        }

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(const string &var) const
        {
            ConstInitialConditionShPtr returnval;

            // Check that var is defined in forcing function list.
            InitialConditionsMap::const_iterator ffIter = m_initialConditions.find(var);

            bool found = false;
            if( ffIter != m_initialConditions.end() )
            {
                returnval = ffIter->second;
                found = true;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable used in obtaining initial condition: ") + var).c_str());
            }

            return returnval;
        }


        bool BoundaryConditions::FoundInitialCondition(const std::string &var)
        {
            bool returnval = false;

            // Check that var is defined in forcing function list.
            InitialConditionsMap::const_iterator ffIter = m_initialConditions.find(var);

            if( ffIter != m_initialConditions.end() )
            {
                returnval = true;
            }

            return returnval;
        }
    }
}
