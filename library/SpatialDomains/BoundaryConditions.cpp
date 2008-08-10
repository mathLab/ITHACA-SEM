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
#include <SpatialDomains/BoundaryConditions.h>
#include <SpatialDomains/ParseUtils.hpp>

#include <string>

namespace Nektar
{
    namespace SpatialDomains
    {
        ParamMap BoundaryConditions::m_Parameters;

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
            
            ReadFunctions(conditions);

            ReadVariables(conditions);
            
            ReadBoundaryRegions(conditions);

            ReadBoundaryConditions(conditions);

            ReadForcingFunctions(conditions);

            ReadInitialConditions(conditions);

            ReadExactSolution(conditions);
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
                            m_Parameters[lhs] = value;
                            expEvaluator.SetParameter(lhs, value);
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }

                // Set ourselves up for evaluation later.
                Equation::SetConstParameters(m_Parameters);
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

                    FunctionMap::iterator fcnIter = m_Functions.find(lhsString);
                    ASSERTL0(fcnIter == m_Functions.end(),
                        (std::string("Function value: ") + lhsString + std::string("already specified.")).c_str());

                    m_Functions[lhsString] = fcnString;
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
                    m_Variables.push_back(variableName);

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
                    m_MeshGraph->GetCompositeList(indxStr, *boundaryRegion);

                    m_BoundaryRegions.push_back(boundaryRegion);
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
                    TiXmlAttribute *attr = conditionElement->FirstAttribute();

                    Variable::iterator iter;
                    std::string attrName;

                    attrData = conditionElement->Attribute("VAR");

                    if (!attrData.empty())
                    {
                        iter = std::find(m_Variables.begin(), m_Variables.end(), attrData);
                        ASSERTL0(iter != m_Variables.end(), (std::string("Cannot find variable: ") + attrData).c_str());
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

                                        cout << "userdefined Neumann data = " << userDefined << endl;

                                     }
                                     else if(attrName=="VALUE")
                                     {
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());
                                        
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must be specified.");

                                        SubstituteFunction(attrData);
                                        
                                        equation = attrData;

                                        cout << "equation in the Neumann VALUE = " << equation << endl;
                                        
                                      }
                                      attr = attr->Next();
                                    }
                                                                        
                                    BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(equation, userDefined));
                                    (*boundaryConditions)[*iter]  = neumannCondition;

                                    cout << "equation in the Neumann SharedPtr = " << equation << endl;
                                    cout << "user defined Neumann = " << userDefined << endl;
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
                            
                                        cout << "userdefined Dirichlet data = " << userDefined << endl;

                                    }
                                    else if(attrName=="VALUE")
                                    {                                        
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                        SubstituteFunction(attrData);

                                        equation = attrData;

                                        cout << "equation Dirichlet in the VALUE = " << equation << endl;
                                    }
                                   attr = attr->Next();
                                }

                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(equation, userDefined));
                                (*boundaryConditions)[*iter]  = dirichletCondition;
                                cout << "equation in the Dirichlet SharedPtr = " << equation << endl;
                                cout << "user defined Dirichlet = " << userDefined << endl;
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
                                std::string equation1, equation2, userDefined;

                                while(attr){

                                attrName1 = attr->Name();
                               
                                if (attrName1=="USERDEFINEDTYPE") {
                                
                                    // Do stuff for the user defined attribute
                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                    SubstituteFunction(attrData1);

                                    userDefined = attrData1;

                                    cout << "user defined Robin data = " << userDefined << endl;

                                 } else if(attrName1 == "A"){

                                    ASSERTL0(attrName1 == "A", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "A attributes must have associated values.");

                                    SubstituteFunction(attrData1);

                                    equation1 = attrData1;

                                    cout << "equation in the A = " << equation1 << endl;

                                    attr = attr->Next();
                                    ASSERTL0(attr, "Unable to read B attribute.");
                                
                                    attrName1= attr->Name();
                                    ASSERTL0(attrName1 == "B", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "B attributes must have associated values.");

                                    SubstituteFunction(attrData1);

                                    equation2 = attrData1;

                                    cout << "equation in the B = " << equation2 << endl;

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

                                for (Variable::iterator varIter = m_Variables.begin();
                                     varIter != m_Variables.end(); ++varIter)
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
                    std::string variableStr = forcingFunction->Attribute("VAR");
                    ASSERTL0(!variableStr.empty(), "The variable must be specified for the forcing function.");

                    std::string fcnStr = forcingFunction->Attribute("VALUE");
                    ASSERTL0(!fcnStr.empty(),
                        (std::string("Forcing function for var: ") + variableStr + std::string(" must be specified.")).c_str());

                    /// Check the RHS against the functions defined in m_Functions.  If the name on the RHS
                    /// is found in the function map, then use the function contained in the function
                    /// map in place of the RHS.
                    SubstituteFunction(fcnStr);

                    ForcingFunctionsMap::iterator forcingFcnsIter = m_ForcingFunctions.find(variableStr);

                    if (forcingFcnsIter != m_ForcingFunctions.end())
                    {
                        m_ForcingFunctions[variableStr]->SetEquation(fcnStr);
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

            if (exactSolutionElement)
            {
                TiXmlElement *exactSolution = exactSolutionElement->FirstChildElement("F");

                // All exact solution functions are initialized to "0" so they only have to be
                // partially specified.  That is, not all variables have to have functions
                // specified.  For those that are missing it is assumed they are "0".
                for (Variable::iterator varIter = m_Variables.begin();
                    varIter != m_Variables.end(); ++varIter)
                {
                    ExactSolutionShPtr exactSolutionShPtr(MemoryManager<ExactSolution>::AllocateSharedPtr("0"));
                    m_ExactSolution[*varIter] = exactSolutionShPtr;
                }

                while (exactSolution)
                {
                    std::string variableStr = exactSolution->Attribute("VAR");
                    ASSERTL0(!variableStr.empty(), "The variable must be specified for the exact solution.");

                    std::string fcnStr = exactSolution->Attribute("VALUE");
                    ASSERTL0(!fcnStr.empty(),
                        (std::string("The exact solution function must be specified for variable: ") + variableStr).c_str());

                    /// Check the RHS against the functions defined in m_Functions.  If the name on the RHS
                    /// is found in the function map, then use the function contained in the function
                    /// map in place of the RHS.
                    SubstituteFunction(fcnStr);

                    ExactSolutionMap::iterator exactSolutionIter = m_ExactSolution.find(variableStr);

                    if (exactSolutionIter != m_ExactSolution.end())
                    {
                        m_ExactSolution[variableStr]->SetEquation(fcnStr);
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Error setting exact solution for variable: ") + variableStr).c_str());
                    }

                    exactSolution = exactSolution->NextSiblingElement("F");
                }
            }
        }

        void BoundaryConditions::ReadInitialConditions(TiXmlElement *conditions)
        {
            TiXmlElement *initialConditionsElement = conditions->FirstChildElement("INITIALCONDITIONS");

            if (initialConditionsElement)
            {
                TiXmlElement *initialCondition = initialConditionsElement->FirstChildElement("F");

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
                    std::string variableStr = initialCondition->Attribute("VAR");
                    ASSERTL0(!variableStr.empty(), "The variable must be specified for the exact solution.");

                    std::string fcnStr = initialCondition->Attribute("VALUE");
                    ASSERTL0(!fcnStr.empty(),
                        (std::string("The initial condition function must be specified for variable: ") + variableStr).c_str());

                    /// Check the RHS against the functions defined in m_Functions.  If the name on the RHS
                    /// is found in the function map, then use the function contained in the function
                    /// map in place of the RHS.
                    SubstituteFunction(fcnStr);

                    InitialConditionsMap::iterator initialConditionFcnsIter =
                        m_InitialConditions.find(variableStr);

                    if (initialConditionFcnsIter != m_InitialConditions.end())
                    {
                        m_InitialConditions[variableStr]->SetEquation(fcnStr);
                    }
                    else
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Error setting initial condition for variable: ") + variableStr).c_str());
                    }

                    initialCondition = initialCondition->NextSiblingElement("F");
                }
            }
        }

        NekDouble BoundaryConditions::GetParameter(const std::string &parmName)
        {
            ParamMap::iterator paramMapIter = m_Parameters.find(parmName);

            ASSERTL0(paramMapIter != m_Parameters.end(),
                (std::string("Unable to find requested parameter: ") + parmName).c_str());

            return paramMapIter->second;
        }


        const std::string &BoundaryConditions::GetFunction(const std::string &lhs)
        {
            FunctionMap::iterator fcnIter = m_Functions.find(lhs);

            ASSERTL1(fcnIter != m_Functions.end(),
                (std::string("Unable to find requested function: ") + lhs).c_str());

            return fcnIter->second;
        }

        Equation BoundaryConditions::GetFunctionAsEquation(const std::string &lhs)
        {
            FunctionMap::iterator fcnIter = m_Functions.find(lhs);

            ASSERTL1(fcnIter != m_Functions.end(),
                (std::string("Unable to find requested function: ") + lhs).c_str());

            return Equation(fcnIter->second);
        }


        bool BoundaryConditions::SubstituteFunction(std::string &str)
        {
            FunctionMap::iterator fcnIter = m_Functions.find(str);
            bool returnval = false;

            if (fcnIter != m_Functions.end())
            {
                str = fcnIter->second;
                returnval = true;
            }

            return returnval;
        }

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(int indx) const
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

        ConstForcingFunctionShPtr BoundaryConditions::GetForcingFunction(const std::string &var) const
        {
            ConstForcingFunctionShPtr returnval;
            
            // Check that var is defined in forcing function list.
            ForcingFunctionsMap::const_iterator ffIter = m_ForcingFunctions.find(var);

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

        ConstExactSolutionShPtr BoundaryConditions::GetExactSolution(int indx) const
        {
            ConstExactSolutionShPtr returnval;

            if (indx >= m_Variables.size() || indx < 0)
            {
                string errStr;
                std::ostringstream strStream(errStr);
                strStream << indx;

                NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to find variable corresponding to index: ") + errStr).c_str());
            }

            return GetExactSolution(m_Variables[indx]);
        }

        ConstExactSolutionShPtr BoundaryConditions::GetExactSolution(const std::string &var) const
        {
            ConstExactSolutionShPtr returnval;
            
            // Check that var is defined in forcing function list.
            ExactSolutionMap::const_iterator exSolnIter = m_ExactSolution.find(var);

            if(exSolnIter != m_ExactSolution.end())
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

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(int indx) const
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

        ConstInitialConditionShPtr BoundaryConditions::GetInitialCondition(const string &var) const
        {
            ConstInitialConditionShPtr returnval;
            
            // Check that var is defined in forcing function list.
            InitialConditionsMap::const_iterator ffIter = m_InitialConditions.find(var);

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

    }
}
