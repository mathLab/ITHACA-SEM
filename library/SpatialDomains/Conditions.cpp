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

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <SpatialDomains/Conditions.h>
#include <tinyxml/tinyxml.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        BoundaryConditions::BoundaryConditions(const LibUtilities::SessionReaderSharedPtr &pSession, const MeshGraphSharedPtr &meshGraph)
            : m_meshGraph(meshGraph), 
              m_session  (pSession)
              
        {
            Read(m_session->GetElement("Nektar/Conditions"));
        }

        BoundaryConditions::~BoundaryConditions()
        {
        }


        /**
         *
         */
        void BoundaryConditions::Read(TiXmlElement *conditions)
        {
            ASSERTL0(conditions, "Unable to find CONDITIONS tag in file.");

            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");

            if(boundaryRegions)
            {
                ReadBoundaryRegions(conditions);

                ReadBoundaryConditions(conditions);
            }
        }


        /**
         *
         */
        void BoundaryConditions::ReadBoundaryRegions(TiXmlElement *conditions)
        {
            TiXmlElement *boundaryRegions = conditions->FirstChildElement("BOUNDARYREGIONS");
            ASSERTL0(boundaryRegions, "Unable to find BOUNDARYREGIONS block.");

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


        /**
         *
         */
        void BoundaryConditions::ReadBoundaryConditions(TiXmlElement *conditions)
        {
            // Read REGION tags
            TiXmlElement *boundaryConditionsElement = conditions->FirstChildElement("BOUNDARYCONDITIONS");
            ASSERTL0(boundaryConditionsElement, "Boundary conditions must be specified.");

            TiXmlElement *regionElement = boundaryConditionsElement->FirstChildElement("REGION");

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

                // Here is the boundary region.
                // m_boundaryRegions[boundaryRegionID];

                TiXmlElement *conditionElement = regionElement->FirstChildElement();
                std::vector<std::string> vars = m_session->GetVariables();

                while (conditionElement)
                {
                    // Check type.
                    std::string conditionType = conditionElement->Value();
                    std::string attrData;

                    // All have var specified, or else all variables are zero.
                    TiXmlAttribute *attr = conditionElement->FirstAttribute();

                    std::vector<std::string>::iterator iter;
                    std::string attrName;

                    attrData = conditionElement->Attribute("VAR");

                    if (!attrData.empty())
                    {
                        iter = std::find(vars.begin(), vars.end(), attrData);
                        ASSERTL0(iter != vars.end(), (std::string("Cannot find variable: ") + attrData).c_str());
                    }

                    if (conditionType == "N")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Neumann and are set to zero.
                            for (std::vector<std::string>::iterator varIter = vars.begin();
                                varIter != vars.end(); ++varIter)
                            {
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(m_session,"00.0"));
                                (*boundaryConditions)[*varIter]  = neumannCondition;
                            }
                        }                       
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string equation, userDefined, filename;

                                while(attr) 
                                {

                                    attrName = attr->Name();

                                    if (attrName=="USERDEFINEDTYPE") 
                                    {

                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                        // Suppose to go here?
                                        m_session->SubstituteExpressions(attrData);

                                        userDefined = attrData;
                                     }
                                     else if(attrName=="VALUE")
                                     {
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must be specified.");

                                        m_session->SubstituteExpressions(attrData);

                                        equation = attrData;
                                      }
                                     else if(attrName=="FILE")
                                     {
                                        ASSERTL0(attrName == "FILE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "FILE attribute must be specified.");

                                        m_session->SubstituteExpressions(attrData);

                                        filename = attrData;
                                      }
                                      attr = attr->Next();
                                }
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(m_session, equation, userDefined, filename));
                                (*boundaryConditions)[*iter]  = neumannCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr neumannCondition(MemoryManager<NeumannBoundaryCondition>::AllocateSharedPtr(m_session, "0"));
                                (*boundaryConditions)[*iter]  = neumannCondition;
                            }
                        }
                    }
                    else if (conditionType == "D")
                    {
                        if (attrData.empty())
                        {
                            // All variables are Dirichlet and are set to zero.
                            for (std::vector<std::string>::iterator varIter = vars.begin();
                                varIter != vars.end(); ++varIter)
                            {
                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(m_session, "0"));
                                (*boundaryConditions)[*varIter] = dirichletCondition;
                            }
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();

                            if (attr)
                            {
                                std::string equation, userDefined, filename;

                                while(attr) {

                                   attrName = attr->Name();

                                    if (attrName=="USERDEFINEDTYPE") {

                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                        m_session->SubstituteExpressions(attrData);

                                        userDefined = attrData;
                                    }
                                    else if(attrName=="VALUE")
                                    {
                                        ASSERTL0(attrName == "VALUE", (std::string("Unknown attribute: ") + attrName).c_str());

                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(), "VALUE attribute must have associated value.");

                                        m_session->SubstituteExpressions(attrData);

                                        equation = attrData;
                                    }
                                    else if(attrName=="FILE")
                                    {
                                       ASSERTL0(attrName == "FILE", (std::string("Unknown attribute: ") + attrName).c_str());

                                       attrData = attr->Value();
                                       ASSERTL0(!attrData.empty(), "FILE attribute must be specified.");

                                       m_session->SubstituteExpressions(attrData);

                                       filename = attrData;
                                     }
                                   attr = attr->Next();
                                }

                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(m_session, equation, userDefined, filename));
                                (*boundaryConditions)[*iter]  = dirichletCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr dirichletCondition(MemoryManager<DirichletBoundaryCondition>::AllocateSharedPtr(m_session, "0"));
                                (*boundaryConditions)[*iter]  = dirichletCondition;
                            }
                        }
                    }
                    else if (conditionType == "R") // Read du/dn +  PRIMCOEFF u = VALUE
                    {
                        if (attrData.empty())
                        {
                            // All variables are Robin and are set to zero.
                            for (std::vector<std::string>::iterator varIter = vars.begin();
                                varIter != vars.end(); ++varIter)
                            {
                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(m_session, "0", "0"));
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
                                std::string attrName1;
                                std::string attrData1;
                                std::string equation1, equation2, userDefined;
                                std::string filename;

                                while(attr){

                                attrName1 = attr->Name();

                                if (attrName1=="USERDEFINEDTYPE") {

                                    // Do stuff for the user defined attribute
                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "USERDEFINEDTYPE attribute must have associated value.");

                                    m_session->SubstituteExpressions(attrData1);

                                    userDefined = attrData1;

                                 }
                                 else if(attrName1 == "VALUE"){

                                    ASSERTL0(attrName1 == "VALUE", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "VALUE attributes must have associated values.");

                                    m_session->SubstituteExpressions(attrData1);

                                    equation1 = attrData1;

                                    attr = attr->Next();
                                    ASSERTL0(attr, "Unable to read PRIMCOEFF attribute.");

                                    attrName1= attr->Name();
                                    ASSERTL0(attrName1 == "PRIMCOEFF", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "PRIMCOEFF attributes must have associated values.");

                                    m_session->SubstituteExpressions(attrData1);

                                    equation2 = attrData1;

                                 }
                                 else if(attrName1=="FILE")
                                 {
                                    ASSERTL0(attrName1 == "FILE", (std::string("Unknown attribute: ") + attrName1).c_str());

                                    attrData1 = attr->Value();
                                    ASSERTL0(!attrData1.empty(), "FILE attribute must be specified.");

                                    m_session->SubstituteExpressions(attrData1);

                                    filename = attrData1;
                                 }
                                 attr = attr->Next();

                                }

                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(m_session, equation1, equation2, userDefined, filename));
                                (*boundaryConditions)[*iter]  = robinCondition;
                            }
                            else
                            {
                                // This variable's condition is zero.
                                BoundaryConditionShPtr robinCondition(MemoryManager<RobinBoundaryCondition>::AllocateSharedPtr(m_session, "0", "0"));
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

                                for (std::vector<std::string>::iterator varIter = vars.begin();
                                     varIter != vars.end(); ++varIter)
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
					else if (conditionType == "JUNCTION")
                    {						
                        if (attrData.empty())
                        {
                            ASSERTL0(false, "JUNCTION not set up correctly");
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();
							
                            if (attr)
                            {
                                int P = 0;
								int D1 = 0;
								
                                while(attr) {
									
									attrName = attr->Name();
									
                                    if (attrName=="P") {
										
                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();	
                                        m_session->SubstituteExpressions(attrData);
                                        P = atoi(attrData.c_str());
                                    }
									else if (attrName=="D1") {
						 
										// Do stuff for the user defined attribute
										attrData = attr->Value();
										m_session->SubstituteExpressions(attrData);
										D1 = atoi(attrData.c_str());

									}
									attr = attr->Next();
                                }
								
                                BoundaryConditionShPtr junctionCondition(MemoryManager<JunctionBoundaryCondition>::AllocateSharedPtr(P, D1));
                                (*boundaryConditions)[*iter]  = junctionCondition;
                            }
                        }
						
                    }
					else if (conditionType == "BIFURCATION")
                    {						
                        if (attrData.empty())
                        {
                            ASSERTL0(false, "BIFURCATION not set up correctly");
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();
							
                            if (attr)
                            {
                                int P = 0;
								int D1 = 0;
								int D2 = 0;
								
								std::string userDefined;
								
                                while(attr) {
									
									attrName = attr->Name();
									
                                    if (attrName=="P") {
										
                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();	
                                        m_session->SubstituteExpressions(attrData);
                                        P = atoi(attrData.c_str());
                                    }
									else if (attrName=="D1") {
										
										// Do stuff for the user defined attribute
										attrData = attr->Value();
										m_session->SubstituteExpressions(attrData);
										D1 = atoi(attrData.c_str());
										
									}
									else if (attrName=="D2") {
										
										// Do stuff for the user defined attribute
										attrData = attr->Value();
										m_session->SubstituteExpressions(attrData);
										D2 = atoi(attrData.c_str());
										
									}
									attr = attr->Next();
                                }
								
                                BoundaryConditionShPtr bifurcationCondition(MemoryManager<BifurcationBoundaryCondition>::AllocateSharedPtr(P, D1, D2));
                                (*boundaryConditions)[*iter]  = bifurcationCondition;
                            }
                        }
                    }
					else if (conditionType == "MERGING")
                    {						
                        if (attrData.empty())
                        {
                            ASSERTL0(false, "MERGING flow condition not set up correctly");
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();
							
                            if (attr)
                            {
                                int P = 0;
								int D1 = 0;
								int D2 = 0;
								
								std::string userDefined;
								
                                while(attr) {
									
									attrName = attr->Name();
									
                                    if (attrName=="P") {
										
                                        // Do stuff for the user defined attribute
                                        attrData = attr->Value();	
                                        m_session->SubstituteExpressions(attrData);
                                        P = atoi(attrData.c_str());
                                    }
									else if (attrName=="D1") {
										
										// Do stuff for the user defined attribute
										attrData = attr->Value();
										m_session->SubstituteExpressions(attrData);
										D1 = atoi(attrData.c_str());
										
									}
									else if (attrName=="D2") {
										
										// Do stuff for the user defined attribute
										attrData = attr->Value();
										m_session->SubstituteExpressions(attrData);
										D2 = atoi(attrData.c_str());
										
									}
									attr = attr->Next();
                                }
								
                                BoundaryConditionShPtr mergingCondition(MemoryManager<MergingBoundaryCondition>::AllocateSharedPtr(P, D1, D2));
                                (*boundaryConditions)[*iter]  = mergingCondition;
                            }
                        }
                    }
					
					

                    conditionElement = conditionElement->NextSiblingElement();
                }

                m_boundaryConditions[boundaryRegionID] = boundaryConditions;
                regionElement = regionElement->NextSiblingElement("REGION");
            }
       }
    }
}
