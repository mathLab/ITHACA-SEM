///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <Collections/Collection.h>

namespace Nektar {
    namespace Collections {
        class CollectionOptimisation
        {
            typedef pair<LibUtilities::ShapeType, int> ElmtOrder;

        public:
            CollectionOptimisation(
                LibUtilities::SessionReaderSharedPtr pSession,
                ImplementationType defaultType = eIterPerExp)
            {
                map<ElmtOrder, ImplementationType> defaults;
                map<ElmtOrder, ImplementationType>::iterator it;

                // Default all elements to eIterPerExp
                defaults[ElmtOrder(LibUtilities::eSegment,       -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::eTriangle,      -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::eQuadrilateral, -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::eTetrahedron,   -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::ePyramid,       -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::ePrism,         -1)] = defaultType;
                defaults[ElmtOrder(LibUtilities::eHexahedron,    -1)] = defaultType;

                map<string, LibUtilities::ShapeType> elTypes;
                map<string, LibUtilities::ShapeType>::iterator it2;
                elTypes["S"] = LibUtilities::eSegment;
                elTypes["T"] = LibUtilities::eTriangle;
                elTypes["Q"] = LibUtilities::eQuadrilateral;
                elTypes["A"] = LibUtilities::eTetrahedron;
                elTypes["P"] = LibUtilities::ePyramid;
                elTypes["R"] = LibUtilities::ePrism;
                elTypes["H"] = LibUtilities::eHexahedron;

                map<string, OperatorType> opTypes;
                for (int i = 0; i < SIZE_OperatorType; ++i)
                {
                    opTypes[OperatorTypeMap[i]] = (OperatorType)i;
                    m_global[(OperatorType)i] = defaults;
                }

                map<string, ImplementationType> impTypes;
                for (int i = 0; i < SIZE_ImplementationType; ++i)
                {
                    impTypes[ImplementationTypeMap[i]] = (ImplementationType)i;
                }

                if(pSession.get()) // turn off file reader if dummy pointer is given
                {
                    TiXmlDocument &doc = pSession->GetDocument();
                    TiXmlHandle docHandle(&doc);
                    TiXmlElement *master
                        = docHandle.FirstChildElement("NEKTAR").Element();
                    ASSERTL0(master, "Unable to find NEKTAR tag in file.");
                    
                    TiXmlElement *xmlCol = master->FirstChildElement("COLLECTIONS");
                    
                    if (xmlCol)
                    {
                        TiXmlElement *elmt = xmlCol->FirstChildElement();
                        
                        while (elmt)
                        {
                            string tagname = elmt->ValueStr();
                            
                            ASSERTL0(boost::iequals(tagname, "OPERATOR"),
                                     "Only OPERATOR tags are supported inside the "
                                     "COLLECTIONS tag.");
                            
                            const char *attr = elmt->Attribute("TYPE");
                            ASSERTL0(attr, "Missing TYPE in OPERATOR tag.");
                            string opType(attr);
                            
                            ASSERTL0(opTypes.count(opType) > 0,
                                     "Unknown OPERATOR type " + opType + ".");
                            
                            OperatorType ot = opTypes[opType];
                            
                            TiXmlElement *elmt2 = elmt->FirstChildElement();
                            
                            while (elmt2)
                            {
                                string tagname = elmt2->ValueStr();
                                ASSERTL0(boost::iequals(tagname, "ELEMENT"),
                                         "Only ELEMENT tags are supported inside the "
                                         "OPERATOR tag.");
                                
                                const char *attr = elmt2->Attribute("TYPE");
                                ASSERTL0(attr, "Missing TYPE in ELEMENT tag.");
                                
                                string elType(attr);
                                it2 = elTypes.find(elType);
                                ASSERTL0(it2 != elTypes.end(),
                                         "Unknown element type "+elType+" in ELEMENT "
                                         "tag");
                                
                                const char *attr2 = elmt2->Attribute("IMPTYPE");
                                ASSERTL0(attr2, "Missing IMPTYPE in ELEMENT tag.");
                                string impType(attr2);
                                ASSERTL0(impTypes.count(impType) > 0,
                                         "Unknown IMPTYPE type " + impType + ".");
                                
                                const char *attr3 = elmt2->Attribute("ORDER");
                                ASSERTL0(attr3, "Missing ORDER in ELEMENT tag.");
                                string order(attr3);

                                if (order == "*")
                                {
                                    m_global[ot][ElmtOrder(it2->second, -1)] 
                                        = impTypes[impType];
                                }
                                else
                                {
                                    vector<unsigned int> orders;
                                    ParseUtils::GenerateSeqVector(order.c_str(), orders);
                                    
                                    for (int i = 0; i < orders.size(); ++i)
                                    {
                                        m_global[ot][ElmtOrder(it2->second, orders[i])] = impTypes[impType];
                                    }
                                }
                                
                                elmt2 = elmt2->NextSiblingElement();
                            }
                            
                            elmt = elmt->NextSiblingElement();
                        }
                    }
                }
#if 0 
                // Print out operator map
                map<OperatorType, map<ElmtOrder, ImplementationType> >::iterator mIt;
                map<ElmtOrder, ImplementationType>::iterator eIt;
                for (mIt = m_global.begin(); mIt != m_global.end(); mIt++)
                {
                    cout << "Operator " << OperatorTypeMap[mIt->first] << ":"
                         << endl;

                    for (eIt = mIt->second.begin(); eIt != mIt->second.end(); eIt++)
                    {
                        cout << "- " << LibUtilities::ShapeTypeMap[eIt->first.first]
                             << " order " << eIt->first.second
                             << " -> " << ImplementationTypeMap[eIt->second]
                             << endl;
                    }
                }
#endif
            }

            OperatorImpMap GetOperatorImpMap(StdRegions::StdExpansionSharedPtr pExp)
            {
                map<OperatorType, map<ElmtOrder, ImplementationType> >::iterator it;
                map<ElmtOrder, ImplementationType>::iterator it2;

                OperatorImpMap ret;
                ElmtOrder searchKey(pExp->DetShapeType(),
                                    pExp->GetBasisNumModes(0));
                ElmtOrder defSearch(pExp->DetShapeType(), -1);

                for (it = m_global.begin(); it != m_global.end(); ++it)
                {
                    ImplementationType impType;

                    it2 = it->second.find(searchKey);

                    if (it2 == it->second.end())
                    {
                        it2 = it->second.find(defSearch);
                        if (it2 == it->second.end())
                        {
                            cout << "shouldn't be here..." << endl;
                            impType = eIterPerExp;
                        }
                        else
                        {
                            impType = it2->second;
                        }
                    }
                    else
                    {
                        impType = it2->second;
                    }

                    ret[it->first] = impType;
                }

                return ret;
            }

        private:
            map<OperatorType, map<ElmtOrder, ImplementationType> > m_global;
        };
    }
}
