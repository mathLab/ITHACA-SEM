////////////////////////////////////////////////////////////////////////////////
//
//  File: Domain.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <sstream>

#include <SpatialDomains/Domain.h>

// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <tinyxml/tinyxml.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        Domain::Domain(MeshGraph *meshGraph):
            m_meshGraph(meshGraph)
        {
        }

        Domain::~Domain()
        {
        }

         void Domain::Read(std::string &infilename)
        {
            SetFileName(infilename);
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") + infilename).c_str());

            Read(doc);
        }

        // \brief Read will read the meshgraph vertices given a TiXmlDocument.
        void Domain::Read(TiXmlDocument &doc)
        {
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.
            TiXmlElement* domain = NULL;
            TiXmlElement* boundary = NULL;

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file");

            domain = master->FirstChildElement("DOMAIN");

            ASSERTL0(domain, "Unable to find DOMAIN tag in file.");

            TiXmlNode *node = domain->FirstChild();
            std::string components(node->ValueStr());
            std::istringstream domainStrm(components);
            std::string entry;

            // May have multiple composites defined.
            if(domainStrm)
            {
                domainStrm >> entry;

                if (entry.length())
                {
                    ASSERTL0(entry[0] == 'C', "Only composites are allowed in a definition of a domain.");

                    // Determine the index associated with the C.  Allow range syntax.

                    std::string::size_type indxBeg = entry.find_first_of('[') + 1;
                    std::string::size_type indxEnd = entry.find_last_of(']') - 1;

                    ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading DOMAIN definition:") + entry).c_str());

                    std::string indxStr = entry.substr(indxBeg, indxEnd - indxBeg + 1);

                    std::istringstream indexStrm(indxStr);
                    int indx1=-1, indx2=-1;
                    if (indexStrm)
                    {
                        // Should read either [a] where a is a nonnegative integer, or
                        // [a-b] where a and b are nonnegative integers, b>a.
                        // Easiest way to know is if a '-' is present we have the latter
                        // case.

                        indexStrm >> indx1;
                        ASSERTL0(indx1 >= 0, (std::string("Error reading collection range: ") + indxStr).c_str());
                        indx2 = indx1;

                        std::string::size_type dashLoc=indxStr.find('-');
                        if (dashLoc != std::string::npos)
                        {
                            // Skip up to and including the '-' character, then read
                            // the other index.  We are safe in doing this because we
                            // already know it is there...somewhere.
                            indexStrm.seekg(dashLoc+1);
                            indexStrm >> indx2;

                            ASSERTL0(indx1 < indx2 && indx2 >= 0,
                                (std::string("Error reading collection range: ") + indxStr).c_str());
                        }

                        for (int i=indx1; i<=indx2; ++i)
                        {
                            Composite composite = m_meshGraph->GetComposite(i);
                            m_domain.push_back(composite);
                        }
                    }
                }
            }

            boundary = master->FirstChildElement("BOUNDARY");
            ASSERTL0(boundary, "Unable to find BOUNDARY tag in file.");

            TiXmlElement *bc = boundary->FirstChildElement();

            // Boundary will have type and composite list.
            while(bc)
            {
                std::string bcType(bc->ValueStr());

                TiXmlNode *node = bc->FirstChild();
                std::string components(node->ValueStr());
                std::istringstream boundaryStrm(components);
                std::string entry;

                // Index of the tag letter into the type enum.
                const char *beginName = BoundaryTypeNameMap;
                // std::find needs the end to be one past the last element.
                const char *endName = BoundaryTypeNameMap+eBoundaryTypeLastElement + 1;
                const char* indx = std::find(beginName, endName, bcType[0]);

                // Not found if the index (ptr) is past the last element.
                ASSERTL0(indx != endName, (std::string("Unable to read boundary type tag: ") + bcType).c_str());

                BoundarySharedPtr boundary(new BoundaryEntry);
                boundary->m_BoundaryType = BoundaryType(indx - BoundaryTypeNameMap);
                m_boundaries.push_back(boundary);

                // May have multiple composites defined.
                if(boundaryStrm)
                {
                    domainStrm >> entry;

                    if (entry.length())
                    {
                        ASSERTL0(entry[0] == 'C', "Only composites are allowed in a definition of a boundary condition.");

                        // Determine the index associated with the C.  Allow range syntax.

                        std::string::size_type indxBeg = entry.find_first_of('[') + 1;
                        std::string::size_type indxEnd = entry.find_last_of(']') - 1;

                        ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading BOUNDARY definition:") + entry).c_str());

                        // Read between the brackets.
                        std::string indxStr = entry.substr(indxBeg, indxEnd - indxBeg + 1);

                        std::istringstream indexStrm(indxStr);
                        int indx1=-1, indx2=-1;
                        if(indexStrm)
                        {
                            // Should read either [a] where a is a nonnegative integer, or
                            // [a-b] where a and b are nonnegative integers, b>a.
                            // Easiest way to know is if a '-' is present we have the latter
                            // case.

                            indexStrm >> indx1;
                            ASSERTL0(indx1 >= 0, (std::string("Error reading collection range: ") + indxStr).c_str());

                            std::string::size_type dashLoc=indxStr.find('-');
                            if (dashLoc != std::string::npos)
                            {
                                // Skip up to and including the '-' character, then read
                                // the other index.  We are safe in doing this because we
                                // already know it is there...somewhere.
                                while(indexStrm.get() != '-');
                                indexStrm >> indx2;

                                ASSERTL0(indx1 < indx2 && indx2 >= 0,
                                    (std::string("Error reading collection range: ") + indxStr).c_str());
                            }

                            for (int i=indx1; i<=indx2; ++i)
                            {
                                Composite composite = m_meshGraph->GetComposite(i);
                                boundary->m_BoundaryComposites.push_back(composite);
                            }
                        }
                    }
                }

                bc = bc->NextSiblingElement();
            }
        }

        void Domain::Write(std::string &outfilename)
        {
        }
    }
};
