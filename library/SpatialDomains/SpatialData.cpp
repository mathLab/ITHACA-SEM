////////////////////////////////////////////////////////////////////////////////
//
//  File: SpatialData.cpp
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
//  Description: Storage for spatially-dependent parameters.
//
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/SpatialData.h>
#include <tinyxml/tinyxml.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * Private default constructor.
         */
/*
        SpatialParameters::SpatialParameters():
                m_session()
        {

        }
*/

        /**
         * Define a set of spatially dependent parameters with \a nq data
         * points.
         * @param   nq          Number of data points in the domain.
         */
        SpatialParameters::SpatialParameters(const LibUtilities::SessionReaderSharedPtr&  pSession, const int nq) :
                m_nq(nq),
                m_spatialMap(),
                m_constMap(),
                m_analyticMap(),
                m_session(pSession)
        {
        }


        /**
         * Define a set of spatially dependent parameters as a deep copy of an
         * existing set of parameters.
         * @param   src         Existing set of parameters to copy.
         */
        SpatialParameters::SpatialParameters(const SpatialParameters& src) :
                m_nq(src.m_nq),
                m_spatialMap(),
                m_constMap(src.m_constMap),
                m_analyticMap(src.m_analyticMap),
                m_session(src.m_session)
        {
            SpatialDataMap::const_iterator x;
            for (x = src.m_spatialMap.begin(); x != src.m_spatialMap.end(); ++x)
            {
                m_spatialMap[x->first] = MemoryManager<SpatialData>
                                        ::AllocateSharedPtr(*(x->second));
            }
        }


        /**
         * Reads spatially-dependent parameters from an XML file.
         * @param   infilename  XML file.
         */
        void SpatialParameters::Read(std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") +
                infilename).c_str());

            Read(doc);
        }


        /**
         * Reads spatially-dependent parameters from an XML object.
         * @param   doc         XML document containing spatial data
         *                      definitions.
         */
        void SpatialParameters::Read(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            TiXmlElement* conditions = NULL;

            /// Look for all data in CONDITIONS block.
            conditions = docHandle.FirstChildElement("NEKTAR")
                                  .FirstChildElement("CONDITIONS").Element();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *spatialData
                                = conditions->FirstChildElement("SPATIALDATA");

            if (!spatialData)
            {
                return;
            }

            // Now read all the different tagged sections
            ReadConstantFunctions(spatialData);
            ReadAnalyticFunctions(spatialData);
        }


        /**
         * Reads constant-valued spatially-dependent parameters. These are
         * defined in the session file with a C tag.
         * @param   spatialData SpatialData XML element.
         */
        void SpatialParameters::ReadConstantFunctions(TiXmlElement *spatialData)
        {
            TiXmlElement *element = spatialData->FirstChildElement("C");

            while (element)
            {
                std::string variableStr = element->Attribute("NAME");
                ASSERTL0(!variableStr.empty(),
                         "The variable must be specified for the constant "
                         "spatial data.");

                std::string fcnStr = element->Attribute("VALUE");
                ASSERTL0(!fcnStr.empty(),
                         (std::string("Constant for var: ") + variableStr
                                + std::string(" must be specified.")).c_str());

                m_constMap[variableStr] = atof(fcnStr.c_str());

                element = element->NextSiblingElement("C");
            }

        }


        /**
         * Reads analytic spatially-dependent function parameters. These are
         * defined in the session file with an A tag.
         * @param   spatialData SpatialData XML element.
         */
        void SpatialParameters::ReadAnalyticFunctions(TiXmlElement *spatialData)
        {
            TiXmlElement *element = spatialData->FirstChildElement("A");

            while (element)
            {
                std::string variableStr = element->Attribute("NAME");
                ASSERTL0(!variableStr.empty(),
                         "The variable must be specified for the analytic "
                         "spatial data.");

                std::string fcnStr = element->Attribute("VALUE");
                ASSERTL0(!fcnStr.empty(),
                         (std::string("Function for var: ") + variableStr
                                + std::string(" must be specified.")).c_str());

                m_analyticMap[variableStr] = fcnStr;

                element = element->NextSiblingElement("A");
            }

        }


        /**
         * Evaluates all spatially-dependent parameters at all quadrature,
         * points specified by the coordinate arrays. This function must be
         * called to populate the list of spatially-dependent data fields.
         * @param   x           x-coordinates of quadrature points.
         * @param   y           y-coordinates of quadrature points.
         * @param   z           z-coordinates of quadrature points.
         */
        void SpatialParameters::EvaluateParameters(
                    const Array<OneD, const NekDouble> x,
                    const Array<OneD, const NekDouble> y,
                    const Array<OneD, const NekDouble> z)
        {
            // First process constant-valued parameters. For this we simply
            // fill the spatial data region with the constant.
            std::map<std::string, NekDouble>::iterator p;
            for (p = m_constMap.begin(); p != m_constMap.end(); p++)
            {
                SpatialDataSharedPtr fn(MemoryManager<SpatialData>
                                                    ::AllocateSharedPtr(m_nq));
                Vmath::Fill(m_nq, p->second, fn->UpdatePhys(), 1);

                m_spatialMap[p->first] = fn;
            }

            // Now process analytic-valued function parameters. For this we
            // evaluate the supplied equation at each quadrature point.
            std::map<std::string, std::string>::iterator q;
            for (q = m_analyticMap.begin(); q != m_analyticMap.end(); q++)
            {
                SpatialDataSharedPtr fn(MemoryManager<SpatialData>
                                                    ::AllocateSharedPtr(m_nq));
                LibUtilities::Equation E(m_session, q->second);
                E.Evaluate(x, y, z, fn->UpdatePhys());

                m_spatialMap[q->first] = fn;
            }
        }
    }
}

