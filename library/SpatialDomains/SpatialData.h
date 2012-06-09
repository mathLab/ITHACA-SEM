////////////////////////////////////////////////////////////////////////////////
//
//  File:  SpatialData.h
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
#ifndef NEKTAR_SPATIALDOMAINS_SPATIALDATA_H
#define NEKTAR_SPATIALDOMAINS_SPATIALDATA_H

#include <string>
#include <map>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

class TiXmlElement;
class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        /// Defines a domain-wide spatially-varying parameter.
        class SpatialData
        {
            public:
                /// Construct new spatially-dependent storage.
                SpatialData(int size) :
                    m_data(size)
                {
                }

                /// Copy existing spatially-dependent storage.
                SpatialData(const SpatialData& src) :
                    m_data(src.m_data)
                {
                }

                /// Returns a const reference to the data.
                SPATIAL_DOMAINS_EXPORT inline const Array<OneD, const NekDouble> &GetPhys()  const;

                /// Returns a non-const reference to the data.
                SPATIAL_DOMAINS_EXPORT inline Array<OneD, NekDouble> &UpdatePhys();

            private:
                /// Default constructor.
                SpatialData();

                /// Storage for spatially-dependent data.
                Array<OneD, NekDouble> m_data;
        };
        /// Shared pointer to a SpatialData object.
        typedef boost::shared_ptr<SpatialData>      SpatialDataSharedPtr;
        /// Vector of SpatialData pointers.
        typedef std::vector< SpatialDataSharedPtr > SpatialDataVector;
        /// Iterator for SpatialData vector.
        typedef std::vector< SpatialDataSharedPtr >::iterator
                                                    SpatialDataVectorIter;
        /// Mapping between name and SpatialData.
        typedef std::map< std::string, SpatialDataSharedPtr > SpatialDataMap;


        /// Constructs a set of spatially-dependent parameters.
        class SpatialParameters
        {
            public:
                /// Define a new set of spatially-dependent parameters.
                SPATIAL_DOMAINS_EXPORT SpatialParameters(const LibUtilities::SessionReaderSharedPtr&  pSession, const int nq);
                /// Copies an existing set of spatially-dependent parameters.
                SPATIAL_DOMAINS_EXPORT SpatialParameters(const SpatialParameters& src);

                /// Reads the set of parameters from an XML file.
                SPATIAL_DOMAINS_EXPORT void Read(std::string &infilename);
                /// Reads the set of parameters from an XML document object.
                SPATIAL_DOMAINS_EXPORT void Read(TiXmlDocument &doc);

                /// Evaluate all the parameters at each quadrature point.
                SPATIAL_DOMAINS_EXPORT void EvaluateParameters(
                        const Array<OneD, const NekDouble> x,
                        const Array<OneD, const NekDouble> y,
                        const Array<OneD, const NekDouble> z);

                /// Determine if a named parameter is defined.
                SPATIAL_DOMAINS_EXPORT inline bool Exists(std::string name);

                /// Retrieve a single spatially-dependent parameter by name.
                SPATIAL_DOMAINS_EXPORT inline SpatialDataSharedPtr& GetData(std::string name);

            private:
                /// Number of quadrature points in the domain.
                int                                 m_nq;

                /// List of SpatialData
                SpatialDataMap                      m_spatialMap;

                /// List of constant-valued parameter definitions.
                std::map<std::string, NekDouble>    m_constMap;

                /// List of analytic parameter definitions.
                std::map<std::string, std::string>  m_analyticMap;

                /// Shared pointer to the current session
                const LibUtilities::SessionReaderSharedPtr&  m_session;


                /// Default constructor.
                // SpatialParameters();

                /// Reads constant-valued parameters from XML object.
                void ReadConstantFunctions(TiXmlElement *spatialData);

                /// Reads analytic parameters from XML object.
                void ReadAnalyticFunctions(TiXmlElement *spatialData);

        };
        /// Shared pointer to a SpatialParameters object.
        typedef boost::shared_ptr<SpatialParameters> SpatialParametersSharedPtr;

        inline const Array<OneD, const NekDouble> &SpatialData::GetPhys()  const
        {
            return m_data;
        }

        inline Array<OneD, NekDouble> &SpatialData::UpdatePhys()
        {
            return m_data;
        }

        inline bool SpatialParameters::Exists(std::string name)
        {
            return (m_spatialMap.find(name) != m_spatialMap.end());
        }

        inline SpatialDataSharedPtr& SpatialParameters::GetData(
                            std::string name)
        {
            return m_spatialMap[name];
        }
    }
}

#endif
