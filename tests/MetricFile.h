///////////////////////////////////////////////////////////////////////////////
//
// File: MetricFile.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Definition of the LInf metric.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRICLINF_H
#define NEKTAR_TESTS_METRICLINF_H

#include <Metric.h>
#include <map>

namespace Nektar
{
    class MetricFile : public Metric
    {
    public:
        virtual ~MetricFile() {}

        static MetricSharedPtr create(TiXmlElement *metric, bool generate)
        {
            return MetricSharedPtr(new MetricFile(metric, generate));
        }

        static std::string type;

    protected:
        MetricFile(TiXmlElement *metric, bool generate);

        std::string CalculateHash(std::string filename);
        
        virtual bool v_Test    (std::istream& pStdout, std::istream& pStderr);
        virtual void v_Generate(std::istream& pStdout, std::istream& pStderr);
        
        /// Stores filenames to perform hash on.
        std::map<std::string, std::string> m_filehash;
    };
}

#endif
