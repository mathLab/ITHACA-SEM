///////////////////////////////////////////////////////////////////////////////
//
// File: MetricL2.h
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
// Description: Definition of the L2 metric.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRICL2_H
#define NEKTAR_TESTS_METRICL2_H

#include <map>
#include <MetricRegex.h>

namespace Nektar
{
    class MetricL2 : public MetricRegex
    {
    public:
        static MetricSharedPtr create(TiXmlElement *metric, bool generate)
        {
            return MetricSharedPtr(new MetricL2(metric, generate));
        }

        static std::string type;
        static std::string defaultTolerance;

    protected:
        MetricL2(TiXmlElement *metric, bool generate);

        std::map<std::string, std::string> m_varTolerance;
        
        virtual void v_Generate(std::istream& pStdout, std::istream& pStderr);
    };
}

#endif
