///////////////////////////////////////////////////////////////////////////////
//
// File: Metric.cpp
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
// Description: Implementation of the metric base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <Metric.h>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    MetricFactory& GetMetricFactory()
    {
        static MetricFactory instance;
        return instance;
    }
    
    /**
     * @brief Constructor.
     */
    Metric::Metric(TiXmlElement *metric, bool generate) :
        m_generate(generate), m_metric(metric)
    {
        if (!metric->Attribute("id"))
        {
            cerr << "Metric has no ID" << endl;
        }
        if (!metric->Attribute("type"))
        {
            cerr << "Metric has no type" << endl;
        }
        m_id = atoi(metric->Attribute("id"));
        m_type = boost::to_upper_copy(string(metric->Attribute("type")));
    }

    /**
     * @brief Test a line of output from an executible.
     */
    bool Metric::Test(std::istream& pStdout, std::istream& pStderr)
    {
        if (m_generate)
        {
            v_Generate(pStdout, pStderr);
            return true;
        }
        else
        {
            return v_Test(pStdout, pStderr);
        }
    }
}
