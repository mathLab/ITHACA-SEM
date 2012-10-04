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
// Description: Implementation of the metric base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <Metric.h>
#include <loki/Singleton.h>

namespace Nektar
{
    MetricFactory& GetMetricFactory()
    {
        typedef Loki::SingletonHolder<MetricFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }
    
    /**
     * @brief Constructor.
     */
    Metric::Metric(int id) : m_id(id)
    {
        
    }
    
    /**
     * @brief Parse the contents of a metric tag. This function is implemented by
     * subclasses.
     */
    void Metric::Parse(TiXmlElement *metric)
    {
        v_Parse(metric);
    }

    /**
     * @brief Test a line of output from an executible.
     */
    bool Metric::TestLine(std::string line)
    {
        return v_TestLine(line);
    }
    
    /**
     * @brief Test which is run after the executible has finished.
     */
    bool Metric::FinishTest()
    {
        return v_FinishTest();
    }
    
    bool Metric::v_TestLine(std::string line)
    {
        return true;
    }
    
    bool Metric::v_FinishTest()
    {
        return true;
    }
}
