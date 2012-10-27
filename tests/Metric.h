///////////////////////////////////////////////////////////////////////////////
//
// File: Metric.h
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
// Description: Definition of the metric base class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRIC_H
#define NEKTAR_TESTS_METRIC_H

#include <tinyxml/tinyxml.h>
#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>

namespace Nektar
{
    class Metric
    {
    public:
        Metric(TiXmlElement *metric);
        
        bool Test  (std::istream& pStdout, std::istream& pStderr);
        bool FinishTest();
        
    protected:
        // Stores the ID of this metric.
        int m_id;
        
//        virtual void v_Parse     (TiXmlElement *metric) = 0;
        virtual bool v_Test      (std::istream& pStdout, std::istream& pStderr);
//        virtual bool v_FinishTest();

//    private:
//        void Parse     (TiXmlElement *metric);
    };

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<Metric> MetricSharedPtr;
    
    /// Datatype of the NekFactory used to instantiate classes derived from the
    /// Advection class.
    typedef LibUtilities::NekFactory<std::string, Metric, TiXmlElement*> MetricFactory;
    MetricFactory& GetMetricFactory();

}


#endif
