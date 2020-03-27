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

#include <tinyxml.h>
#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include <TestException.hpp>

namespace fs = boost::filesystem;

std::string PortablePath(const boost::filesystem::path& path);

namespace Nektar
{
    /**
     * @brief Check to see whether the given string @p s is empty (or null).
     */
    inline bool EmptyString(const char *s)
    {
        if (!s)
        {
            return true;
        }
        return std::string(s) == "";
    }

    class Metric
    {
    public:
        Metric(TiXmlElement *metric, bool generate);
        
        /// Perform the test, given the standard output and error streams
        bool Test     (std::istream& pStdout, std::istream& pStderr);
        /// Perform the test, given the standard output and error streams
        void Generate (std::istream& pStdout, std::istream& pStderr);
        /// Return metric type
        std::string GetType()
        {
            return m_type;
        }
        /// Return metric ID
        int GetID()
        {
            return m_id;
        }

    protected:
        /// Stores the ID of this metric.
        int m_id;
        /// Stores the type of this metric (uppercase).
        std::string m_type;
        /// Determines whether to generate this metric or not.
        bool m_generate;
        /// Pointer to XML structure containing metric definition.
        TiXmlElement *m_metric;
        
        virtual bool v_Test     (std::istream& pStdout, 
                                 std::istream& pStderr) = 0;
        virtual void v_Generate (std::istream& pStdout, 
                                 std::istream& pSrderr) = 0;
    };

    /// A shared pointer to an EquationSystem object
    typedef std::shared_ptr<Metric> MetricSharedPtr;

    /// Datatype of the NekFactory used to instantiate classes derived from the
    /// Advection class.
    class MetricFactory
    {
    public:
        typedef MetricSharedPtr (*CreatorFunction)(TiXmlElement *, bool);

        std::string RegisterCreatorFunction(std::string key, CreatorFunction func)
        {
            m_map[key] = func;
            return key;
        }

        MetricSharedPtr CreateInstance(
            std::string key, TiXmlElement *elmt, bool generate)
        {
            return m_map[key](elmt, generate);
        }

    private:
        std::map<std::string, CreatorFunction> m_map;
    };

    MetricFactory& GetMetricFactory();
}

#endif
