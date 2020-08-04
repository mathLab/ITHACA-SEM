////////////////////////////////////////////////////////////////////////////////
//
//  File: Log.hpp
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
//  Description: Simple logging class for NekMesh modules.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_MODULE_LOG_HPP
#define NEKMESH_MODULE_LOG_HPP

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>

namespace Nektar
{
namespace NekMesh
{

/**
 * @brief Log level outputs.
 */
enum LogLevel
{
    FATAL,     ///< Fatal errors that trigger exceptions.
    WARNING,   ///< Warnings (i.e. recoverable errors).
    INFO,      ///< Normal output.
    VERBOSE,   ///< Verbose output.
    TRACE      ///< Debug-level output.
};

/**
 * @brief Custom exception class for all NekMesh-related errors.
 */
class NekMeshError : public std::runtime_error
{
public:
    NekMeshError(const std::string& message) : std::runtime_error(message)
    {
    }
};

/**
 * @brief Base class for log output handling.
 *
 * This class defines the abstract concept of log writing, where subclasses
 * implement concrete implementations for e.g. writing to stdout or a file.
 */
class LogOutput
{
public:
    /// Default constructor.
    LogOutput() = default;
    virtual ~LogOutput() = default;

    /**
     * @brief Write a log message.
     *
     * Subclasses should override this function to write out the log message
     * defined in @p msg.
     *
     * @param msg   Log message to write.
     */
    virtual void Log(const std::string &msg) = 0;

protected:
    /**
     * @brief Finalise the log.
     *
     * When called, this function should any finalising required for the log
     * (e.g. closing a file).
     */
    virtual void Finalise() = 0;
};

class StreamOutput : public LogOutput
{
public:
    StreamOutput(std::ostream &os) : LogOutput(), m_os(os)
    {
    }

    virtual ~StreamOutput() = default;

    void Log(const std::string &msg) override
    {
        m_os << msg << std::flush;
    }

private:
    std::ostream &m_os;

    void Finalise() override
    {
    }
};

class Logger
{
    typedef std::ostream&  (*ManipFn)(std::ostream&);
    typedef std::ios_base& (*FlagsFn)(std::ios_base&);

public:
    Logger() = default;

    Logger(std::shared_ptr<LogOutput> logOutput, LogLevel level) :
        m_logOutput(logOutput), m_level(level)
    {
    }

    Logger(const Logger &log)
    {
        m_prefix = log.m_prefix;
        m_level = log.m_level;
        m_curLevel = log.m_curLevel;
        m_logOutput = log.m_logOutput;
    }

    Logger &operator=(const Logger &log)
    {
        m_prefix = log.m_prefix;
        m_level = log.m_level;
        m_curLevel = log.m_curLevel;
        m_logOutput = log.m_logOutput;
        return *this;
    }

    template <typename T>
    Logger &operator<<(const T& val)
    {
        m_buffer << val;
        return *this;
    }

    Logger &operator<<(ManipFn manip)
    {
        manip(m_buffer);

        if (m_curLevel <= m_level)
        {
            std::stringstream tmp;
            const std::string reset("\033[0m");
            tmp << GetPrefixString() << m_buffer.str() << reset;
            m_logOutput->Log(tmp.str());
        }

        // For fatal exceptions, store temporary message.
        std::string msg;
        if (m_curLevel == FATAL)
        {
            msg = m_buffer.str();
        }

        // Reset buffer & log level.
        m_buffer.str(std::string());
        m_buffer.clear();

        if (m_curLevel == FATAL)
        {
            m_curLevel = INFO;

            // Throw exception.
            throw NekMeshError(msg);
        }
        m_curLevel = INFO;

        return *this;
    }

    Logger& operator<<(FlagsFn manip)
    {
        manip(m_buffer);
        return *this;
    }

    Logger &operator()(LogLevel level)
    {
        m_curLevel = level;
        return *this;
    }

    void Progress(const int position, const int goal, const std::string message,
                  int lastprogress = -1)
    {
        float progress = position / float(goal);
        int  numeq = static_cast<int>(ceil(progress *49));

        // Avoid lots of output.
        if (lastprogress == numeq)
        {
            return;
        }

        // carriage return
        std::stringstream ss;
        ss << "\r" << GetPrefixString()
           << message << ": "
           << std::setw(3) << ceil(100 * progress) << "% [";

        for (int j = 0; j < numeq; j++)
        {
            ss << "=";
        }
        for (int j = numeq; j < 49; j++)
        {
            ss << " ";
        }
        ss << "]";

        m_logOutput->Log(ss.str());
    }

    void Newline()
    {
        m_logOutput->Log("\n");
    }

    void Overwrite()
    {
        m_logOutput->Log("\r\e[0K");
    }

    void SetPrefix(const std::string &prefix)
    {
        m_prefix = prefix;
    }

private:
    /// Log output method.
    std::shared_ptr<LogOutput> m_logOutput;
    /// A buffer that will hold message output until std::endl or std::flush is
    /// seen in the ostream.
    std::stringstream m_buffer;
    /// The level after which to suppress output, i.e. if set to INFO then
    /// VERBOSE and DEBUG messages will not be output.
    LogLevel m_level = INFO;
    /// The level currently being used for logging.
    LogLevel m_curLevel = INFO;
    /// A prefix that will be enclosed with `[` and `]` in the output message.
    std::string m_prefix = "";

    std::string GetPrefixString()
    {
        std::stringstream ss;

        const std::string bold("\x1B[1m");
        const std::string red("\033[0;31m");
        const std::string green("\033[1;32m");
        const std::string yellow("\033[1;33m");
        const std::string cyan("\033[0;36m");
        const std::string magenta("\033[0;35m");
        const std::string reset("\033[0m");

        std::string msgcolour = "";
        std::string msgprefix = "";

        switch(m_curLevel)
        {
            case FATAL:
                msgcolour = red;
                msgprefix = "ERROR: ";
                break;
            case WARNING:
                msgcolour = yellow;
                msgprefix = "WARNING: ";
                break;
            default:
                break;
        }

        if (m_prefix != "")
        {
            ss << msgcolour << bold << std::setw(20) << std::left
               << ("[" + m_prefix + "]") << reset;
        }

        ss << msgcolour <<  msgprefix;

        return ss.str();
    }
};

}
}

#endif
