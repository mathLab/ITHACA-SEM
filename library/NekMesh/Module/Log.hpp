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

#ifdef _WIN32
#include <io.h>
#define ISSTDOUTTTY _isatty(_fileno(stdout))
#define ISSTDERRTTY _isatty(_fileno(stderr))
#else
#include <unistd.h>
#define ISSTDOUTTTY isatty(fileno(stdout))
#define ISSTDERRTTY isatty(fileno(stderr))
#endif

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

    /**
     * @brief Returns true if this is an interactive (tty) terminal
     * (e.g. stdout).
     */
    bool IsTty() const
    {
        return m_isTty;
    }

protected:
    /// True if this log output is tty-based.
    bool m_isTty = true;

    /**
     * @brief Finalise the log.
     *
     * When called, this function should any finalising required for the log
     * (e.g. closing a file).
     */
    virtual void Finalise() = 0;
};

/**
 * @brief Log output for std::ostream objects.
 */
class StreamOutput : public LogOutput
{
public:
    /**
     * @brief Construct a stream output from the stream @p os.
     *
     * This attempts to determine whether the @p os is either std::cout or
     * std::cerr, so that it can distinguish between tty-type output or
     * file/buffer-based output streams.
     */
    StreamOutput(std::ostream &os) : LogOutput(), m_os(os)
    {
        // Bit of a hack to figure out if we are a tty.
        if ((&os == &std::cout && !ISSTDOUTTTY) ||
            (&os == &std::cerr && !ISSTDERRTTY))
        {
            m_isTty = false;
        }
    }

    /// Default destructor.
    virtual ~StreamOutput() = default;

    /**
     * @brief Writes a log message to @p m_os.
     */
    void Log(const std::string &msg) override
    {
        m_os << msg << std::flush;
    }

private:
    /// Reference to output stream.
    std::ostream &m_os;

    /// Finalise routine. By default we do nothing for streams and assume this
    /// will be handled by the caller.
    void Finalise() override
    {
    }
};

/**
 * @brief Convenience namespace for some basic ANSI terminal codes.
 */
namespace ansi
{
/// Bold text
const std::string bold = "\x1B[1m";
/// Red text
const std::string red = "\033[0;31m";
/// Green text
const std::string green = "\033[1;32m";
/// Yellow text
const std::string yellow = "\033[1;33m";
/// Cyan text
const std::string cyan = "\033[0;36m";
/// Magenta text
const std::string magenta = "\033[0;35m";
/// Reset/remove all formatting
const std::string reset = "\033[0m";
}

/**
 * @brief Basic logging class for the NekMesh library.
 *
 * Usage is like a normal std::ostream, but there is support for a number of
 * additional things:
 *
 * - Various logging levels, defined by the #LogLevel enum.
 * - Automatically preprend a prefix to output messages, which helps in
 *   distinguishing which modules are generating output.
 * - Messages are sent to a #LogOutput class, which handles the physical manner
 *   of handling the string message. This allows for backends that print to
 *   stdout, write to files, or redirect to Python for the NekPy bindings.
 * - Logging levels are defined by the Logger's operator() function, and are
 *   reset to the default LogLevel::INFO level when the logger sees a std::endl
 *   or std::flush.
 * - LogLevel::FATAL messages will trigger a #NekMeshException.
 *
 * Example usage:
 *
 * ```
 * // Select the logger's output level.
 * LogLevel outLevel = INFO;
 * // Create a stream output object to write to std::cout.
 * std::shared_ptr<LogOutput> out = std::make_shared<StreamOutput>(std::cout);
 * // Create the logger.
 * Logger log(out, outLevel);
 * // Write some output
 * log << "This is an 'INFO' level output and will be displayed" << std::endl;
 * log(VERBOSE) << "This is an 'VERBOSE' level output and will not be "
 *              << "displayed: the log level is set at INFO." << std::endl;
 * log(FATAL)   << "This 'FATAL' level output will trigger an exception."
 *              << std::endl;
 * ```
 */
class Logger
{
    typedef std::ostream&  (*ManipFn)(std::ostream&);
    typedef std::ios_base& (*FlagsFn)(std::ios_base&);

public:
    /// Default constructor.
    Logger() = default;

    /**
     * @brief Create a Logger given an output object @p logOutput and a minimum
     * logging level @p level. Log messages above this level will be suppressed.
     */
    Logger(std::shared_ptr<LogOutput> logOutput, LogLevel level) :
        m_logOutput(logOutput), m_level(level)
    {
    }

    /**
     * @brief Copy constructor.
     */
    Logger(const Logger &log)
    {
        m_prefix = log.m_prefix;
        m_level = log.m_level;
        m_curLevel = log.m_curLevel;
        m_logOutput = log.m_logOutput;
        m_prefixLen = log.m_prefixLen;
    }

    /**
     * @brief Assignment operator.
     */
    Logger &operator=(const Logger &log)
    {
        m_prefix = log.m_prefix;
        m_level = log.m_level;
        m_curLevel = log.m_curLevel;
        m_logOutput = log.m_logOutput;
        m_prefixLen = log.m_prefixLen;
        return *this;
    }

    /**
     * @brief Stream operator to act in the same manner as std::ostream.
     *
     * Writes the value @p val internally to a std::stringstream.
     */
    template <typename T>
    Logger &operator<<(const T& val)
    {
        m_buffer << val;
        return *this;
    }

    /**
     * @brief Stream operator to act in the same manner as std::ostream.
     *
     * Applies the manipulator @p manip internally to a std::stringstream.
     */
    Logger &operator<<(ManipFn manip)
    {
        manip(m_buffer);

        if (m_curLevel <= m_level)
        {
            std::stringstream ss;
            ss << GetPrefixString() << m_buffer.str();

            if (m_logOutput->IsTty())
            {
                ss << ansi::reset;
            }

            m_logOutput->Log(ss.str());
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

    /**
     * @brief Stream operator to act in the same manner as std::ostream.
     *
     * Applies the flags function @p manip internally to a std::stringstream.
     */
    Logger& operator<<(FlagsFn manip)
    {
        manip(m_buffer);
        return *this;
    }

    /**
     * @brief Sets the log level of the current message to @p level.
     */
    Logger &operator()(LogLevel level)
    {
        m_curLevel = level;
        return *this;
    }

    /**
     * @brief Prints a progress bar.
     */
    void Progress(const int position, const int goal, const std::string message,
                  int lastprogress = -1)
    {
        std::stringstream ss;
        if (m_logOutput->IsTty())
        {
            float progress = position / float(goal);
            int  numeq = static_cast<int>(ceil(progress *49));

            // Avoid lots of output.
            if (lastprogress == numeq)
            {
                return;
            }

            // carriage return
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
        }
        else
        {
            // print only every 2 percent
            if (int(ceil(double(100 * position / goal))) % 2 ==  0)
            {
                ss << ".";
            }
        }

        m_logOutput->Log(ss.str());
    }

    /**
     * @brief Logs a newline character. Helpful in combination with the
     * #Progress function.
     */
    void Newline()
    {
        m_logOutput->Log("\n");
    }

    /**
     * @brief Overwrites current line if TTY mode is enabled.
     */
    void Overwrite()
    {
        if (m_logOutput->IsTty())
        {
            m_logOutput->Log("\r\033[0K");
        }
    }

    /**
     * @brief Sets the prefix for this Logger that is automatically preprended
     * to all messages.
     *
     * For example, the code
     *
     * ```c++
     * log.SetPrefix("TestModule");
     * log << "Test message" << std::endl;
     * ```
     *
     * would produce the output
     *
     * ```
     * [TestModule]       Test message
     * ```
     */
    void SetPrefix(const std::string &prefix)
    {
        m_prefix = prefix;
    }

    /**
     * @brief Sets the prefix length used to calculate the correct spacing for
     * prefix strings.
     *
     * This is a cosmetic feature intended to reducing spacing between the
     * prefix and the message. For example, the code
     *
     * ```c++
     * log.SetPrefix("TestModule");
     * log.SetPrefixLen(10);
     * log << "Test message" << std::endl;
     * log.SetPrefixLen(15);
     * log << "Test message" << std::endl;
     * ```
     *
     * would produce the output:
     *
     * ```
     * [TestModule] Test message
     * [TestModule]      Test message
     * ```
     */
    void SetPrefixLen(size_t prefixLen)
    {
        m_prefixLen = prefixLen;
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
    /// Default width for prefix.
    size_t m_prefixLen = 0;

    /**
     * @brief Get the prefix string. This is a helper function that also applies
     * basic formatting for the message.
     */
    std::string GetPrefixString()
    {
        std::stringstream ss;
        std::string msgcolour = "", msgprefix = "";

        switch(m_curLevel)
        {
            case FATAL:
                msgcolour = ansi::red;
                msgprefix = "ERROR: ";
                break;
            case WARNING:
                msgcolour = ansi::yellow;
                msgprefix = "WARNING: ";
                break;
            default:
                break;
        }

        size_t prefixLen = m_prefixLen == 0 ? 20 : m_prefixLen + 3;

        if (m_prefix != "")
        {
            if (m_logOutput->IsTty())
            {
                ss << msgcolour << ansi::bold;
            }

            ss << std::setw(prefixLen) << std::left << ("[" + m_prefix + "]");

            if (m_logOutput->IsTty())
            {
                ss << ansi::reset;
            }
        }

        if (m_logOutput->IsTty())
        {
            ss << msgcolour;
        }

        ss << msgprefix;

        return ss.str();
    }
};

}
}

#endif
