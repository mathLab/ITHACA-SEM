///////////////////////////////////////////////////////////////////////////////
//
// File ErrorUtil.hpp
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
// Description: error related utilities
//
///////////////////////////////////////////////////////////////////////////////
#ifndef ERRORUTIL_HPP
#define ERRORUTIL_HPP

#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/LibUtilitiesDeclspec.h>

#if defined(NEKTAR_USE_MPI)
#include <mpi.h>
#endif

#ifndef _WIN32
#include <execinfo.h>
#endif

namespace Nektar
{

class ErrorUtil
{
public:
    class NekError : public std::runtime_error
    {
    public:
        NekError(const std::string& message) : std::runtime_error(message)
        {
        }
    };

    enum ErrType
    {
        efatal,
        ewarning
    };

    inline static void SetErrorStream(std::ostream& o)
    {
        m_outStream = &o;
    }

    inline static void SetPrintBacktrace(bool b)
    {
        m_printBacktrace = b;
    }

    inline static bool HasCustomErrorStream()
    {
        return m_outStream != &std::cerr;
    }

    inline static void Error(ErrType       type,
                             const char   *routine,
                             int           lineNumber,
                             const char   *msg,
                             unsigned int  level,
                             bool          DoComm = false)
    {
        boost::ignore_unused(DoComm);

        // The user of outStream is primarily for the unit tests.  The unit
        // tests often generate errors on purpose to make sure invalid usage is
        // flagged appropriately.  Printing the error messages to cerr made the
        // unit test output hard to parse.

        std::string baseMsg = "Level " + std::to_string(level) +
            " assertion violation\n";
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
        baseMsg += "Where   : " + std::string(routine) + "[" +
            std::to_string(lineNumber) + "]\nMessage : ";
#else
        boost::ignore_unused(routine, lineNumber);
#endif
        baseMsg += std::string(msg);

        // Default rank is zero. If MPI used and initialised, populate with
        // the correct rank. Messages are only printed on rank zero.
        int rank = 0;
#if defined(NEKTAR_USE_MPI) && !defined(NEKTAR_USE_CWIPI)
        int flag = 0;
        if(DoComm)
        {
            MPI_Initialized(&flag);
            if(flag)
            {
                MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            }
        }
#else
        boost::ignore_unused(DoComm);
#endif

        std::string btMessage("");
#if defined(NEKTAR_FULLDEBUG)
#ifndef _WIN32
        if (m_printBacktrace)
        {
            void *btArray[40];
            int btSize;
            char **btStrings;

            btSize = backtrace(btArray, 40);
            btStrings = backtrace_symbols(btArray, btSize);

            for (int i = 0 ; i < btSize ; ++i)
            {
                btMessage +=  std::string(btStrings[i]) + "\n";
            }
            free(btStrings);
        }
#endif
#endif

        switch (type)
        {
        case efatal:
            if (!rank)
            {
                if (m_printBacktrace)
                {
                    (*m_outStream) << btMessage;
                }
                (*m_outStream) << "Fatal   : " << baseMsg << std::endl;
            }

#if defined(NEKTAR_USE_MPI) && !defined(NEKTAR_USE_CWIPI)
            if(DoComm)
            {
                if (flag)
                {
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
#endif
            throw NekError(baseMsg);
            break;
        case ewarning:
            if (!rank)
            {
                if (m_printBacktrace)
                {
                    (*m_outStream) << btMessage;
                }
                (*m_outStream) << "Warning : " << baseMsg << std::endl;
            }
            break;
        default:
            (*m_outStream) << "Unknown warning type: " << baseMsg << std::endl;
        }
    }

    inline static void Error(ErrType type, const char *routine, int lineNumber, const std::string& msg, unsigned int level)
    {
        Error(type, routine, lineNumber, msg.c_str(), level);
    }

    inline static void Error(ErrType type, const char *routine, int lineNumber, const char *msg)
    {
        Error(type, routine, lineNumber, msg, 0);
    }

private:
    LIB_UTILITIES_EXPORT static std::ostream *m_outStream;
    LIB_UTILITIES_EXPORT static bool          m_printBacktrace;
};

/// Assert Level 0 -- Fundamental assert which
/// is used whether in FULLDEBUG, DEBUG or OPT
/// compilation mode.  This level assert is
/// considered code critical, even under
/// optimized compilation.

#define NEKERROR(type, msg)                                     \
    Nektar::ErrorUtil::Error(type, __FILE__, __LINE__, msg, 0);


#define ROOTONLY_NEKERROR(type, msg)                                    \
    Nektar::ErrorUtil::Error(type, __FILE__, __LINE__, msg, 0, true);

#define ASSERTL0(condition,msg)                                         \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::efatal, __FILE__, __LINE__, msg, 0);     \
    }

#define WARNINGL0(condition,msg)                                        \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::ewarning, __FILE__, __LINE__, msg, 0);   \
    }

/// Assert Level 1 -- Debugging which is used whether in FULLDEBUG or
/// DEBUG compilation mode.  This level assert is designed for aiding
/// in standard debug (-g) mode
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)

#define ASSERTL1(condition,msg)                                         \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::efatal, __FILE__, __LINE__, msg, 1);     \
    }

#define WARNINGL1(condition,msg)                                        \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::ewarning, __FILE__, __LINE__, msg, 1);   \
    }

#else //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
#define ASSERTL1(condition,msg)
#define WARNINGL1(condition,msg)
#endif //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)


/// Assert Level 2 -- Debugging which is used FULLDEBUG compilation
/// mode.  This level assert is designed to provide addition safety
/// checks within the code (such as bounds checking, etc.).
#ifdef NEKTAR_FULLDEBUG

#define ASSERTL2(condition,msg)                                         \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::efatal, __FILE__, __LINE__, msg, 2);     \
    }
#define WARNINGL2(condition,msg)                                        \
    if(!(condition))                                                    \
    {                                                                   \
        Nektar::ErrorUtil::Error(                                       \
            Nektar::ErrorUtil::ewarning, __FILE__, __LINE__, msg, 2);   \
    }

#else //NEKTAR_FULLDEBUG
#define ASSERTL2(condition,msg)
#define WARNINGL2(condition,msg)
#endif //NEKTAR_FULLDEBUG

}

#endif //ERRORUTIL_HPP

