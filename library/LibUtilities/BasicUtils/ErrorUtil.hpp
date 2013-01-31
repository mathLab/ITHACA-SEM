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
// Description: error related utilities
//
///////////////////////////////////////////////////////////////////////////////
#ifndef ERRORUTIL_HPP
#define ERRORUTIL_HPP

#include <iostream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace ErrorUtil
{
    static boost::optional<std::ostream&> outStream;

    inline static void SetErrorStream(std::ostream& o)
    {
        outStream = o;
    }
    
    inline static bool HasCustomErrorStream()
    {
        return outStream;
    }

    enum ErrType
    {
        efatal,
        ewarning
    };

    class NekError : public std::runtime_error
    {
    public:
        NekError(const std::string& message) : std::runtime_error(message) {}
    };
        
    inline static void Error(ErrType type, const char *routine, int lineNumber, const char *msg, unsigned int level)
    {
        // The user of outStream is primarily for the unit tests.
        // The unit tests often generate errors on purpose to make sure
        // invalid usage is flagged appropriately.  Printing the error
        // messages to cerr made the unit test output hard to parse.

        std::string baseMsg = std::string("Level ") +
            boost::lexical_cast<std::string>(level) +
            std::string(" assertion violation\n") +
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            std::string("Where   : ") + boost::lexical_cast<std::string>(routine) +  std::string("[") +  boost::lexical_cast<std::string>(lineNumber) +  std::string("]\n") + std::string("Message : ") +
#endif
            msg;

        switch(type)
        {
            case efatal:
                if( outStream )
                {
                    (*outStream) << "Fatal   : " << baseMsg << std::endl;
                }
                else
                {
                    std::cerr << std::endl << "Fatal   : " << baseMsg << std::endl;
                }
                throw NekError(baseMsg);
                break;

            case ewarning:
                if( outStream )
                {
                    (*outStream) << "Warning: " << baseMsg << std::endl;
                }
                else
                {
                    std::cerr << "Warning: " << baseMsg << std::endl;
                }
                break;

            default:
                std::cerr << "Unknown warning type: " << baseMsg << std::endl;
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
} // end of namespace

/// Assert Level 0 -- Fundamental assert which
/// is used whether in FULLDEBUG, DEBUG or OPT
/// compilation mode.  This level assert is
/// considered code critical, even under
/// optimized compilation.

#define NEKERROR(type, msg) \
    ErrorUtil::Error(type, __FILE__, __LINE__, msg, 0);

#define ASSERTL0(condition,msg) \
    if(!(condition)) \
{ \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg, 0); \
}


/// Assert Level 1 -- Debugging which is used whether in FULLDEBUG or
/// DEBUG compilation mode.  This level assert is designed for aiding
/// in standard debug (-g) mode
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)

#define ASSERTL1(condition,msg) \
    if(!(condition)) \
{ \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg, 1); \
}

#else //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
#define ASSERTL1(condition,msg)
#endif //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)


/// Assert Level 2 -- Debugging which is used FULLDEBUG compilation
/// mode.  This level assert is designed to provide addition safety
/// checks within the code (such as bounds checking, etc.).
#ifdef NEKTAR_FULLDEBUG

#define ASSERTL2(condition,msg) \
    if(!(condition)) \
{ \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg, 2); \
}

#else //NEKTAR_FULLDEBUG
#define ASSERTL2(condition,msg)
#endif //NEKTAR_FULLDEBUG

#endif //ERRORUTIL_HPP

