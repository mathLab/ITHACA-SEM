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

namespace ErrorUtil
{

    enum ErrType
    {
        efatal,
        ewarning
    };

#ifdef ENABLE_NEKTAR_EXCEPTIONS
    static void ExceptionError(ErrType type, const char *routine, int lineNumber, const char *msg)
    {
        std::string errorMessage = std::string(routine) + "[" +
                boost::lexical_cast<std::string>(lineNumber) + "]:" +
                std::string(msg);
        std::cerr << errorMessage << std::endl;
        throw std::runtime_error(errorMessage);
    }
#endif

    static void Error(ErrType type, const char *routine, int lineNumber, const char *msg)
    {
        switch(type)
        {
            case efatal:
                std::cerr << routine << "[" << lineNumber << "]:" << msg << std::endl;
                exit(1);
                break;
            case ewarning:
                std::cerr << routine << ": " << msg << std::endl;
                break;
            default:
                std::cerr << "Unknown warning type" << std::endl;
        }
    }
} // end of namespace

/// Assert Level 0 -- Fundamental assert which
/// is used whether in FULLDEBUG, DEBUG or OPT
/// compilation mode.  This level assert is
/// considered code critical, even under
/// optimized compilation.

#define NEKERROR(type, msg) \
    ErrorUtil::Error(type, __FILE__, __LINE__, msg);

#ifdef ENABLE_NEKTAR_EXCEPTIONS

#define ASSERTL0(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 0 Assert Violation\n"); \
    ErrorUtil::ExceptionError(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#else

#define ASSERTL0(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 0 Assert Violation\n"); \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#endif

/// Assert Level 1 -- Debugging which is used whether in FULLDEBUG or
/// DEBUG compilation mode.  This level assert is designed for aiding
/// in standard debug (-g) mode
#ifdef DEBUG
#ifdef ENABLE_NEKTAR_EXCEPTIONS

#define ASSERTL1(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 1 Assert Violation\n"); \
    ErrorUtil::ExceptionError(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#else

#define ASSERTL1(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 1 Assert Violation\n"); \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#endif

#else
#define ASSERTL1(condition,msg)
#endif


/// Assert Level 2 -- Debugging which is used FULLDEBUG compilation
/// mode.  This level assert is designed to provide addition safety
/// checks within the code (such as bounds checking, etc.).
#ifdef FULLDEBUG

#ifdef ENABLE_NEKTAR_EXCEPTIONS

#define ASSERTL2(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 2 Assert Violation\n"); \
    ErrorUtil::ExceptionError(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#else

#define ASSERTL2(condition,msg) \
    if(!(condition)) \
{ \
    fprintf(stderr,"Level 2 Assert Violation\n"); \
    ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__, msg); \
}

#endif

#else
#define ASSERTL2(condition,msg)
#endif

#endif //ERRORUTIL_HPP
/***
$Log: ErrorUtil.hpp,v $
Revision 1.1  2006/06/01 11:07:52  kirby
*** empty log message ***

Revision 1.4  2006/05/16 20:40:44  jfrazier
Renamed ERROR to NEKERROR to prevent conflict in Visual Studio.

Revision 1.3  2006/05/09 16:40:57  jfrazier
Added ERROR macro definition.

Revision 1.2  2006/05/07 18:51:05  bnelson
Format changes for coding standard.

Revision 1.1  2006/05/04 18:57:41  kirby
*** empty log message ***

Revision 1.7  2006/04/14 14:51:17  jfrazier
Fixed a problem which is most likely a preprocessor problem.  The file and line
number were inconsistent between release and debug builds.

**/
