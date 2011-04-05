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
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace ErrorUtil
{
    LIB_UTILITIES_EXPORT void SetErrorStream(std::ostream& o);
    LIB_UTILITIES_EXPORT bool HasCustomErrorStream();

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
        
    LIB_UTILITIES_EXPORT void Error(ErrType type, const char *routine, int lineNumber, const char *msg, unsigned int level);
    LIB_UTILITIES_EXPORT void Error(ErrType type, const char *routine, int lineNumber, const std::string& msg, unsigned int level);
    LIB_UTILITIES_EXPORT void Error(ErrType type, const char *routine, int lineNumber, const char *msg);
    
    
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

/***
$Log: ErrorUtil.hpp,v $
Revision 1.9  2008/01/20 04:02:49  bnelson
Fixed compiler warnings.

Revision 1.8  2007/07/10 22:22:19  jfrazier
Added an additional indicator of the error severity (fatal or warning).

Revision 1.7  2007/06/10 23:36:59  bnelson
A previous change to ErrorUtil::Error added an additional level parameter, breaking code which called it directly outside of ErrorUtil.  A new method without this parameter was added.

Revision 1.6  2007/06/06 04:57:21  bnelson
Removed all printing to the console when error reporting is done via exceptions.  This makes it easier to see errors in the UnitTests.

Revision 1.5  2007/05/27 18:28:23  bnelson
*** empty log message ***

Revision 1.4  2007/05/22 02:02:35  bnelson
Changed Array::size to Array::num_elements.

Fixed some compiler errors in assertions.

Revision 1.3  2007/05/14 23:22:38  bnelson
Errors throw a new exception type NekError.

Revision 1.2  2006/08/14 02:18:02  bnelson
Added the option to throw exceptions when an error is encountered.

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
