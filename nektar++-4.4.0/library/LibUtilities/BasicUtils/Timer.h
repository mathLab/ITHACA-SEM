///////////////////////////////////////////////////////////////////////////////
//
// File: Timer.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Time getting class
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_TIMER_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_TIMER_H

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

namespace Nektar
{
    class Timer
    {    
        public:
            #ifdef _WIN32
                typedef LARGE_INTEGER CounterType;
            #elif  defined(__APPLE__)
                typedef timeval CounterType;
            #else
                typedef timespec CounterType;
            #endif

        public:
            LIB_UTILITIES_EXPORT Timer();
            LIB_UTILITIES_EXPORT ~Timer();

            LIB_UTILITIES_EXPORT void Start();
            LIB_UTILITIES_EXPORT void Stop();
            LIB_UTILITIES_EXPORT CounterType Elapsed();

            /// \brief Returns amount of seconds per iteration in
            ///        a test with n iterations.
            LIB_UTILITIES_EXPORT NekDouble TimePerTest(unsigned int n);

        private:
            Timer(const Timer& rhs);
            Timer& operator=(const Timer& rhs);

            CounterType m_start;
            CounterType m_end;
            CounterType m_resolution;
    };
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_TIMER_H
