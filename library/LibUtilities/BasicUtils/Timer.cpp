///////////////////////////////////////////////////////////////////////////////
//
// File: Timer.cpp
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

#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
    Timer::Timer() :
        m_start(),
        m_end(),
        m_resolution()
    {
    }

    Timer::~Timer()
    {
    }

    void Timer::Start()
    {
        #ifdef _WIN32
            QueryPerformanceCounter(&m_start);
        #elif defined(__APPLE__)
            gettimeofday(&m_start, 0);
        #else 
            clock_gettime(CLOCK_REALTIME, &m_start);
        #endif
    }

    void Timer::Stop()
    {
        #ifdef _WIN32
            QueryPerformanceCounter(&m_end);
        #elif defined(__APPLE__)
            gettimeofday(&m_end, 0);
        #else 
            clock_gettime(CLOCK_REALTIME, &m_end);
        #endif
    }

    Timer::CounterType Timer::Elapsed()
    {
        #ifdef _WIN32
            CounterType result;
            result.QuadPart = m_end.QuadPart - m_start.QuadPart;
            return result;
        #elif defined(__APPLE__)
            CounterType result = m_end;
            
            if( result.tv_usec < m_start.tv_usec) 
            {
                result.tv_sec -= 1;
                result.tv_usec += 1000000;
            }
            
            result.tv_sec -= m_start.tv_sec;
            result.tv_usec -= m_start.tv_usec;
            
            return result;
        #else
            CounterType result = m_end;

            if( result.tv_nsec < m_start.tv_nsec) 
            {
                result.tv_sec -= 1;
                result.tv_nsec += 1000000000;
            }
            
            result.tv_sec -= m_start.tv_sec;
            result.tv_nsec -= m_start.tv_nsec;

            return result;
        #endif
    }

    NekDouble Timer::TimePerTest(unsigned int n)
    {
        #ifdef _WIN32
            CounterType frequency;
            QueryPerformanceFrequency(&frequency);
            return Elapsed().QuadPart/static_cast<NekDouble>(n) * 1.0/frequency.QuadPart;
        #elif defined(__APPLE__)
            CounterType elapsed = Elapsed();
            NekDouble result = elapsed.tv_sec/static_cast<NekDouble>(n) +
                ( elapsed.tv_usec/static_cast<NekDouble>(n) * 1.0e-6);
            return result;
        #else
            CounterType elapsed = Elapsed();
            NekDouble result = elapsed.tv_sec/static_cast<NekDouble>(n) +
                ( elapsed.tv_nsec/static_cast<NekDouble>(n) * 1.0e-9);
            return result;
        #endif
    }

}     
