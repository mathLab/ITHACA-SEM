
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

    double Timer::TimePerTest(unsigned int n)
    {
        #ifdef _WIN32
            CounterType frequency;
            QueryPerformanceFrequency(&frequency);
            return Elapsed().QuadPart/static_cast<double>(n) * 1.0/frequency.QuadPart;
        #elif defined(__APPLE__)
            CounterType elapsed = Elapsed();
            double result = elapsed.tv_sec/static_cast<double>(n) +
                ( elapsed.tv_usec/static_cast<double>(n) * 1.0e-6);
            return result;
        #else
            CounterType elapsed = Elapsed();
            double result = elapsed.tv_sec/static_cast<double>(n) +
                ( elapsed.tv_nsec/static_cast<double>(n) * 1.0e-9);
            return result;
        #endif
    }

}     
