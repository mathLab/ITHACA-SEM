
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
        #ifdef _WIN32
            QueryPerformanceFrequency(&m_resolution);
        #else
            clock_getres(CLOCK_REALTIME, &m_resolution);
        #endif
    }

    void Timer::Start()
    {
        #ifdef _WIN32
            QueryPerformanceCounter(&m_start);
        #else
            clock_gettime(CLOCK_REALTIME, &m_start);
        #endif
    }

    void Timer::Stop()
    {
        #ifdef _WIN32
            QueryPerformanceCounter(&m_end);
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
        #else
            CounterType elapsed = Elapsed();
            double result = elapsed.tv_sec/static_cast<double>(n) +
                ( elapsed.tv_nsec/static_cast<double>(n) * 1.0e-9);
            return result;
        #endif
    }

}     
