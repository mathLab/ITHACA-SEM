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
namespace LibUtilities
{

void Timer::Start()
{
    m_start = Clock::now();
}

void Timer::Stop()
{
    m_end = Clock::now();
}

Timer::Seconds Timer::Elapsed()
{
    return std::chrono::duration_cast<Seconds>(m_end - m_start);
}

NekDouble Timer::TimePerTest(unsigned int n)
{
    return Elapsed().count() / static_cast<NekDouble>(n);
}

void Timer::AccumulateRegion(std::string region)
{
    auto search = m_elapsedRegion.find(region);
    if (search == m_elapsedRegion.end())
    {
        m_elapsedRegion.insert({region, this->Elapsed()});
    }
    else
    {
        search->second += this->Elapsed();
    }
}

// Timer::Seconds Timer::ElapsedRegion(std::string region)
// {
//     auto search = m_elapsedRegion.find(region);
//     if (search == m_elapsedRegion.end())
//     {
//         NEKERROR(ErrorUtil::efatal, "Region not registered.");
//     }
//     return search->second;
// }

void Timer::PrintElapsedRegions(LibUtilities::CommSharedPtr comm)
{
    if (comm->GetRank() == 0)
    {
        std::cout
            << "-------------------------------------------\n"
            << "Region\t\t Elapsed time Ave (s)"
            << "\t Min (s)"
            << "\t Max (s)\n";
    }
    for (auto item = m_elapsedRegion.begin();
            item != m_elapsedRegion.end(); ++item)
    {
        auto elapsedAve = item->second.count();
        comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
        elapsedAve /= comm->GetSize();
        auto elapsedMin = item->second.count();
        comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
        auto elapsedMax = item->second.count();
        comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);

        if (comm->GetRank() == 0)
        {
            std::cout << item->first << '\t'
                << elapsedAve << '\t'
                << elapsedMin << '\t'
                << elapsedMax << '\n';
        }
    }
}

std::map<std::string, Timer::Seconds> Timer::m_elapsedRegion{};

}
} // end Nektar namespace