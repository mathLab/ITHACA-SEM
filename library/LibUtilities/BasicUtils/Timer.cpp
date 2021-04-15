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
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <tuple>

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

void Timer::AccumulateRegion(std::string region, int iolevel)
{
    // search for region
    auto search = m_elapsedRegion.find(region);
    if (search == m_elapsedRegion.end())
    {
        m_elapsedRegion.insert({region,
         std::make_tuple<Timer::Seconds, size_t>(this->Elapsed(),1, iolevel)});
        // update field width
        m_maxStringWidth = std::max(
            static_cast<decltype(m_maxStringWidth)>(region.size()),
            m_maxStringWidth);
    }
    else
    {
        std::get<0>(search->second) += this->Elapsed();
        std::get<1>(search->second) += 1;
    }
}

void Timer::PrintElapsedRegions(LibUtilities::CommSharedPtr comm,
                                std::ostream &o,
                                int iolevel)
{
    std::vector<std::string> labels{
        "Region",
        "Elapsed time Avg (s)",
        "Min (s)",
        "Max (s)",
        "Count"};


    if (comm->GetRank() == 0 &&
        m_elapsedRegion.begin() != m_elapsedRegion.end())
    {
        o << "-------------------------------------------\n";
        o << std::setw(m_maxStringWidth+2) << labels[0] << '\t'
          << std::setw(10) << labels[1] << '\t'
          << std::setw(10) << labels[2] << '\t'
          << std::setw(10) << labels[3] << '\t'
          << std::setw(10) << labels[4] << '\n';
    }
    // first write out execute time
    auto item = m_elapsedRegion.find("Execute");
    if(item != m_elapsedRegion.end())
    {
        auto elapsedAve = std::get<0>(item->second).count();
        comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
        elapsedAve /= comm->GetSize();
        auto elapsedMin = std::get<0>(item->second).count();
        comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
        auto elapsedMax = std::get<0>(item->second).count();
        comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);
        
        if (comm->GetRank() == 0)
        {
            o << std::setw(m_maxStringWidth+2) << item->first << '\t'
              << std::setw(10) << elapsedAve << '\t'
              << std::setw(10) << elapsedMin << '\t'
              << std::setw(10) << elapsedMax << '\t'
              << std::setw(10) << std::get<1>(item->second) << '\n';
        }
    }            

    // write out all other timings order alphabetically on string
    for (auto item = m_elapsedRegion.begin();
            item != m_elapsedRegion.end(); ++item)
    {
        if(std::get<2>(item->second) < iolevel)
        {
            if(boost::iequals(item->first,"Execute"))
            {
                continue;
            }

            auto elapsedAve = std::get<0>(item->second).count();
            comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
            elapsedAve /= comm->GetSize();
            auto elapsedMin = std::get<0>(item->second).count();
            comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
            auto elapsedMax = std::get<0>(item->second).count();
            comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);

            if (comm->GetRank() == 0)
            {
                o << std::setw(m_maxStringWidth+2) << item->first << '\t'
                  << std::setw(10) << elapsedAve << '\t'
                  << std::setw(10) << elapsedMin << '\t'
                  << std::setw(10) << elapsedMax << '\t'
                  << std::setw(10) << std::get<1>(item->second) << '\n';
            }
        }
    }
}
// static members init
std::map<std::string, std::tuple<Timer::Seconds, size_t, int>>
    Timer::m_elapsedRegion{};

unsigned short Timer::m_maxStringWidth = 10;

}
} // end Nektar namespace
