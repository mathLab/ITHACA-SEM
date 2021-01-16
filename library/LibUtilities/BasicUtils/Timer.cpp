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
#include <LibUtilities/Communication/CommSerial.h>

#include <iomanip>

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
    // search for region
    auto search = m_elapsedRegion.find(region);
    if (search == m_elapsedRegion.end())
    {
        m_elapsedRegion.insert({region,
            std::make_pair<Timer::Seconds, size_t>(this->Elapsed(),1)});
        // update width field width
        m_maxStringWidth = std::max(
            static_cast<decltype(m_maxStringWidth)>(region.size()),
            m_maxStringWidth);
    }
    else
    {
        search->second.first += this->Elapsed();
        search->second.second += 1;
    }
}

void Timer::PrintElapsedRegions()
{
    std::string  def("default");
    char *argv = new char [def.length()+1];
    std::strcpy(argv,def.c_str());
    LibUtilities::CommSharedPtr comm = 
        MemoryManager<LibUtilities::CommSerial>:: AllocateSharedPtr(1,&argv);

    PrintElapsedRegions(comm);
}
    
void Timer::PrintElapsedRegions(LibUtilities::CommSharedPtr comm,
                                std::ostream &o)
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
    for (auto item = m_elapsedRegion.begin();
            item != m_elapsedRegion.end(); ++item)
    {
        auto elapsedAve = item->second.first.count();
        comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
        elapsedAve /= comm->GetSize();
        auto elapsedMin = item->second.first.count();
        comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
        auto elapsedMax = item->second.first.count();
        comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);

        if (comm->GetRank() == 0)
        {
            o << std::setw(m_maxStringWidth+2) << item->first << '\t'
              << std::setw(10) << elapsedAve << '\t'
              << std::setw(10) << elapsedMin << '\t'
              << std::setw(10) << elapsedMax << '\t'
              << std::setw(10) << item->second.second << '\n';
        }
    }
}

// static members init
std::unordered_map<std::string, std::pair<Timer::Seconds, size_t>>
    Timer::m_elapsedRegion{};

unsigned short Timer::m_maxStringWidth = 10;

}
} // end Nektar namespace
