////////////////////////////////////////////////////////////////////////////////
//
// Progressbar.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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
// Description: Print a simple progress bar
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PROGRESSBAR_HPP
#define NEKTAR_LIBUTILITIES_PROGRESSBAR_HPP

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#ifdef _WIN32
#include <io.h>
#define ISTTY _isatty(_fileno(stdout))
#else
#include <unistd.h>
#define ISTTY isatty(fileno(stdout))
#endif

namespace Nektar
{
namespace LibUtilities
{

/**
 * @brief Prints a progressbar
 *
 * @param position  State of the current process
 * @param goal      Goal of the current process
 * @param message   Short Description of the current process
 *
 * This function plots a simple progressbar to the console or log file to
 * visualize the current state of an ongoing process. Make sure you minimize
 * calling this routine. Ideally, this should be called only when the
 * percentage is increased by an integer.
 */
inline int PrintProgressbar(
    const int position, const int goal, const std::string message,
    int lastprogress = -1)
{
    std::cout.unsetf ( std::ios::floatfield );
    if (ISTTY)
    {
        float progress = position / float(goal);
        int  numeq = static_cast<int>(ceil(progress *49));
        if(lastprogress == numeq)
        {
            return numeq;
        }
        else
        {
            // carriage return
            std::cout << "\r";

            std::cout << message << ": ";

            std::cout << std::setw(3) << ceil(100 * progress) << "% [";
            for (int j = 0; j < numeq; j++)
            {
                std::cout << "=";
            }
            for (int j = numeq; j < 49; j++)
            {
                std::cout << " ";
            }
            std::cout << "]" << std::flush;
            return numeq;
        }
    }
    else
    {
        // print only every 2 percent
        if (int(ceil(double(100 * position / goal))) % 2 ==  0)
        {
            std::cout << "." <<  std::flush;
        }
        return -1;
    }
}

}
}

#endif // NEKTAR_LIBUTILITIES_PROGRESSBAR_HPP
