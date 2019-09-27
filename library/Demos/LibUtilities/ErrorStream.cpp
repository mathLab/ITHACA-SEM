////////////////////////////////////////////////////////////////////////////////
//
//  File: ErrorStream.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Ensure that error stream tester and exception handling works.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <sstream>
#include <vector>

#include <boost/core/ignore_unused.hpp>

using namespace Nektar;

int main(int argc, char *argv[])
{
    boost::ignore_unused(argc, argv);

    // Set up a stringstream to catch any error output.
    std::stringstream ss;
    ErrorUtil::SetErrorStream(ss);
    ErrorUtil::SetPrintBacktrace(false);

    // Set up output that will be overwritten with exception error (any
    // non-empty string).
    std::string output = "some error";

    // Run a piece of code designed to generate an error: create a SessionReader
    // with missing file.
    try
    {
        char *fake_argv[3] = {
            (char *)"ErrorStream", (char *)"missing.xml", NULL };
        LibUtilities::SessionReaderSharedPtr session =
            LibUtilities::SessionReader::CreateInstance(2, fake_argv);
        session->InitSession();
    }
    catch (const ErrorUtil::NekError &e)
    {
        // Format exception error to the appropriate message that should be
        // output by the ErrorUtil class.
        std::stringstream outf;
        outf << "Fatal   : " << e.what() << std::endl;
        output = outf.str();
    }

    // Check to see whether we get the correct output.
    if (output == ss.str())
    {
        std::cout << "Caught exception error message matches stringstream"
                  << std::endl;
        return 0;
    }

    std::cout << "Difference between exception message and stringstream."
              << std::endl << "Stringstream contains:" << std::endl
              << ss.str() << std::endl << std::endl
              << "Exception contains:" << std::endl
              << output << std::endl;
    return 1;
}
