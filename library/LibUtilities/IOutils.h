///////////////////////////////////////////////////////////////////////////////
//
// File IOUtils.h
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
// Description: Header file for Input/Output Utilities
//
///////////////////////////////////////////////////////////////////////////////

#ifndef IOUTILS_H
#define IOUTILS_H

#include <string>
#include <iostream>
#include <fstream>

namespace IOutils
{

    int FindToken(const std::string& line, const std::string& identifier);

    std::string FindSection(std::ifstream& in, const std::string& token);

}

#endif //IOUTILS_H

/***
$Log: IOutils.h,v $
Revision 1.1  2006/05/04 18:57:42  kirby
*** empty log message ***

Revision 1.5  2006/02/12 21:51:42  sherwin

Added licence

Revision 1.4  2006/02/12 15:06:12  sherwin

Changed .h files to .hpp

**/
