///////////////////////////////////////////////////////////////////////////////
//
// File IOUtils.cpp
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
// Description: Input/Output Utilities
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/IOutils.h>
#include <LibUtilities/ErrorUtil.hpp>

namespace IOutils
{

  int FindToken(const std::string &line, const std::string &identifier)
  {
    int i, start;
    std::string Token;

    if((start = line.find_first_of("<",0))!=std::string::npos) // found '<'
    {
      for(i =start; i < line.size(); ++i)
      {
    if(line[i] == '>')
    {
      Token += line[i];
      break;
    }

    if((line[i] != '\t')&&(line[i]!=' '))// remove whitespace
    {
      Token += toupper(line[i]);
    }
      }

      if(Token.find(identifier,0) != std::string::npos)
      {
      return true;
      }
    }
    return false;
  }

  /// \brief  search file 'in' for a line with the token 'token'.
  std::string FindSection(std::ifstream &in, const std::string token)
  {
    std::string line;

    while(getline(in,line)) // Find Nodes section
    {
      if(IOutils::FindToken(line,token))
      {
    break;
      }
    }

    // check to see if at end of file or have token
    if(!IOutils::FindToken(line,token))
    {
      in.clear();
      in.seekg(0,std::ios::beg); // rewind file;

      while(getline(in,line)) // Find Nodes section
      {
    if(IOutils::FindToken(line,token))
    {
      break;
    }
      }

      if(!IOutils::FindToken(line,token))
      {
    line = "Cannot find Section "+token;
    ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,line.c_str());
      }
    }
    return line;
  }
}

/***
$Log: IOutils.cpp,v $
Revision 1.1  2006/05/04 18:57:42  kirby
*** empty log message ***

Revision 1.9  2006/03/13 11:17:03  sherwin

First compiing version of Demos in SpatialDomains and LocalRegions. However they do not currently seem to execute properly

Revision 1.8  2006/02/26 21:13:45  bnelson
Fixed a variety of compiler errors caused by updates to the coding standard.

Revision 1.7  2006/02/12 21:51:42  sherwin

Added licence

Revision 1.6  2006/02/12 15:06:12  sherwin

Changed .h files to .hpp

**/
