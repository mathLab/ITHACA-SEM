///////////////////////////////////////////////////////////////////////////////
//
// File Loc2GloMap.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Local to Global mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LOC2GLOMAP_H
#define LOC2GLOMAP_H

#include <vector>

namespace Nektar
{
  namespace MultiRegions
  {
  
    typedef boost::shared_ptr<int} IntPtr;
    typedef std::vector< IntPtr > IntVector;
    typedef std::vector< IntPtr >::iterator IntVectorIter;
    
    class Loc2GloMap
    {
    public:
      Loc2GloMap();
      ~Loc2GloMap();

    protected:
      
      
    private:
      int *m_nbndry;
      IntVector m_locid;
      IntVector  m_gmap;
      
    };
    
  } // end of namespace
} // end of namespace

#endif
