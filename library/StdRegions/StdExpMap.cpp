///////////////////////////////////////////////////////////////////////////////
//
// File StdExpMap.cpp
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
// Description: Determine mapping information between 1D and 2D
// expansions as well as between 2D and 3D Expansions
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpMap.h>

namespace Nektar
{
    namespace StdRegions
    {
	
	StdExpMap::StdExpMap():
	    m_len(0),
	    m_map(NULL)
	{
	}
	
	StdExpMap::StdExpMap(const int len)
	{
	    m_len = len;
	    ASSERTL2(len > 0,"called with zero length");
	    m_map = MemoryManager::AllocateSharedArray<int>(m_len);
	}
	
	
	StdExpMap::~StdExpMap()
	{
	}
	
	void StdExpMap::SetMapMemory(const int len)
	{	    
	    if(m_len != len)
	    {
		//MemoryManager::DeallocateArray<int> (m_map,m_len);
		m_len = len;
		m_map = MemoryManager::AllocateSharedArray<int> (m_len);
	    }
	}
	
    } // end of namespace
} // end of namespace

/** 
 * $Log: StdExpMap.cpp,v $
 * Revision 1.2  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.4  2006/04/01 21:59:26  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.3  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 *
 **/ 
