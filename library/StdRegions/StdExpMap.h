///////////////////////////////////////////////////////////////////////////////
//
// File StdExpMap.h
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

#ifndef NEKTAR_LIB_STDREGIONS_STDEXPMAP_H
#define NEKTAR_LIB_STDREGIONS_STDEXPMAP_H

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace StdRegions
    {
    
    class StdExpMap
    {
        
    public:
        StdExpMap();
        StdExpMap(const int len);
        ~StdExpMap();
        
        void SetMapMemory(const int len);

        void SetMap(const int id, const int val)
        {
            ASSERTL1(id < m_len,"id is larger than length of map");
            
            m_map[id] = val;
        }
        
        inline int GetLen()
        {
            return m_len;
        }
        
        inline const Array<OneD, const int>& GetMap() const
        {
            return m_map;
        }

        inline const Array<OneD, const int>& GetSign() const
        {
            return m_sign;
        }


        inline const Array<OneD, int>& UpdateSign() 
        {
            return m_sign;
        }
        
        inline NekDouble GetSign(const int i) 
        {
        
            return (NekDouble) m_sign[i];
        }
        
        int operator[](const int i) const
        {
            return m_map[i];
        }
        
        int& operator[](const int i)
        {
            return m_map[i];
        }
        
    protected:
        
    private:
        
        int m_len;
        Array<OneD, int> m_map;
        Array<OneD, int> m_sign;
    };
    
    } // end of namespace
} // end of namespace

#endif //STDEXPMAP_H

/**
 * $Log: StdExpMap.h,v $
 * Revision 1.10  2008/04/06 06:04:14  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.9  2007/07/20 02:16:52  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.8  2007/05/15 05:18:23  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.7  2007/03/29 19:35:08  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.6  2007/03/20 16:58:42  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.5  2007/01/28 18:34:17  sherwin
 * More modifications to make Demo Project1D compile
 *
 * Revision 1.4  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.3  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.2  2006/06/01 13:43:19  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.6  2006/04/25 20:23:33  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.5  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.4  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 *
 **/


