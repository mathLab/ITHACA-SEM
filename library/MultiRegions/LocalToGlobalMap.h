///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap.h
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

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
	
	class LocalToGlobalMap
	{
        public:
            LocalToGlobalMap();
	    LocalToGlobalMap(const int totdata, int *map);
	    
            ~LocalToGlobalMap();
	    
	    inline int GetMap(const int n, const int i) const
	    {
		return m_locToContMap[n][i];
	    }

	    inline void LocalToCont(ExpList &loc, double *cont)
	    {
		int i;
		for(i = 0; i < loc.GetNexp(); ++i)
		{
		    Vmath::Scatr(loc.GetExp(i).GetNcoeffs(), loc.GetExp(i).GetCoeffs(),
				 &m_locToContMap[i][0],cont);
		}
	    }
	    
	    inline void ContToLocal(const double *cont, ExpList &loc)
	    {
		int i;
		for(i = 0; i < loc.GetNexp(); ++i)
		{
		    Vmath::Gathr(loc.GetExp(i).GetNCoeffs,cont,&m_locToContMap[i][0],
				 loc.GetExp(i).GetCoeffs());
		}

	    }
	    
	    inline void Assemble(ExpList &loc, double *cont)
	    {
		int i;

		Vmath::Zero(m_totGloLen,cont,1);
		for(i = 0; i < loc.GetNexp(); ++i)
		{
		    Vmath::Assmb(loc.GetExp(i).GetNcoeffs(),loc.GetExp(i).GetCoeffs(),
				 &m_locToContMap[i][0],cont);
		}
	    }
	    
	    inline int GetTotGloLen()
	    {
		return m_totGloLen;
	    }

        protected:
            int m_totLocLen;    //< length of local mapping 
	    int m_totGloLen;    //< length of global dofs
            boost::shared_ptr<int> m_locToContMap;
        private:
	};
	
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** $Log: LocalToGlobalMap.h,v $
/** Revision 1.1  2006/07/02 17:16:17  sherwin
/**
/** Modifications to make MultiRegions work for a connected domain in 2D (Tris)
/**
/** Revision 1.3  2006/06/05 00:14:33  bnelson
/** Fixed a compiler error (couldn't find boost::shared_ptr) and a couple of formatting updates for the standard.
/** */

