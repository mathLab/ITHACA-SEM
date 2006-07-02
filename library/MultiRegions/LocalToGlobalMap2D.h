///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap2D.h
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
// Description: Local to Global mapping routines in 2D, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP2D_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP2D_H

#include <MultiRegions/LocalToGlobalMap.h>
#include <SpatialDomains/MeshGraph2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	class LocalToGlobalMap2D: 
	public LocalToGlobalMap
	{
        public:
            LocalToGlobalMap2D(){};
	    LocalToGlobalMap2D(int loclen, std::vector<StdRegions::StdExpansionVector> &exp_shapes, SpatialDomains::MeshGraph2D &graph2D);

            virtual ~LocalToGlobalMap2D();
	    
	    inline double GetSign(int i) 
	    {
		if(m_sign_change)
		{
		    return m_sign.get()[i];
		}
		else
		{
		    return 1.0;
		}
	    }
	    
	    inline void LocalToCont(double *loc, double *cont)
	    {
		if(m_sign_change)
		{
		    Vmath::Vmul(m_totLocLen,m_sign.get(),1,loc,1,loc,1);
		}
		
		Vmath::Scatr(m_totLocLen,loc,m_locToContMap.get(),cont);
		
		// additional sign change to make sure local vector is unaffected
		if(m_sign_change)
		{
		    Vmath::Vmul(m_totLocLen,m_sign.get(),1,loc,1,loc,1);
		}
	    }
	    
	    inline void ContToLocal(double *cont, double *loc)
	    {
		Vmath::Gathr(m_totLocLen,cont,m_locToContMap.get(),loc);
		
		if(m_sign_change)
		{
		    Vmath::Vmul(m_totLocLen,m_sign.get(),1,loc,1,loc,1);
		}
	    }
	    
	    inline void Assemble(double *loc, double *cont)
	    {
		Vmath::Zero(m_totGloLen,cont,1);
		if(m_sign_change)
		{
		    Vmath::Vmul(m_totLocLen,m_sign.get(),1,loc,1,loc,1);
		}

		Vmath::Assmb(m_totLocLen,loc,m_locToContMap.get(),cont);

		// additional sign change to make sure local vector is unaffected
		if(m_sign_change)
		{
		    Vmath::Vmul(m_totLocLen,m_sign.get(),1,loc,1,loc,1);
		}
	    }


        protected:
	    
        private:
	    boost::shared_ptr<double> m_sign;
	    bool  m_sign_change;
	};
	
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP2D_H


/** $Log: Loc2GloMap.h,v $
/** Revision 1.3  2006/06/05 00:14:33  bnelson
/** Fixed a compiler error (couldn't find boost::shared_ptr) and a couple of formatting updates for the standard.
/** */

