///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.cpp
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
// Description: Expansion list definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>

namespace Nektar
{
    namespace MultiRegions
    {

	ExpList::ExpList(void):
	    m_transState(eNotSet),
	    m_physState(false)
	    
	{
	}
	
	ExpList::~ExpList()
	{
	}
	
	
	/** \brief Integrate the physical point list \a inarray over region
	    and return the value
	    
	    Inputs:\n
	    
	    - \a inarray: definition of function to be returned at quadrature point 
	    of expansion. 
	    
	    Outputs:\n
	    
	    - returns \f$ \sum_{i=1}^{n_{el}} \int_{\Omega_i} u(\xi_1)d \xi_1 \f$ 
	*/
	double ExpList::Integral(const NekDoubleSharedArray &inarray)
	{
	    int    i;
	    int    cnt = 0;
	    double sum = 0.0;
	    NekDoubleSharedArray offset;
	    
	    for(i = 0; i < GetExpSize(); ++i)
	    {
		sum += GetExp(i)->Integral(offset.reset(inarray[cnt]));
		cnt += GetExp(i)->GetTotPoints();
	    }

	    return sum; 
	}
	
	
	void ExpList::IProductWRTBase(const NekDoubleSharedArray &inarray, 
				      NekDoubleSharedArray &outarray)
	{
	    int    i;
	    int    cnt  = 0;
	    int    cnt1 = 0;
	  
	    for(i = 0; i < GetExpSize(); ++i)
	    {
		    GetExp(i)->IProductWRTBase(inarray+cnt,outarray+cnt1);
		    cnt  += GetExp(i)->GetTotPoints();
		    cnt1 += GetExp(i)->GetNcoeffs();
	    }
	    m_transState = eLocal;
	}
	
	void ExpList::IProductWRTBase(const ExpList &S1)
	{
	    int i;
	    
	    for(i = 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->IProductWRTBase(((ExpList) S1).GetExp(i)->GetPhys(),
					   GetExp(i)->GetCoeffs());
	    }
	}
	
       void ExpList::IProductWRTBase(const ExpList &S1, 
				     NekDoubleSharedArray &outarray)
	{
	    int cnt = 0;
	    int i;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->IProductWRTBase(&(((ExpList)S1).
					     GetExp(i)->GetPhys())[0],
					   outarray+cnt);
		cnt += GetExp(i)->GetNcoeffs();
	    }
	}
	
	void ExpList::PhysDeriv(NekDoubleSharedArray out_0,
				NekDoubleSharedArray out_1,
				NekDoubleSharedArray out_2)
	{
	    int  cnt = 0;
	    int  i;

	    
	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    

	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->PhysDeriv(&GetExp(i)->GetPhys()[0],
				     out_0+cnt,out_1+cnt,out_2+cnt);
		cnt  += GetExp(i)->GetTotPoints();
	    }
	}

      
	void ExpList::FwdTrans(const NekDoubleSharedArray &inarray, NekDoubleSharedArray &outarray)
	{
	    int cnt  = 0;
	    int cnt1 = 0;
	    int i;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->FwdTrans(inarray+cnt, outarray + cnt1);
		cnt  += GetExp(i)->GetTotPoints();
		cnt1 += GetExp(i)->GetNcoeffs();
	    }
	    
	    m_transState = eLocal;
	}

	void ExpList::FwdTrans(const NekDoubleSharedArray &inarray, ExpList &Sout)
	{
	    int cnt  = 0;
	    int i;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->FwdTrans(inarray+cnt, Sout.GetExp(i)->GetCoeffs());
		cnt  += GetExp(i)->GetTotPoints();
	    }
	    
	    m_transState = eLocal;
	}


	void ExpList::FwdTrans(const ExpList &Sin)
	{
	    int i;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->FwdTrans( ((ExpList )Sin).GetExp(i)->GetPhys(),
				     GetExp(i)->GetCoeffs());
	    }
	    
	    m_transState = eLocal;
	}
		
	void ExpList::BwdTrans(const NekDoubleSharedArray &inarray,
			       NekDoubleSharedArray &outarray)
	{
	    int  i;
	    int  cnt = 0;
	    int  cnt1 = 0;
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->BwdTrans(inarray + cnt, outarray+cnt1);
		cnt   += GetExp(i)->GetNcoeffs();
		cnt1  += GetExp(i)->GetTotPoints();
	    }
	    
	    m_physState = true;
	}

	void ExpList::BwdTrans(const ExpList &Sin)
	{
	    int  i;
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->BwdTrans(((ExpList) Sin).GetExp(i)->GetCoeffs(),
				    GetExp(i)->GetPhys());
	    }
	    
	    m_physState = true;
	}
	
	void ExpList::GetCoords(NekDoubleSharedArray &coord_0,
				NekDoubleSharedArray &coord_1,
				NekDoubleSharedArray &coord_2)
	{
	    int    i, j, cnt = 0;
	    NekDoubleSharedArray &E_coords[3];
	    
	    
	    switch(GetExp(0)->GetCoordim())
	    {
	    case 1:
		for(i= 0; i < GetExpSize(); ++i)
		{
		    
		    GetExp(i)->GetCoords(coord_0+cnt);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;
	    case 2: 
		ASSERTL0(coord_1 != NullNekDoubleSharedArray, 
			 "output coord_1 is not defined");
    
		for(i= 0; i < GetExpSize(); ++i)
		{
		    GetExp(i)->GetCoords(coord_0+cnt,coord_1+cnt);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;

	    case 2: 
		ASSERTL0(coord_1 != NullNekDoubleSharedArray, 
			 "output coord_1 is not defined");
		ASSERTL0(coord_2 != NullNekDoubleSharedArray, 
			 "output coord_2 is not defined");
    
		for(i= 0; i < GetExpSize(); ++i)
		{
		    GetExp(i)->GetCoords(coord_0+cnt,coord_1+cnt,coord_2+cnt);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;

	}
	
	void ExpList::WriteToFile(std::ofstream &out)
	{
	    int i;
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    
	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    
	    m_exp[0]->WriteToFile(out,1);
	    
	    for(i= 1; i < GetExpSize(); ++i)
	    {
		GetExp(i)->WriteToFile(out,0);
	    }
	}
	
	double  ExpList::Linf(const NekDoubleSharedArray &sol)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    double err = 0.0;
	    int    i,cnt = 0;
	    
	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
		err  = std::max(err,GetExp(i)->Linf(sol+cnt));
		cnt  += GetExp(i)->GetTotPoints();
	    }
	    
	    return err;
	}
	
	double  ExpList::L2(const NekDoubleSharedArray &sol)
	{
	    double err = 0.0,errl2;
	    int    i,cnt = 0;
	    
	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
		errl2 = GetExp(i)->L2(sol+cnt);
		err += errl2*errl2;
		cnt  += GetExp(i)->GetTotPoints();
	    }
	    
	    return sqrt(err);
	}
	
    } //end of namespace
} //end of namespace

