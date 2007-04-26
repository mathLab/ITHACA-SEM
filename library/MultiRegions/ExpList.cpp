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
	    m_ncoeffs(0),
	    m_npoints(0),
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
	NekDouble ExpList::Integral(ConstNekDoubleSharedArray inarray)
	{
	    int    i;
	    int    cnt = 0;
	    NekDouble sum = 0.0;

	    for(i = 0; i < GetExpSize(); ++i)
	    {
		sum += GetExp(i)->Integral(inarray + cnt);
		cnt += GetExp(i)->GetTotPoints();
	    }
            
	    return sum; 
	}
		
	void ExpList::IProductWRTBase(ConstNekDoubleSharedArray inarray, 
				      NekDoubleSharedArray &outarray)
	{
	    int    i;
	    int    cnt  = 0;
	    int    cnt1 = 0;
            NekDoubleSharedArray elmt_outarray;

	    for(i = 0; i < GetExpSize(); ++i)
	    {
                elmt_outarray = outarray + cnt1;
                GetExp(i)->IProductWRTBase(inarray+cnt,elmt_outarray);
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
		GetExp(i)->IProductWRTBase(((ExpList &)S1).GetExp(i)->GetPhys(),
                                           GetExp(i)->UpdateCoeffs());
	    }
	}
	
       void ExpList::IProductWRTBase(const ExpList &S1, 
				     NekDoubleSharedArray &outarray)
       {
           int cnt = 0;
           int i;
           NekDoubleSharedArray elmt_outarray;

           for(i= 0; i < GetExpSize(); ++i)
           {
               elmt_outarray = outarray + cnt; 
               GetExp(i)->IProductWRTBase(((ExpList &)S1).GetExp(i)->GetPhys(),
                                          elmt_outarray);
               cnt += GetExp(i)->GetNcoeffs();
           }
       }
	
        void ExpList::PhysDeriv(ConstNekDoubleSharedArray inarray,
                                NekDoubleSharedArray &out_d0, 
                                NekDoubleSharedArray &out_d1, 
                                NekDoubleSharedArray &out_d2)
	{
	    int  cnt = 0;
	    int  i;
            NekDoubleSharedArray e_out_d0;
            NekDoubleSharedArray e_out_d1;
            NekDoubleSharedArray e_out_d2;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
                e_out_d0 = out_d0 + cnt;
                if(out_d1)
                {
                    e_out_d1 = out_d1 + cnt;
                }
                
                if(out_d2)
                {
                    e_out_d2 = out_d2 + cnt;
                }
                
		GetExp(i)->PhysDeriv(inarray,e_out_d0,e_out_d1,e_out_d2);
		cnt  += GetExp(i)->GetTotPoints();
	    }
	}

      
	void ExpList::FwdTrans(ConstNekDoubleSharedArray inarray, 
                               NekDoubleSharedArray    &outarray)
	{
	    int cnt  = 0;
	    int cnt1 = 0;
	    int i;
            NekDoubleSharedArray e_outarray;
            
	    for(i= 0; i < GetExpSize(); ++i)
	    {
                e_outarray = outarray + cnt1;
                GetExp(i)->FwdTrans(inarray+cnt, e_outarray);
		cnt  += GetExp(i)->GetTotPoints();
		cnt1 += GetExp(i)->GetNcoeffs();
	    }
	    
	    m_transState = eLocal;
	}

	void ExpList::FwdTrans(ConstNekDoubleSharedArray inarray, 
                               ExpList &Sout)
	{
            FwdTrans(inarray,Sout.UpdateCoeffs());
	}


	void ExpList::FwdTrans(const ExpList &Sin)
	{
	    int cnt  = 0;
	    int i;
            NekDoubleSharedArray e_coeffs;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
                e_coeffs = m_coeffs + cnt;
		GetExp(i)->FwdTrans( ((ExpList &)Sin).GetExp(i)->GetPhys(),
				     e_coeffs);
		cnt += GetExp(i)->GetNcoeffs();
	    }
	    
	    m_transState = eLocal;
	}
		
	void ExpList::BwdTrans(ConstNekDoubleSharedArray inarray,
			       NekDoubleSharedArray &outarray)
	{
	    int  i;
	    int  cnt  = 0;
	    int  cnt1 = 0;
	    NekDoubleSharedArray e_outarray;

	    for(i= 0; i < GetExpSize(); ++i)
	    {
                e_outarray = outarray + cnt1;
		GetExp(i)->BwdTrans(inarray + cnt, e_outarray);
		cnt   += GetExp(i)->GetNcoeffs();
		cnt1  += GetExp(i)->GetTotPoints();
	    }
	    
	    m_physState = true;
	}

	void ExpList::BwdTrans(const ExpList &Sin)
	{
	    int  i;
	    int  cnt  = 0;
            ConstNekDoubleSharedArray inarray = Sin.GetCoeffs();
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
		GetExp(i)->BwdTrans(inarray + cnt,GetExp(i)->UpdatePhys());
		cnt   += GetExp(i)->GetNcoeffs();
	    }
	    
	    m_physState = true;
	}
	
        void ExpList::GetCoords(NekDoubleArrayVector &coords)
        {
	    switch(GetExp(0)->GetCoordim())
            {
            case 1:
                GetCoords(coords[0]);
                break;
            case 2:
                GetCoords(coords[0],coords[1]);
                break;
            case 3:
                GetCoords(coords[0],coords[1],coords[2]);
                break;
            }
        }

	void ExpList::GetCoords(NekDoubleSharedArray &coord_0,
				NekDoubleSharedArray &coord_1,
				NekDoubleSharedArray &coord_2)
	{
	    int    i, j, cnt = 0;
	    NekDoubleSharedArray e_coord_0;
	    NekDoubleSharedArray e_coord_1;
	    NekDoubleSharedArray e_coord_2;

	    switch(GetExp(0)->GetCoordim())
	    {
	    case 1:
		for(i= 0; i < GetExpSize(); ++i)
		{
                    e_coord_0 = coord_0 + cnt;
		    GetExp(i)->GetCoords(e_coord_0);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;
	    case 2: 
		ASSERTL0(coord_1 != NullNekDoubleSharedArray, 
			 "output coord_1 is not defined");
    
		for(i= 0; i < GetExpSize(); ++i)
		{
                    e_coord_0 = coord_0 + cnt;
                    e_coord_1 = coord_1 + cnt;
		    GetExp(i)->GetCoords(e_coord_0,e_coord_1);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;
	    case 3: 
		ASSERTL0(coord_1 != NullNekDoubleSharedArray, 
			 "output coord_1 is not defined");
		ASSERTL0(coord_2 != NullNekDoubleSharedArray, 
			 "output coord_2 is not defined");
    
		for(i= 0; i < GetExpSize(); ++i)
		{
                    e_coord_0 = coord_0 + cnt;
                    e_coord_1 = coord_1 + cnt;
                    e_coord_2 = coord_2 + cnt;

		    GetExp(i)->GetCoords(e_coord_0,e_coord_1,e_coord_2);
		    cnt  += GetExp(i)->GetTotPoints();
		}
		break;
            }
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
	
	NekDouble  ExpList::Linf(ConstNekDoubleSharedArray sol)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    NekDouble err = 0.0;
	    int       i,cnt = 0;
	    
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
	
	NekDouble  ExpList::L2(ConstNekDoubleSharedArray sol)
	{
	    NekDouble err = 0.0,errl2;
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

