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
	    m_ncoeffs(0),
	    m_npoints(0),
	    m_transState(eNotSet),
	    m_physState(false)
	{
            m_exp = MemoryManager::AllocateSharedPtr<StdRegions::StdExpansionVector>();
	}
	

        ExpList::ExpList(const ExpList &in):
            m_ncoeffs(in.m_ncoeffs),
            m_npoints(in.m_npoints),
            m_transState(eNotSet),
            m_physState(false),
            m_exp(in.m_exp)
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
        NekDouble ExpList::PhysIntegral(void)
        {
            ASSERTL2(m_physState == true,
                     "local physical space is not true ");

            return PhysIntegral(m_phys);
        }
        
	NekDouble ExpList::PhysIntegral(ConstNekDoubleSharedArray inarray)
	{
	    int       i;
	    int       cnt = 0;
	    NekDouble sum = 0.0;

	    for(i = 0; i < GetExpSize(); ++i)
	    {
		sum += (*m_exp)[i]->Integral(inarray + cnt);
		cnt += (*m_exp)[i]->GetTotPoints();
	    }
            
	    return sum; 
	}
		
	void ExpList::IProductWRTBase(const ExpList &Sin)
	{
            ASSERTL2(Sin.GetPhysState == true,
                     "Sin physical space is not true ");

            IProductWRTBase(Sin.GetPhys(),m_coeffs);
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
                (*m_exp)[i]->IProductWRTBase(inarray+cnt,elmt_outarray);
                cnt  += (*m_exp)[i]->GetTotPoints();
                cnt1 += (*m_exp)[i]->GetNcoeffs();
	    }
	    m_transState = eLocal;
	}
	

        void ExpList::PhysDeriv(ExpList &out_d0, 
                                ExpList &out_d1, 
                                ExpList &out_d2)
        {
            ASSERTL2(m_physState == true,
                     "local physical space is not true ");
            PhysDeriv(m_phys,
                      out_d0.UpdatePhys(), 
                      out_d0.UpdatePhys(), 
                      out_d0.UpdatePhys());
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
                
		(*m_exp)[i]->PhysDeriv(inarray,e_out_d0,e_out_d1,e_out_d2);
		cnt  += (*m_exp)[i]->GetTotPoints();
	    }
	}

        void ExpList::FwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetPhysState == true,
                     "Sin physical space is not true ");

            FwdTrans(Sin.GetPhys(),m_coeffs);
	    m_transState = eLocal;
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
                (*m_exp)[i]->FwdTrans(inarray+cnt, e_outarray);
		cnt  += (*m_exp)[i]->GetTotPoints();
		cnt1 += (*m_exp)[i]->GetNcoeffs();
	    }
	}



        void ExpList::BwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetTransState() == eLocal ||
                     Sin.GetTransState() == eLocalCont, 
                     "Error input state is not in transformed space");

            BwdTrans(Sin.GetCoeffs(),m_phys);
	    m_physState = true;
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
		(*m_exp)[i]->BwdTrans(inarray + cnt, e_outarray);
		cnt   += (*m_exp)[i]->GetNcoeffs();
		cnt1  += (*m_exp)[i]->GetTotPoints();
	    }	    
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
		    (*m_exp)[i]->GetCoords(e_coord_0);
		    cnt  += (*m_exp)[i]->GetTotPoints();
		}
		break;
	    case 2: 
		ASSERTL0(coord_1 != NullNekDoubleSharedArray, 
			 "output coord_1 is not defined");
    
		for(i= 0; i < GetExpSize(); ++i)
		{
                    e_coord_0 = coord_0 + cnt;
                    e_coord_1 = coord_1 + cnt;
		    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1);
		    cnt  += (*m_exp)[i]->GetTotPoints();
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

		    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1,e_coord_2);
		    cnt  += (*m_exp)[i]->GetTotPoints();
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
	    
	    (*m_exp)[0]->WriteToFile(out,1);
	    
	    for(i= 1; i < GetExpSize(); ++i)
	    {
		(*m_exp)[i]->WriteToFile(out,0);
	    }
	}
	
	NekDouble  ExpList::Linf(const ExpList &Sol)
	{
            ASSERTL2(Sol.GetPhysState == true,
                     "local physical space is not true ");

	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    NekDouble err = 0.0;
	    int       i,cnt = 0;
	    ConstNekDoubleSharedArray soln = Sol.GetPhys();
	    ConstNekDoubleSharedArray phys = m_phys;

	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(phys+cnt);
		err  = std::max(err,(*m_exp)[i]->Linf(soln + cnt));
		cnt  += (*m_exp)[i]->GetTotPoints();
	    }
	    
	    return err;
	}
	
	NekDouble  ExpList::L2(const ExpList &Sol)
	{
            ASSERTL2(Sol.GetPhysState == true,
                     "local physical space is not true ");

	    NekDouble err = 0.0,errl2;
	    int    i,cnt = 0;
	    ConstNekDoubleSharedArray soln = Sol.GetPhys();
	    ConstNekDoubleSharedArray phys = m_phys;
	    
	    if(m_physState == false)
	    {
		BwdTrans(*this);
	    }
	    
	    for(i= 0; i < GetExpSize(); ++i)
	    {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(phys+cnt);
		errl2 = (*m_exp)[i]->L2(soln+cnt);
		err += errl2*errl2;
		cnt  += (*m_exp)[i]->GetTotPoints();
	    }
	    
	    return sqrt(err);
	}
	
    } //end of namespace
} //end of namespace

