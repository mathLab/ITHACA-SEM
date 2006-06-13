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

	ExpList::ExpList()
	{
	    m_transState = eNotSet; 
	    m_physState  = false;
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
	double ExpList::Integral(const double *inarray)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    
	    int    cnt = 0;
	    double sum = 0.0;
	    
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    sum += (*def)->Integral(inarray+cnt);
		    cnt += (*def)->GetPointsTot();
		}
	    }
	    return sum; 
	}
	
	
	void ExpList::IProductWRTBase(const double *inarray, double *outarray)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    int    cnt  = 0;
	    int    cnt1 = 0;
	  
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    (*def)->IProductWRTBase(inarray+cnt,outarray+cnt1);
		    cnt  += (*def)->GetPointsTot();
		    cnt1 += (*def)->GetNcoeffs();
		}
	    }
	    m_transState = eLocal;
	}
	
	void ExpList::IProductWRTBase(ExpList &S1, ExpList &S2)
	{
	    IProductWRTBase(S1.GetPhys(),S2.GetCoeffs());
	}
	
	void ExpList::IProductWRTBase(ExpList &S1, double * outarray)
	{
	    IProductWRTBase( S1.GetPhys(),outarray);
	}
	
	void ExpList::Deriv(const int n, double **outarray)
	{
	    Deriv(n,m_phys,outarray);
	}
      
	void ExpList::Deriv(const int n, const double *inarray,
			      double **outarray)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    int    cnt = 0;
	    
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    (*def)->Deriv(n,inarray+cnt,outarray+cnt);
		    cnt  += (*def)->GetPointsTot();
		}
	    }
	}	    

	void ExpList::FwdTrans(const double *inarray)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    int    cnt = 0;
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    (*def)->FwdTrans(inarray+cnt);
		    cnt  += (*def)->GetPointsTot();
		}
		}
	    
	    m_transState = eLocal;
	}
	
	
	void ExpList::BwdTrans(double *outarray)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    int    cnt = 0;
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    (*def)->BwdTrans(outarray+cnt);
		    cnt  += (*def)->GetPointsTot();
		}
	    }
	    m_physState = true;
	}
	
	void ExpList::GetCoords(double **coords)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    int    i, cnt = 0;
	    double *E_coords[3];
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    for(i = 0 ; i < (*def)->GetCoordim(); ++i)
		    {
			    E_coords[i] = coords[i]+cnt;
		    }
		    
		    (*def)->GetCoords(E_coords);
		    cnt  += (*def)->GetPointsTot();
		}
	    }
	}
	
	void ExpList::WriteToFile(std::ofstream &out)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    m_exp_shapes[0][0]->WriteToFile(out,1);
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    (*def)->WriteToFile(out,0);
		}
	    }
	}
	
	double  ExpList::Linf(const double *sol)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    double err = 0.0;
	    int    cnt = 0;
	    
	    if(m_physState == false)
	    {
		v_BwdTrans(m_phys);
	    }
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    err  = std::max(err,(*def)->Linf(sol+cnt));
		    cnt  += (*def)->GetPointsTot();
		}
	    }
	    
	    return err;
	}
	
	
	double  ExpList::L2(const double *sol)
	{
	    std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
	    StdRegions::StdExpansionVectorIter def;
	    double err = 0.0,errl2;
	    int    cnt = 0;
	    
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    for(sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end(); ++sdef)
	    {
		for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		{
		    errl2 = (*def)->L2(sol+cnt);
		    err += errl2*errl2;
		    cnt  += (*def)->GetPointsTot();
		}
	    }
	    
	    return sqrt(err);
	}
	
    } //end of namespace
} //end of namespace

