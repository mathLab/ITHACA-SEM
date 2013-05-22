///////////////////////////////////////////////////////////////////////////////
//
// File FentonKarma2b.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) n6 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Courtemanche-Ramirez-Nattel ionic atrial cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <CardiacEPSolver/CellModels/FentonKarma2b.h>

using namespace std;

namespace Nektar
{
    std::string FentonKarma2b::className
    = GetCellModelFactory().RegisterCreatorFunction(
                                                    "FentonKarma2b",
                                                    FentonKarma2b::create,
                                                    "Phenomenological Model.");
    
    
    /**
     *
     */
    FentonKarma2b::FentonKarma2b(
                                 const LibUtilities::SessionReaderSharedPtr& pSession,
                                 const MultiRegions::ExpListSharedPtr& pField)
    : CellModel(pSession, pField)
    {
        C_m = 1; // picoF
        tauvplus=10;
        tauv1minus=100;
        tauy2minus=20;
        tauwplus=800;
        tauwminus=45;
        taud=0.15;
        tau0=1.5;
        taur=31;
        tausi=53;
        k1=10;
        k2=1;
        vcsi=0.7;
        vc=0.25;
        vr=0.6;
        vv=0.05;
        vfi=0.11;
        
        
        
        m_nvar = 3;
        
        // List gates and concentrations
        m_gates.push_back(1); // u
        m_gates.push_back(2); // w

        
        
        
        
    }
    
    
    
    /**
     *
     */
    FentonKarma2b::~FentonKarma2b()
    {
        
    }
    
    
    
    void FentonKarma2b::v_Update(
                                 const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                 Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                 const NekDouble time )
    {
        ASSERTL0(inarray.get() != outarray.get(),
                 "Must have different arrays for input and output.");
        
        // Variables
        //  0   V    membrane potential
        //  1   -    unused
        //  2   u    u gate
        //  3   w    w gate

        
        int n = m_nq;
        int i = 0;
        NekDouble alpha, beta, tauvminus;
        Vmath::Zero(n, outarray[0], 1);
	
	// create temporary Nektar++ arrays
        Array <OneD,NekDouble> myarray1(m_nq,0.0); //Creates an array of size nq and fills it with 0.0
        Array <OneD,NekDouble> myarray2(m_nq,0.0);
        Array <OneD,NekDouble> myarray3(m_nq,0.0);
        Array <OneD,NekDouble> myarray4(m_nq,0.0);
        Array <OneD,NekDouble> myarray5(m_nq,0.0);
        Array <OneD,NekDouble> myarray6(m_nq,0.0);
	Array <OneD,NekDouble> myarray7(m_nq,0.0);
        Array <OneD,NekDouble> myarray8(m_nq,0.0);
        Array <OneD,NekDouble> myarray9(m_nq,0.0);
        Array <OneD,NekDouble> myarray10(m_nq,0.0);
	Array <OneD,NekDouble> myarray11(m_nq,0.0);
        Array <OneD,NekDouble> myarray12(m_nq,0.0);
        Array <OneD,NekDouble> myarray13(m_nq,0.0);


           
        // for p. Have tried several ways of implementing heaviside function; here using logistic function as analytic approx
	Array<OneD, NekDouble> &tmp_xp = myarray1;
        Array<OneD, NekDouble> &tmp_xq = myarray2;
        Array<OneD, NekDouble> &tmp_xr = myarray3;
	Array<OneD, NekDouble> &tmp_xp1 = myarray4;
        Array<OneD, NekDouble> &tmp_xq1 = myarray5;
        Array<OneD, NekDouble> &tmp_xr1 = myarray6;

	Vmath::Sadd(n, -vc,inarray[0],1, tmp_xp, 1);
	Vmath::Smul(n, 1000.0,  tmp_xp, 1, tmp_xp, 1);
	Vmath::Vexp(n, tmp_xp, 1, tmp_xp,1);
        Vmath::Sadd(n, 1.0,tmp_xp,1, tmp_xp1, 1);
	Vmath::Vdiv(n, tmp_xp, 1, tmp_xp1, 1, tmp_xp,1);
        	

	// for q
	Vmath::Sadd(n, -vr,inarray[0],1, tmp_xq, 1);
	Vmath::Smul(n, 1000.0,  tmp_xq, 1, tmp_xq, 1);
	Vmath::Vexp(n, tmp_xq, 1, tmp_xq,1);
        Vmath::Sadd(n, 1.0,tmp_xq,1, tmp_xq1, 1);
        Vmath::Vdiv(n, tmp_xq, 1, tmp_xq1, 1, tmp_xq,1);

	// for r
	Vmath::Sadd(n, -vv,inarray[0],1, tmp_xr, 1);
	Vmath::Smul(n, 1000.0,  tmp_xr, 1, tmp_xr, 1);
	Vmath::Vexp(n, tmp_xr, 1, tmp_xr,1);
        Vmath::Sadd(n, 1.0,tmp_xr,1, tmp_xr1, 1);
        Vmath::Vdiv(n, tmp_xr, 1, tmp_xr1, 1, tmp_xr,1);
      
        //Ifi
        Array<OneD, NekDouble> &tmp_I_fi =  myarray7;
        Array<OneD, NekDouble> &tmp_I_fi2 = myarray8;
        Array<OneD, NekDouble> &tmp_I_fi3 = myarray9;
      
        
        Vmath::Smul(n, -1.0,  inarray[0], 1, tmp_I_fi, 1);
        Vmath::Sadd(n, 1.0,tmp_I_fi,1, tmp_I_fi, 1);
        Vmath::Sadd(n, -vfi,inarray[0],  1, tmp_I_fi2, 1);
        Vmath::Vmul(n,inarray[1], 1, tmp_xp, 1, tmp_I_fi3,1);
        Vmath::Smul(n, -1.0, tmp_I_fi3, 1, tmp_I_fi3,1);
        Vmath::Vmul(n,tmp_I_fi3, 1, tmp_I_fi, 1, tmp_I_fi3,1);
        Vmath::Vmul(n,tmp_I_fi3, 1, tmp_I_fi2, 1, tmp_I_fi3,1);
        Vmath::Sdiv(n, taud, tmp_I_fi3, 1, tmp_I_fi3,1 );
        Vmath::Vsub(n, outarray[0], 1, tmp_I_fi3, 1, outarray[0], 1);
        
        //Iso
        Array<OneD, NekDouble> &tmp_I_so = myarray10;
        Array<OneD, NekDouble> &tmp_I_so1 = myarray11;
        
        Vmath::Smul(n, -1.0,  tmp_xr, 1, tmp_I_so, 1);
        Vmath::Sadd(n, 1.0, tmp_I_so,1, tmp_I_so,1);
        Vmath::Smul(n, -k2, inarray[1],1,tmp_I_so1,1);
        Vmath::Sadd(n, 1.0, tmp_I_so1, 1, tmp_I_so1, 1);
        Vmath::Vmul(n, tmp_I_so, 1, tmp_I_so1, 1, tmp_I_so1, 1);
        Vmath::Vmul(n, tmp_I_so1, 1, inarray[0], 1, tmp_I_so1, 1);
        Vmath::Sdiv(n, tau0, tmp_I_so1, 1, tmp_I_so1, 1);
        Vmath::Sdiv(n, taur, tmp_xr, 1, tmp_I_so, 1);
        Vmath::Vadd(n, tmp_I_so1,1, tmp_I_so, 1, tmp_I_so1, 1); 
        Vmath::Vsub(n, outarray[0], 1, tmp_I_so1, 1, outarray[0], 1);
        
        
        
        //Isi
        Array<OneD, NekDouble> &tmp_I_si = myarray12;
        Array<OneD, NekDouble> &tmp_I_si1 =myarray13;
        
        Vmath::Sadd(n, -vcsi, inarray[0], 1, tmp_I_si, 1);
        Vmath::Smul(n, k1, tmp_I_si,1,tmp_I_si,1);
        Vmath::Smul(n, 2.0, tmp_I_si,1,tmp_I_si,1);
        Vmath::Vexp(n, tmp_I_si, 1, tmp_I_si,1);
        Vmath::Sadd(n, -1.0,tmp_I_si, 1,tmp_I_si1,1);
        Vmath::Sadd(n, 1.0, tmp_I_si, 1,tmp_I_si,1);
        Vmath::Vdiv(n, tmp_I_si1, 1, tmp_I_si,1, tmp_I_si,1);
        Vmath::Sadd(n, 1.0, tmp_I_si,  1, tmp_I_si,1);
        Vmath::Vmul(n, inarray[2], 1, tmp_I_si,1, tmp_I_si,1);
        Vmath::Smul(n, -1.0, tmp_I_si,1,tmp_I_si,1);
        Vmath::Sdiv(n, tausi,tmp_I_si,1, tmp_I_si,1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_si, 1, outarray[0], 1);
        
	// Process gating variables
        const NekDouble * v;
        const NekDouble * x;
        NekDouble * x_tau;
        NekDouble * x_new;
        
         
        // FK u, inarray1
        for (i = 0,  v = &inarray[0][0], x = &inarray[1][0], x_new = &outarray[1][0], x_tau = &m_gates_tau[0][0];
             
             i < n; ++i, ++v, ++x_new, ++x_tau)
        {
	    tauvminus= (*v >=vr) ? tauy2minus : tauv1minus ;
            alpha = (*v >=vc) ? 0.0 : (1/tauvminus);
            beta= (*v >= vc) ? (1/tauvplus) : 0.0; 
            *x_tau = 1.0/(alpha + beta);
            *x_new = alpha*(*x_tau);
        }
        
        
        //FK w, inarray2 
        for (i = 0,  v = &inarray[0][0], x = &inarray[2][0], x_new = &outarray[2][0], x_tau = &m_gates_tau[1][0] ;
             i < n; ++i, ++v, ++x_new, ++x_tau)
        {
            alpha = (*v >=vc) ? 0.0 : (1/tauwminus);
            beta= (*v >= vc) ? (1/tauwplus) : 0.0; 
            *x_tau = 1.0/(alpha + beta);
            *x_new = alpha*(*x_tau);
        }
         
        
        
    }
    
    void FentonKarma2b::v_PrintSummary(std::ostream &out)
    {
        out << "	Cell model      : FentonKarma2b" << std::endl;
    }
    
    
    void FentonKarma2b::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, 0.0,  m_cellSol[0],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[1],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[2],  1);
      
        
    }
}