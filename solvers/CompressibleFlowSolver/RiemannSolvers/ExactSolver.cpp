///////////////////////////////////////////////////////////////////////////////
//
// File: ExactSolver.cpp
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
// Description: Exact Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/ExactSolver.h>

namespace Nektar
{
    std::string ExactSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Exact",
            ExactSolver::create,
            "Exact Riemann solver");

    ExactSolver::ExactSolver() : CompressibleSolver()
    {

    }
    
    void ExactSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        // Exact Riemann Solver (GOTTLIEB AND GROTH - 1987; TORO - 1998) 
        NekDouble gamma = m_params["gamma"]();
        
        // gamma operation
        NekDouble f1=(gamma-1.0)/(2.0*gamma);
        NekDouble f2=2.0/(gamma-1.0);
        NekDouble f3=(gamma+1.0)/(2.0*gamma);
        NekDouble f4=(gamma+1.0)/4.0;

        // Right and left variables
        NekDouble uL = rhouL/rhoL;
        NekDouble vL = rhovL/rhoL;
        NekDouble wL = rhowL/rhoL;
        NekDouble uR = rhouR/rhoR;
        NekDouble vR = rhovR/rhoR;
        NekDouble wR = rhowR/rhoR;
        NekDouble pL = (gamma - 1.0) * (EL - 0.5 * (rhouL*uL + rhovL*vL + rhowL*wL));
        NekDouble pR = (gamma - 1.0) * (ER - 0.5 * (rhouR*uR + rhovR*vR + rhowR*wR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);

        // Five possible configuration: rcs=1,scr=2,scs=3,rcr=4,rcvcr=5
        // r = rarefaction
        // c = contact discontinuity
        // s = shock
        // v = vacuum

        NekDouble pratio,z;
        NekDouble u_ncr,u_rcn,u_scn,u_ncs;
        NekDouble u_rcvr = uL +(f2*cL)+(f2*cR);
        int pattern = -1;
    
        if(uR<u_rcvr)
        {
            if(uR>=uL)
            {
                if(pR>=pL)
                { 
                    pratio=pL/pR;
                    u_ncr=uL+(f2*cR)*(1.0-pow(pratio,f1));
                    if(uR>=u_ncr)
                        pattern=4;
                    else
                        pattern=2;
                }
                else if(pR<pL)
                {
                    pratio=pR/pL;
                    u_rcn=uL+(f2*cL)*(1.0-pow(pratio,f1));
                    if(uR>=u_rcn)
                        pattern=4;
                    else
                        pattern=1;
                }
            }
            else if(uR<uL)
            {
                if(pR>=pL)
                {
                    pratio=pR/pL;
                    u_scn=uL-((cL/gamma)*(pratio-1.0)/sqrt((f3*pratio)-f1));
                    if(uR>=u_scn)
                        pattern=2;
                    else
                        pattern=3;
                }
                else if(pR<pL)
                {
                    pratio=pL/pR;
                    u_ncs=uL-((cR/gamma)*(pratio-1.0)/sqrt((f3*pratio)-f1));
                    if(uR>=u_ncs)
                        pattern=1;
                    else
                        pattern=3;
                }
            }
        }
        else if(uR>=u_rcvr)
            pattern=5;
    
        // Initial Guess for u_int
        NekDouble aL = (gamma*pL)/cL;
        NekDouble aR = (gamma*pR)/cR;
        if(pL>=pR)
            z = (f2/f2)*(cR/cL)*pow((pL/pR),f1);
        else if(pL<pR)
            z = (f2/f2)*(cR/cL)*pow((pL/pR),f1);
        NekDouble u_int = ((z*(uL+(f2*cL)))+(uR-f2*cR))/(1.0+z); 

        /*   PATTERN RAREFACTION WAVE - CONTACT SURFACE - SHOCK  */
    
        NekDouble c_intR,c_intL,derp_intR,derp_intL,wR2,wL2,p_int,p_intL,p_intR,u_intL,u_intR;
        NekDouble swave1,swave2,swave3,swave4,swave5; 
        NekDouble u0,u1,u2,u3,u4;
        NekDouble uaux,aaux,paux;
        NekDouble chi = 0.0;
        NekDouble EPSILON = 1.0e-6;
    
        if(pattern==1)
        {
            p_intR = 1.0;
            p_intL = p_intR*(1.0+10.0*EPSILON);
            while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
            { 
                c_intL    = cL-((u_int-uL)/f2);
                p_intL    = pL*pow((c_intL/cL),(1.0/f1));
                derp_intL = (-gamma*p_intL)/c_intL;
                wR2       = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
                p_intR    = pR+(aR*wR2*(u_int-uR));
                derp_intR = (2.0*aR*pow(wR2,3.0))/(1.0+(wR2*wR2));
                u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
            }
            p_int  = (p_intL+p_intR)/2.0;
            c_intR = cR*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pR))/(gamma+1.0+(gamma-1.0)*(pR/p_int)));
            c_intL = cL-((u_int-uL)/f2);
            wR2    = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
            swave1 = (wR2*cR)+uR;
            swave3 = u_int;
            swave4 = u_int-c_intL;
            swave5 = uL-cL;
            if(chi>=swave1)
            {
                u0 = gamma*(pR/(cR*cR));
                u1 = u0*uR;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR + vR*vR + wR*wR)/(cR*cR)));
            }
            else if((chi<swave1)&&(chi>=swave3))
            {
                u0 = gamma*(p_int/(c_intR*c_intR));
                u1 = u0*u_int;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int + vR*vR + wR*wR)/(c_intR*c_intR)));
            }
            else if((chi<swave3)&&(chi>=swave4))
            {
                u0 = gamma*(p_int/(c_intL*c_intL));
                u1 = u0*u_int;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int + vL*vL + wL*wL)/(c_intL*c_intL)));
            }
            else if((chi<swave4)&&(chi>=swave5))
            {
                uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
                aaux = cL+((1.0/f2)*(uL-uaux));
                paux = pL*pow((aaux/cL),(1.0/f1));
                u0   = (gamma*paux)/(aaux*aaux);
                u1   = u0*uaux;
                u2   = u0*vL;
                u3   = u0*wL;
                u4   = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL+wL*wL)*0.5);
            }
            else if(chi<swave5)
            {
                u0 = gamma*(pL/(cL*cL));
                u1 = u0*uL;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL+wL*wL)/(cL*cL)));
            }    
        }
    
        /*        PATTERN:   SHOCK - CONTACT SURFACE - RAREFACTION WAVE    */
    
        else if(pattern==2)
        {
            p_intR = 1.0;
            p_intL = p_intR*(1.0+10.0*EPSILON);
            while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
            {    
                wL2       = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
                p_intL    = pL+(aL*wL2*(u_int-uL));
                derp_intL = (2.0*aL*pow(wL2,3.0))/(1.0+(wL2*wL2));
                c_intR    = cR+((u_int-uR)/f2);
                p_intR    = pR*pow((c_intR/cR),(1.0/f1));
                derp_intR = (gamma*p_intR)/c_intR;
                u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
            }
            p_int  = (p_intL+p_intR)/2.0;
            c_intL = cL*sqrt((gamma+1.0+((gamma-1.0)*(p_int/pL)))/(gamma+1.0+((gamma-1.0)*(pL/p_int))));
            wL2    = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
            swave1 = uR+cR;
            swave2 = u_int+c_intR;
            swave3 = u_int;
            swave5 = (wL2*cL)+uL;
            if(chi>=swave1)
            {
                u0 = gamma*(pR/(cR*cR));
                u1 = u0*uR;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR+wR*wR)/(cR*cR)));
            }
            else if((chi<swave1)&&(chi>=swave2))
            { 
                uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
                aaux = cR+((1.0/f2)*(uaux-uR));
                paux = pR*pow((aaux/cR),(1.0/f1));
                u0 = (gamma*paux)/(aaux*aaux);
                u1 = u0*uaux;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR+wR*wR)*0.5);
            }
            else if((chi<swave2)&&(chi>=swave3))
            { 
                u0 = gamma*(p_int/(c_intR*c_intR));
                u1 = u0*u_int;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (p_int/(gamma-1.0))+(u0*(u_int*u_int+vR*vR+wR*wR)*0.5);
            }
            else if((chi<swave3)&&(chi>=swave5))
            { 
                u0 = gamma*(p_int/(c_intL*c_intL));
                u1 = u0*u_int;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (p_int/(gamma-1.0))+(u0*(u_int*u_int+vL*vL+wL*wL)*0.5);
            }
            else if(chi<swave5)
            { 
                u0 = gamma*(pL/(cL*cL));
                u1 = u0*uL;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL+wL*wL)/(cL*cL)));
            }
        }

        /*       PATTERN: SHOCK - CONTACT SURFACE -SHOCK      */
    
        else if(pattern==3)
        {
            p_intR = 1.0;
            p_intL = p_intR*(1.0+10.0*EPSILON);
            while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
            {      
                wL2       = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
                p_intL    = pL+(aL*wL2*(u_int-uL));
                derp_intL = (2.0*aL*pow(wL2,3.0))/(1.0+(wL2*wL2));
                wR2       = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
                p_intR    = pR+(aR*wR2*(u_int-uR));
                derp_intR = (2.0*aR*pow(wR2,3.0))/(1.0+(wR2*wR2));
                u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
            }
            p_int  = (p_intL+p_intR)/2.0;
            c_intL = cL*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pL))/(gamma+1.0+(gamma-1.0)*(pL/p_int)));
            c_intR = cR*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pR))/(gamma+1.0+(gamma-1.0)*(pR/p_int)));
            wR2    = f4*((u_int-uR)/cR)+sqrt(1+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
            wL2    = f4*((u_int-uL)/cL)-sqrt(1+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
            swave1 = (wR2*cR)+uR;
            swave3 = u_int;
            swave5 = (wL2*cL)+uL;
            if(chi>=swave1)
            {
                u0 = gamma*(pR/(cR*cR));
                u1 = u0*uR;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR+wR*wR)/(cR*cR)));
            }
            else if((chi<swave1)&&(chi>=swave3))
            {
                u0 = gamma*(p_int/(c_intR*c_intR));
                u1 = u0*u_int;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (p_int/(gamma-1.0))+(0.5*u0*(u_int*u_int+vR*vR+wR*wR));
            }
            else if((chi<swave3)&&(chi>=swave5))
            {
                u0 = gamma*(p_int/(c_intL*c_intL));
                u1 = u0*u_int;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (p_int/(gamma-1.0))+(0.5*u0*(u_int*u_int+vL*vL+wL*wL));
            }
            else if(chi<swave5)
            {
                u0 = gamma*(pL/(cL*cL));
                u1 = u0*uL;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL)/(cL*cL+vL*vL+wL*wL)));
            }
        }

        /*    PATTERN:  RAREFACTION WAVE - CONTAC SURFACE - RAREFACTION WAVE   */
    
        else if(pattern==4)
        { 
            p_intR=1.0;
            p_intL=p_intR*(1.0+10.0*EPSILON);
            while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
            {
                c_intL=cL-((u_int-uL)/f2);
                p_intL=pL*pow((c_intL/cL),(1.0/f1));
                derp_intL=(-gamma*p_intL)/c_intL;
                c_intR=cR+((u_int-uR)/f2);
                p_intR=pR*pow((c_intR/cR),(1.0/f1));
                derp_intR=(gamma*p_intR)/c_intR;
                u_int=u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
            }
            p_int=(p_intL+p_intR)/2.0;
            c_intL=cL-((u_int-uL)/f2);
            c_intR=cR+((u_int-uR)/f2);
            swave1=uR+cR;
            swave2=u_int+c_intR;
            swave3=u_int;
            swave4=u_int-c_intL;
            swave5=uL-cL;
	  
            if(chi>=swave1)
	    {
                u0 = gamma*(pR/(cR*cR));
                u1 = u0*uR;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR+wR*wR)/(cR*cR)));
	    }
            else if((chi<swave1)&&(chi>=swave2))
	    {
                uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
                aaux = cR+((1.0/f2)*(uaux-uR));
                paux = pR*pow((aaux/cR),(1.0/f1));
                u0 = (gamma*paux)/(aaux*aaux);
                u1 = u0*uaux;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR+wR*wR)*.5);  
	    }
            else if((chi<swave2)&&(chi>=swave3))
	    {
                u0 = gamma*(p_int/(c_intR*c_intR));
                u1 = u0*u_int;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int+vR*vR+wR*wR)/(c_intR*c_intR)));
	    }
            else if((chi<swave3)&&(chi>=swave4))
	    {
                u0 = gamma*(p_int/(c_intL*c_intL));
                u1 = u0*u_int;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int+vL*vL+wL*wL)/(c_intL*c_intL)));
	    }
            else if((chi<swave4)&&(chi>=swave5))
	    {
                uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
                aaux = cL+((1.0/f2)*(uL-uaux));
                paux = pL*pow((aaux/cL),(1.0/f1));
                u0 = (gamma*paux)/(aaux*aaux);
                u1 = u0*uaux;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL+wL*wL)*.5);  
	    }
            else if(chi<swave5)
	    {
                u0 = gamma*(pL/(cL*cL));
                u1 = u0*uL;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL+wL*wL)/(cL*cL)));
	    }
        }
    
        /*PATTERN:   RAREFACTION WAVE - CONTACT SURFACE - VACUUM -
          CONTACT SURFACE - RAREFACTION WAVE    */
        else if(pattern==5)
        {
            p_int  = 0.0;
            u_intR = uR-(f2*cR);
            u_intL = uL+(f2*cL);
            c_intR = 0.0;
            c_intL = 0.0;
            swave1 = uR+cR;
            swave2 = u_intR+c_intR;
            swave4 = u_intL-c_intL;
            swave5 = uL-cL;
            if(chi>=swave1)
            {
                u0 = gamma*(pR/(cR*cR));
                u1 = u0*uR;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR+wR*wR)/(cR*cR)));
            }
            else if((chi<swave1)&&(chi>=swave2))
            {
                uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
                aaux = cR+((1.0/f2)*(uaux-uR));
                paux = pR*pow((aaux/cR),(1.0/f1));
                u0 = (gamma*paux)/(aaux*aaux);
                u1 = u0*uaux;
                u2 = u0*vR;
                u3 = u0*wR;
                u4 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR+wR*wR)*.5);
            }
            else if((chi<swave2)&&(chi>=swave4))
            {
                u0 = 0.0;
                u1 = 0.0;
                u2 = 0.0;
                u3 = 0.0;
                u4 = 0.0;
            }
            else if((chi<swave4)&&(chi>=swave5))
            {
                uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
                aaux = cL+((1.0/f2)*(uL-uaux));
                paux = pL*pow((aaux/cL),(1.0/f1));
                u0 = (gamma*paux)/(aaux*aaux);
                u1 = u0*uaux;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL+wL*wL)*.5);
            }
            else if(chi<swave5)
            {
                u0 = gamma*(pL/(cL*cL));
                u1 = u0*uL;
                u2 = u0*vL;
                u3 = u0*wL;
                u4 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL+wL*wL)/(cL*cL)));
            }
        }
    
        rhof  = u1;
        rhouf = u1*u1/u0 + (gamma-1.0)*(u4 - 0.5*(u1*u1/u0+u2*u2/u0+u3*u3/u0));
        rhovf = u1*u2/u0;
        rhowf = u1*u3/u0;
        Ef    = u1/u0 * (gamma*u4 - (gamma-1.0)*0.5*(u1*u1/u0+u2*u2/u0+u3*u3/u0));
    }
}
