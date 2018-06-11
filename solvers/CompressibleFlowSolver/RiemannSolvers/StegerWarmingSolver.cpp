///////////////////////////////////////////////////////////////////////////////
//
// File: StegerWarmingSolver.cpp
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
// Description: StegerWarming Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/StegerWarmingSolver.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
//#include <>

namespace Nektar
{
    std::string StegerWarmingSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "StegerWarming",
			StegerWarmingSolver::create,
            "StegerWarming Riemann solver");

    StegerWarmingSolver::StegerWarmingSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {
    }
    
    /**
     * @brief StegerWarming Riemann solver
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.
     * @param rhouL     x-momentum component left state.
     * @param rhouR     x-momentum component right state.
     * @param rhovL     y-momentum component left state.
     * @param rhovR     y-momentum component right state.
     * @param rhowL     z-momentum component left state.
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.
     * @param ER        Energy right state.
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component
     * @param rhovf     Computed Riemann flux for y-momentum component
     * @param rhowf     Computed Riemann flux for z-momentum component
     * @param Ef        Computed Riemann flux for energy.
     */
    void StegerWarmingSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        NekDouble nx,ny,nz;
        NekDouble f1,f2,f3,f4,f5;
        nx = 1.0;
        ny = 0.0;
        nz = 0.0;

        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;
        fsw =  1.0;
        flux_sw_pn(
            rhoL,   rhouL,  rhovL,  rhowL,  EL,
            nx  ,   ny   ,  nz,
            f1  ,   f2   ,  f3   ,  f4   ,  f5,
            efix_StegerWarming,   fsw);

        rhof    = f1;
        rhouf   = f2;
        rhovf   = f3;
        rhowf   = f4;
        Ef      = f5;

        fsw =-1.0;
        flux_sw_pn(
            rhoR,   rhouR,  rhovR,  rhowR,  ER,
            nx  ,   ny   ,  nz,
            f1  ,   f2   ,  f3   ,  f4   ,  f5,
            efix_StegerWarming,   fsw);

        rhof    += f1;
        rhouf   += f2;
        rhovf   += f3;
        rhowf   += f4;
        Ef      += f5;

        if(false)
        {
            // int nvariables=5;
            // Array<OneD, NekDouble> PointFwd(nvariables,0.0),PointBwd(nvariables,0.0);
            // Array<OneD, NekDouble> PointFlux(nvariables,0.0);
            // Array<OneD, NekDouble> PointNormal(3,0.0);

            // DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
            //     ::AllocateSharedPtr(nvariables, nvariables);
            // DNekMatSharedPtr PointBJac = MemoryManager<DNekMat>
            //     ::AllocateSharedPtr(nvariables, nvariables);

            // PointNormal[0] = nx;
            // PointNormal[1] = ny;
            // PointNormal[2] = nz;
            
            
            // PointFwd[0] = rhoL;
            // PointFwd[1] = rhouL;
            // PointFwd[2] = rhovL;
            // PointFwd[3] = rhowL;
            // PointFwd[4] = EL;

            // PointBwd[0] = rhoR;
            // PointBwd[1] = rhouR;
            // PointBwd[2] = rhovR;
            // PointBwd[3] = rhowR;
            // PointBwd[4] = ER;
            // v_PointFluxJacobian(nvariables,PointFwd,PointBwd,PointNormal,PointFJac,PointBJac);


            // DNekMat &M = (*PointBJac);
            // NekVector<NekDouble> VectPrim(nvariables,PointBwd,eWrapper);
            // NekVector<NekDouble> VectFlux(nvariables,PointFlux,eWrapper);
            // VectFlux = M * VectPrim;
            // // for(int i =0;i<nvariables;i++)
            // // {
            // //     PointFlux[i] = 0.0;
            // //     for(int j =0;j<nvariables;j++)
            // //     {
            // //         PointFlux[i] += M(j,i)*VectPrim[j];
            // //     }
            // // }


            // std::cout   <<std::scientific<<std::setw(12)<<std::setprecision(5)
            //             <<"abs(PointFlux[0]-f1)   =   "<<abs(PointFlux[0]-f1)<<"    "<<PointFlux[0]<<"    "<<f1<<std::endl
            //             <<"abs(PointFlux[1]-f2)   =   "<<abs(PointFlux[1]-f2)<<"    "<<PointFlux[1]<<"    "<<f2<<std::endl
            //             <<"abs(PointFlux[2]-f3)   =   "<<abs(PointFlux[2]-f3)<<"    "<<PointFlux[2]<<"    "<<f3<<std::endl
            //             <<"abs(PointFlux[3]-f4)   =   "<<abs(PointFlux[3]-f4)<<"    "<<PointFlux[3]<<"    "<<f4<<std::endl
            //             <<"abs(PointFlux[4]-f5)   =   "<<abs(PointFlux[4]-f5)<<"    "<<PointFlux[4]<<"    "<<f5<<std::endl;

        }
    }


    void StegerWarmingSolver::flux_sw_pn(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  nx  , double  ny   , double  nz   ,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef,
        NekDouble efix, NekDouble fsw)
    {
        NekDouble   ro = rhoL;
        NekDouble   vx = rhouL / rhoL;
        NekDouble   vy = rhovL / rhoL;
        NekDouble   vz = rhowL / rhoL;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * vx + rhovL * vy + rhowL * vz)) / rhoL;


        NekDouble   ps = m_eos->GetPressure(rhoL, eL);
        NekDouble   c  = m_eos->GetSoundSpeed(rhoL,eL);
        NekDouble   T  = m_eos->GetTemperature(rhoL, eL);
        NekDouble   h  = m_eos->GetEnthalpy(T);
        NekDouble   c2 = c*c;
        NekDouble   v2 = vx*vx + vy*vy + vz*vz;
        NekDouble   h0 = h + 0.5*v2;
        NekDouble   e0 = eL + 0.5*v2;
        
        NekDouble sml_ssf= 1.0E-12;


        NekDouble   vn = nx*vx + ny*vy + nz*vz;
        // NekDouble   sn = std::numeric_limits::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        NekDouble   sn = std::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        NekDouble   osn = 1.0/sn;
        NekDouble   nxa = nx * osn;
        NekDouble   nya = ny * osn;
        NekDouble   nza = nz * osn;
        NekDouble   vna = vn * osn;
        NekDouble   l1 = vn;
        NekDouble   l4 = vn + sn*c;
        NekDouble   l5 = vn - sn*c;

        NekDouble   eps = efix*sn;
        NekDouble   eps2 = eps*eps;
        NekDouble   al1 = sqrt(l1*l1 + eps2);
        NekDouble   al4 = sqrt(l4*l4 + eps2);
        NekDouble   al5 = sqrt(l5*l5 + eps2);

        l1 = 0.5*(l1 + fsw*al1);
        l4 = 0.5*(l4 + fsw*al4);
        l5 = 0.5*(l5 + fsw*al5);

        
        NekDouble   c2r = c2 / m_eos->GetGamma();
        NekDouble   x1 = c2r * ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 );
        NekDouble   x2 = c2r * ( l4 - l5 )/( 2.0 * c );

        rhof    = (l1 - x1 ) * ro;
        rhouf   = (l1*vx - x1*vx + nxa*x2 ) * ro;
        rhovf   = (l1*vy - x1*vy + nya*x2 ) * ro;
        rhowf   = (l1*vz - x1*vz + nza*x2 ) * ro;
        Ef      = (l1*e0 - x1*h0 + vna*x2 ) * ro;
    }
    

    void StegerWarmingSolver::v_PointFluxJacobian(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &Bwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
                  DNekMatSharedPtr       &BJac)
    {
        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;

        fsw = 1.0;
        PointFluxJacobian_pn(Fwd,normals,FJac,efix_StegerWarming,fsw);
  
        fsw = -1.0;
        PointFluxJacobian_pn(Bwd,normals,BJac,efix_StegerWarming,fsw);
        return;
    }

    // Currently duplacate in compressibleFlowSys
    // if fsw=+-1 calculate the steger-Warming flux vector splitting flux Jacobian
    // if fsw=0   calculate the Jacobian of the exact flux 
    // efix is the numerical flux entropy fix parameter
    void StegerWarmingSolver::PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw1)
    {
            NekDouble ro,vx,vy,vz,ps,gama,ae ;
            NekDouble a,a2,h,h0,v2,vn,eps,eps2;
            NekDouble nx,ny,nz;
            NekDouble sn,osn,nxa,nya,nza,vna;
            NekDouble l1,l4,l5,al1,al4,al5,x1,x2,x3,y1;
            NekDouble c1,d1,c2,d2,c3,d3,c4,d4,c5,d5;
            NekDouble sml_ssf= 1.0E-12;

            NekDouble fsw = fsw1;

            NekDouble fExactorSplt = 2.0-abs(fsw); // if fsw=+-1 calculate 

            NekDouble   rhoL  = Fwd[0];
            NekDouble   rhouL = Fwd[1];
            NekDouble   rhovL = Fwd[2];
            NekDouble   rhowL = Fwd[3];
            NekDouble   EL    = Fwd[4];

            ro = rhoL;
            vx = rhouL / rhoL;
            vy = rhovL / rhoL;
            vz = rhowL / rhoL;

            // Internal energy (per unit mass)
            NekDouble eL =
                    (EL - 0.5 * (rhouL * vx + rhovL * vy + rhowL * vz)) / rhoL;
            ps = m_eos->GetPressure(rhoL, eL);
            gama = m_eos->GetGamma();

            ae = gama - 1.0;
            v2 = vx*vx + vy*vy + vz*vz;
            a2 = gama*ps/ro;
            h = a2/ae;

            h0 = h + 0.5*v2;
            a = sqrt(a2);

            nx = normals[0];
            ny = normals[1];
            nz = normals[2];
            vn = nx*vx + ny*vy + nz*vz;

            sn = std::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
            osn = 1.0/sn;

            nxa = nx * osn;
            nya = ny * osn;
            nza = nz * osn;
            vna = vn * osn;

            // if(vn<0)
            // {
            //     a = -a;
            // }

            l1 = vn;
            l4 = vn + sn*a;
            l5 = vn - sn*a;

            eps = efix*sn;
            eps2 = eps*eps;

            al1 = sqrt(l1*l1 + eps2);
            al4 = sqrt(l4*l4 + eps2);
            al5 = sqrt(l5*l5 + eps2);

            l1 = 0.5*(fExactorSplt*l1 + fsw*al1);
            l4 = 0.5*(fExactorSplt*l4 + fsw*al4);
            l5 = 0.5*(fExactorSplt*l5 + fsw*al5);

            x1 = 0.5*(l4 + l5);
            x2 = 0.5*(l4 - l5);
            x3 = x1 - l1;
            y1 = 0.5*v2;
            c1 = ae*x3/a2;
            d1 = x2/a;

            int nsf = 0;
            (*FJac)(nsf  ,nsf  ) = c1*y1 - d1*vna + l1;
            (*FJac)(nsf  ,nsf+1) = -c1*vx + d1*nxa;
            (*FJac)(nsf  ,nsf+2) = -c1*vy + d1*nya;
            (*FJac)(nsf  ,nsf+3) = -c1*vz + d1*nza;
            (*FJac)(nsf  ,nsf+4) = c1;
            c2 = c1*vx + d1*nxa*ae;
            d2 = x3*nxa + d1*vx;
            (*FJac)(nsf+1,nsf  ) = c2*y1 - d2*vna;
            (*FJac)(nsf+1,nsf+1) = -c2*vx + d2*nxa + l1;
            (*FJac)(nsf+1,nsf+2) = -c2*vy + d2*nya;
            (*FJac)(nsf+1,nsf+3) = -c2*vz + d2*nza;
            (*FJac)(nsf+1,nsf+4) = c2;
            c3 = c1*vy + d1*nya*ae;
            d3 = x3*nya + d1*vy;
            (*FJac)(nsf+2,nsf  ) = c3*y1 - d3*vna;
            (*FJac)(nsf+2,nsf+1) = -c3*vx + d3*nxa;
            (*FJac)(nsf+2,nsf+2) = -c3*vy + d3*nya + l1;
            (*FJac)(nsf+2,nsf+3) = -c3*vz + d3*nza;
            (*FJac)(nsf+2,nsf+4) = c3;
            c4 = c1*vz + d1*nza*ae;
            d4 = x3*nza + d1*vz;
            (*FJac)(nsf+3,nsf  ) = c4*y1 - d4*vna;
            (*FJac)(nsf+3,nsf+1) = -c4*vx + d4*nxa;
            (*FJac)(nsf+3,nsf+2) = -c4*vy + d4*nya;
            (*FJac)(nsf+3,nsf+3) = -c4*vz + d4*nza + l1;
            (*FJac)(nsf+3,nsf+4) = c4;
            c5 = c1*h0 + d1*vna*ae;
            d5 = x3*vna + d1*h0;
            (*FJac)(nsf+4,nsf  ) = c5*y1 - d5*vna;
            (*FJac)(nsf+4,nsf+1) = -c5*vx + d5*nxa;
            (*FJac)(nsf+4,nsf+2) = -c5*vy + d5*nya;
            (*FJac)(nsf+4,nsf+3) = -c5*vz + d5*nza;
            (*FJac)(nsf+4,nsf+4) = c5 + l1;

    }
}
