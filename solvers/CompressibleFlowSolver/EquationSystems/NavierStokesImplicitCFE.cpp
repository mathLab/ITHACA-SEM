///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesImplicitCFE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesImplicitCFE.h>

using namespace std;

namespace Nektar
{
    string NavierStokesImplicitCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesImplicitCFE", NavierStokesImplicitCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesImplicitCFE::NavierStokesImplicitCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          CompressibleFlowSystem(pSession, pGraph),
          NavierStokesCFE(pSession, pGraph),
          CFSImplicit(pSession, pGraph)
    {
    }

    NavierStokesImplicitCFE::~NavierStokesImplicitCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesImplicitCFE::v_InitObject()
    {
        CFSImplicit::v_InitObject(); 

        NavierStokesCFE::InitObject_Explicit();
        
        m_GetdFlux_dDeriv_Array = Array<OneD, GetdFlux_dDeriv> (m_spacedim);
        switch (m_spacedim)
        {
        case 2:
            /* code */
            m_GetdFlux_dDeriv_Array[0] = std::bind(
            &NavierStokesImplicitCFE::GetdFlux_dQx_2D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            m_GetdFlux_dDeriv_Array[1] = std::bind(
            &NavierStokesImplicitCFE::GetdFlux_dQy_2D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);
            break;
        case 3:
            /* code */
            m_GetdFlux_dDeriv_Array[0] = std::bind(
            &NavierStokesImplicitCFE::GetdFlux_dQx_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            m_GetdFlux_dDeriv_Array[1] = std::bind(
            &NavierStokesImplicitCFE::GetdFlux_dQy_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);
            m_GetdFlux_dDeriv_Array[2] = std::bind(
            &NavierStokesImplicitCFE::GetdFlux_dQz_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            break;

        default:

            break;
        }
    }

    void NavierStokesImplicitCFE::v_DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
        size_t nvariables = inarray.size();
        size_t npoints    = GetNpoints();
        size_t ncoeffs    = GetNcoeffs();
        size_t nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff{nvariables};
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>{ncoeffs, 0.0};
        }

        // get artificial viscosity
        if (m_shockCaptureType == "Physical" && m_CalcPhysicalAV)
        {
            GetPhysicalAV(inarray);
        }

        if (m_is_diffIP)
        {
            if (m_bndEvaluateTime < 0.0)
            {
                NEKERROR(ErrorUtil::efatal, "m_bndEvaluateTime not setup");
            }
            m_diffusion->DiffuseCoeffs(nvariables, m_fields, inarray,
                                        outarrayDiff, m_bndEvaluateTime,
                                        pFwd, pBwd);
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(ncoeffs,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff{nvariables - 1};
            Array<OneD, Array<OneD, NekDouble> > inFwd{nvariables - 1};
            Array<OneD, Array<OneD, NekDouble> > inBwd{nvariables-1};

            for (int i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
                inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
                inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables-2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayOfArray ||
                pBwd == NullNekDoubleArrayOfArray)
            {
                inFwd = NullNekDoubleArrayOfArray;
                inBwd = NullNekDoubleArrayOfArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd,    inFwd[0]);
                m_varConv->GetPressure(pBwd,    inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables-2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables-2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff,
                                inFwd, inBwd);

            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }

            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusionCoeff(
                    inarray, outarray);
            }
        }
    }


    /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesImplicitCFE::GetdFlux_dQx_2D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble E=U[3]*orho;
        NekDouble q2=u*u+v*v;
        NekDouble gamma=m_gamma;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr=  1.0/Pr;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE 
        // NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        tmpArray[0+0*nrow]=tmp*(-FourThird*u*nx-v*ny);
        tmpArray[0+1*nrow]=tmp*(FourThird*nx);
        tmpArray[0+2*nrow]=tmp*ny;
        tmpArray[0+3*nrow]=0.0;
        tmpArray[1+0*nrow]=tmp*(-v*nx+TwoThird*u*ny);
        tmpArray[1+1*nrow]=tmp*(-TwoThird*ny);
        tmpArray[1+2*nrow]=tmp*nx;
        tmpArray[1+3*nrow]=0.0;
        tmpArray[2+0*nrow]=(FourThird*u*u+v*v+tmp2*(E-q2))*nx+OneThird*u*v*ny;
        tmpArray[2+0*nrow]=-tmp*(*OutputMatrix)(2,0);
        tmpArray[2+1*nrow]=(FourThird-tmp2)*u*nx-TwoThird*v*ny;
        tmpArray[2+1*nrow]=tmp*(*OutputMatrix)(2,1);
        tmpArray[2+2*nrow]=(1-tmp2)*v*nx+u*ny;
        tmpArray[2+2*nrow]=tmp*(*OutputMatrix)(2,2);
        tmpArray[2+3*nrow]=tmp*tmp2*nx;
    }

     /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dy,drhou_dy,drhov_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesImplicitCFE::GetdFlux_dQy_2D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble E=U[3]*orho;
        NekDouble q2=u*u+v*v;
        NekDouble gamma=m_gamma;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr=  1.0/Pr;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE 
        // NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
           
        tmpArray[0+0*nrow]=tmp*(TwoThird*v*nx-u*ny);
        tmpArray[0+1*nrow]=tmp*ny;
        tmpArray[0+2*nrow]=tmp*(-TwoThird)*nx;
        tmpArray[0+3*nrow]=0.0;
        tmpArray[1+0*nrow]=tmp*(-u*nx-FourThird*v*ny);
        tmpArray[1+1*nrow]=tmp*nx;
        tmpArray[1+2*nrow]=tmp*(FourThird*ny);
        tmpArray[1+3*nrow]=0.0;
        tmpArray[2+0*nrow]=OneThird*u*v*nx+(FourThird*v*v+u*u+tmp2*(E-q2))*ny;
        tmpArray[2+0*nrow]=-tmp*(*OutputMatrix)(2,0);
        tmpArray[2+1*nrow]=(1-tmp2)*u*ny+v*nx;
        tmpArray[2+1*nrow]=tmp*(*OutputMatrix)(2,1);
        tmpArray[2+2*nrow]=(FourThird-tmp2)*v*ny-TwoThird*u*nx;
        tmpArray[2+2*nrow]=tmp*(*OutputMatrix)(2,2);
        tmpArray[2+3*nrow]=tmp*tmp2*ny;
    }

    /**
     * @brief return part of viscous Jacobian derived with 
     * Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=0)= dF_dQx;
     */
    void NavierStokesImplicitCFE::GetdFlux_dQx_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble gamma=m_gamma;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE 
        // NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        tmpArray[0+0*nrow]=tmpx*(-FourThird*u)+tmpy*(-v)+tmpz*(-w);
        tmpArray[0+1*nrow]=tmpx*FourThird;
        tmpArray[0+2*nrow]=tmpy;
        tmpArray[0+3*nrow]=tmpz;
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpx*(-v)+tmpy*(TwoThird*u);
        tmpArray[1+1*nrow]=tmpy*(-TwoThird);
        tmpArray[1+2*nrow]=tmpx;
        tmpArray[1+3*nrow]=0.0;
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpx*(-w)+tmpz*(TwoThird*u);
        tmpArray[2+1*nrow]=tmpz*(-TwoThird);
        tmpArray[2+2*nrow]=0.0;
        tmpArray[2+3*nrow]=tmpx;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=-tmpx*(FourThird*u*u+v*v+w*w+tmp2*(E-q2)) 
                            +tmpy*(-OneThird*u*v)+tmpz*(-OneThird*u*w);
        tmpArray[3+1*nrow]=tmpx*(FourThird-tmp2)*u+tmpy*(-TwoThird*v)
                            +tmpz*(-TwoThird*w);
        tmpArray[3+2*nrow]=tmpx*(1.0-tmp2)*v+tmpy*u;
        tmpArray[3+3*nrow]=tmpx*(1.0-tmp2)*w+tmpz*u;
        tmpArray[3+4*nrow]=tmpx*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian derived with 
     * Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=1)= dF_dQy;
     */
    void NavierStokesImplicitCFE::GetdFlux_dQy_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble gamma=m_gamma;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE 
        // NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        tmpArray[0+0*nrow]=tmpx*(TwoThird*v)+tmpy*(-u);
        tmpArray[0+1*nrow]=tmpy;
        tmpArray[0+2*nrow]=tmpx*(-TwoThird);
        tmpArray[0+3*nrow]=0.0;
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpx*(-u)+tmpy*(-FourThird*v)+tmpz*(-w);
        tmpArray[1+1*nrow]=tmpx;
        tmpArray[1+2*nrow]=tmpy*FourThird;
        tmpArray[1+3*nrow]=tmpz;
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpy*(-w)+tmpz*(TwoThird*v);
        tmpArray[2+1*nrow]=0.0;
        tmpArray[2+2*nrow]=tmpz*(-TwoThird);
        tmpArray[2+3*nrow]=tmpy;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=tmpx*(-OneThird*u*v)-tmpy*(u*u+FourThird*v*v+w*w 
                            +tmp2*(E-q2))+tmpz*(-OneThird*v*w);
        tmpArray[3+1*nrow]=tmpx*v+tmpy*(1-tmp2)*u;
        tmpArray[3+2*nrow]=tmpx*(-TwoThird*u)+tmpy*(FourThird-tmp2)*v
                            +tmpz*(-TwoThird*w);
        tmpArray[3+3*nrow]=tmpy*(1-tmp2)*w+tmpz*v;
        tmpArray[3+4*nrow]=tmpy*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian derived with 
     * Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=2)= dF_dQz;
     */
    void NavierStokesImplicitCFE::GetdFlux_dQz_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble gamma=m_gamma;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE 
        // NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        tmpArray[0+0*nrow]=tmpx*(TwoThird*w)+tmpz*(-u);
        tmpArray[0+1*nrow]=tmpz;
        tmpArray[0+2*nrow]=0.0;
        tmpArray[0+3*nrow]=tmpx*(-TwoThird);
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpy*(TwoThird*w)+tmpz*(-v);
        tmpArray[1+1*nrow]=0.0;
        tmpArray[1+2*nrow]=tmpz;
        tmpArray[1+3*nrow]=tmpy*(-TwoThird);
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpx*(-u)+tmpy*(-v)+tmpz*(-FourThird*w);
        tmpArray[2+1*nrow]=tmpx;
        tmpArray[2+2*nrow]=tmpy;
        tmpArray[2+3*nrow]=tmpz*FourThird;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=tmpx*(-OneThird*u*w)+tmpy*(-OneThird*v*w)
                            -tmpz*(u*u+v*v+FourThird*w*w+tmp2*(E-q2));
        tmpArray[3+1*nrow]=tmpx*w+tmpz*(1-tmp2)*u;
        tmpArray[3+2*nrow]=tmpy*w+tmpz*(1-tmp2)*v;
        tmpArray[3+3*nrow]=tmpx*(-TwoThird*u)+tmpy*(-TwoThird*v)
                            +tmpz*(FourThird-tmp2)*w;
        tmpArray[3+4*nrow]=tmpz*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhoE]
     * Output: 3*4 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
     */
    void NavierStokesImplicitCFE::GetdFlux_dU_2D(
        const Array<OneD, NekDouble>                    &normals,
        const NekDouble                                 mu,
        const NekDouble                                 dmu_dT,
        const Array<OneD, NekDouble>                    &U,
        const Array<OneD, const Array<OneD, NekDouble>> &qfield,
              DNekMatSharedPtr                          &OutputMatrix)
    {
        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble U1=U[0];
        NekDouble U2=U[1];
        NekDouble U3=U[2];
        NekDouble U4=U[3];
        NekDouble dU1_dx=qfield[0][0];
        NekDouble dU2_dx=qfield[0][1];
        NekDouble dU3_dx=qfield[0][2];
        NekDouble dU4_dx=qfield[0][3];
        NekDouble dU1_dy=qfield[1][0];
        NekDouble dU2_dy=qfield[1][1];
        NekDouble dU3_dy=qfield[1][2];
        NekDouble dU4_dy=qfield[1][3];
        NekDouble gamma=m_gamma;
        NekDouble Cv=m_Cv;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;

        NekDouble orho1,orho2,orho3,orho4;
        NekDouble oCv=1./Cv;
        orho1=1.0/U1;
        orho2=orho1*orho1;
        orho3=orho1*orho2;
        orho4=orho2*orho2;

        //Assume Fn=mu*Sn
        //Sn=Sx*nx+Sy*ny

        NekDouble TwoThrid=2./3.;
        NekDouble FourThird=2.0*TwoThrid;
        NekDouble u=U2*orho1;
        NekDouble v=U3*orho1;
        NekDouble du_dx=orho1*(dU2_dx-u*dU1_dx);
        NekDouble dv_dx=orho1*(dU3_dx-v*dU1_dx);
        NekDouble du_dy=orho1*(dU2_dy-u*dU1_dy);
        NekDouble dv_dy=orho1*(dU3_dy-v*dU1_dy);
        NekDouble s12=FourThird*du_dx-TwoThrid*dv_dy;
        NekDouble s13=du_dy+dv_dx;
        NekDouble s22=s13;
        NekDouble s23=FourThird*dv_dy-TwoThrid*du_dx;
        NekDouble snx=s12*nx+s22*ny;
        NekDouble sny=s13*nx+s23*ny;
        NekDouble snv=snx*u+sny*v;
        NekDouble qx=-gamma*mu*oPr*(orho1*dU4_dx-U[3]*orho2*dU1_dx
                        -u*(orho1*dU2_dx-U[1]*orho2*dU1_dx)
                        -v*(orho1*dU3_dx-U[2]*orho2*dU1_dx));
        NekDouble qy=-gamma*mu*oPr*(orho1*dU4_dy-U[3]*orho2*dU1_dy
                        -u*(orho1*dU2_dy-U[1]*orho2*dU1_dy)
                        -v*(orho1*dU3_dy-U[2]*orho2*dU1_dy));
        NekDouble qn=qx*nx+qy*ny;

        //Term1 mu's derivative with U: dmu_dU*Sn
        Array<OneD,NekDouble> tmp(3,0.0);
        tmp[0]=snx;
        tmp[1]=sny;
        tmp[2]=snv-qn/mu;
        Array<OneD,NekDouble> dT_dU (4,0.0);
        dT_dU[0]=oCv*(-orho2*U4+orho3*U2*U2+orho3*U3*U3);
        dT_dU[1]=-oCv*orho2*U2;
        dT_dU[2]=-oCv*orho2*U3;
        dT_dU[3]=oCv*orho1;
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<4;j++)
            {
                tmpArray[i+j*nrow]=dmu_dT*dT_dU[j]*tmp[i];
            }
        }

        //Term 2 +mu*dSn_dU
        NekDouble du_dx_dU1,du_dx_dU2;
        NekDouble du_dy_dU1,du_dy_dU2;
        NekDouble dv_dx_dU1,dv_dx_dU3;
        NekDouble dv_dy_dU1,dv_dy_dU3;
        NekDouble ds12_dU1,ds12_dU2,ds12_dU3;
        NekDouble ds13_dU1,ds13_dU2,ds13_dU3;
        NekDouble ds22_dU1,ds22_dU2,ds22_dU3;
        NekDouble ds23_dU1,ds23_dU2,ds23_dU3;
        NekDouble dsnx_dU1,dsnx_dU2,dsnx_dU3;
        NekDouble dsny_dU1,dsny_dU2,dsny_dU3;
        NekDouble dsnv_dU1,dsnv_dU2,dsnv_dU3;

        du_dx_dU1=-orho2*dU2_dx+2*orho3*U2*dU1_dx;
        du_dx_dU2=-orho2*dU1_dx;
        du_dy_dU1=-orho2*dU2_dy+2*orho3*U2*dU1_dy;
        du_dy_dU2=-orho2*dU1_dy;
        dv_dx_dU1=-orho2*dU3_dx+2*orho3*U3*dU1_dx;
        dv_dx_dU3=du_dx_dU2;
        dv_dy_dU1=-orho2*dU3_dy+2*orho3*U3*dU1_dy;
        dv_dy_dU3=du_dy_dU2;
        ds12_dU1=FourThird*du_dx_dU1-TwoThrid*dv_dy_dU1;
        ds12_dU2=FourThird*du_dx_dU2;
        ds12_dU3=-TwoThrid*dv_dy_dU3;
        ds13_dU1=du_dy_dU1+dv_dx_dU1;
        ds13_dU2=du_dy_dU2;
        ds13_dU3=dv_dx_dU3;
        ds22_dU1=ds13_dU1;
        ds22_dU2=ds13_dU2;
        ds22_dU3=ds13_dU3;
        ds23_dU1=FourThird*dv_dy_dU1-TwoThrid*du_dx_dU1;
        ds23_dU2=-TwoThrid*du_dx_dU2;
        ds23_dU3=FourThird*dv_dy_dU3;
        dsnx_dU1=ds12_dU1*nx+ds22_dU1*ny;
        dsnx_dU2=ds12_dU2*nx+ds22_dU2*ny;
        dsnx_dU3=ds12_dU3*nx+ds22_dU3*ny;
        dsny_dU1=ds13_dU1*nx+ds23_dU1*ny;
        dsny_dU2=ds13_dU2*nx+ds23_dU2*ny;
        dsny_dU3=ds13_dU3*nx+ds23_dU3*ny;
        dsnv_dU1=u*dsnx_dU1+v*dsny_dU1-orho2*U2*snx-orho2*U3*sny;
        dsnv_dU2=u*dsnx_dU2+v*dsny_dU2+orho1*snx;
        dsnv_dU3=u*dsnx_dU3+v*dsny_dU3+orho1*sny;
        tmpArray[0+0*nrow]=tmpArray[0+0*nrow]+mu*dsnx_dU1;
        tmpArray[0+1*nrow]=tmpArray[0+1*nrow]+mu*dsnx_dU2;
        tmpArray[0+2*nrow]=tmpArray[0+2*nrow]+mu*dsnx_dU3;
        tmpArray[1+0*nrow]=tmpArray[1+0*nrow]+mu*dsny_dU1;
        tmpArray[1+1*nrow]=tmpArray[1+1*nrow]+mu*dsny_dU2;
        tmpArray[1+2*nrow]=tmpArray[1+2*nrow]+mu*dsny_dU3;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]+mu*dsnv_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]+mu*dsnv_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]+mu*dsnv_dU3;

        //Consider +qn's effect (does not include mu's effect)
        NekDouble dqx_dU1,dqx_dU2,dqx_dU3,dqx_dU4;
        NekDouble dqy_dU1,dqy_dU2,dqy_dU3,dqy_dU4;
        NekDouble tmpx=-nx*mu*gamma*oPr;
        dqx_dU1=tmpx*(-orho2*dU4_dx+2*orho3*U4*dU1_dx+2*orho3*U2*dU2_dx
                -3*orho4*U2*U2*dU1_dx+2*orho3*U3*dU3_dx-3*orho4*U3*U3*dU1_dx);
        dqx_dU2=tmpx*(-orho2*dU2_dx+2*orho3*U2*dU1_dx);
        dqx_dU3=tmpx*(-orho2*dU3_dx+2*orho3*U3*dU1_dx);
        dqx_dU4=-tmpx*orho2*dU1_dx;
        NekDouble tmpy=-ny*mu*gamma*oPr;
        dqy_dU1=tmpy*(-orho2*dU4_dy+2*orho3*U4*dU1_dy+2*orho3*U2*dU2_dy
                -3*orho4*U2*U2*dU1_dy+2*orho3*U3*dU3_dy-3*orho4*U3*U3*dU1_dy);
        dqy_dU2=tmpy*(-orho2*dU2_dy+2*orho3*U2*dU1_dy);
        dqy_dU3=tmpy*(-orho2*dU3_dy+2*orho3*U3*dU1_dy);
        dqy_dU4=-tmpy*orho2*dU1_dy;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]-dqx_dU1-dqy_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]-dqx_dU2-dqy_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]-dqx_dU3-dqy_dU3;
        tmpArray[2+3*nrow]=tmpArray[2+3*nrow]-dqx_dU4-dqy_dU4;
    }

     /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhow,rhoE]
     * Output: 4*5 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
     */
    void NavierStokesImplicitCFE::GetdFlux_dU_3D(
        const Array<OneD, NekDouble>                    &normals,
        const NekDouble                                 mu,
        const NekDouble                                 dmu_dT,
        const Array<OneD, NekDouble>                    &U,
        const Array<OneD, const Array<OneD, NekDouble>> &qfield,
              DNekMatSharedPtr                          &OutputMatrix)
    {
        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();

        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble U1=U[0];
        NekDouble U2=U[1];
        NekDouble U3=U[2];
        NekDouble U4=U[3];
        NekDouble U5=U[4];
        NekDouble dU1_dx=qfield[0][0];
        NekDouble dU2_dx=qfield[0][1];
        NekDouble dU3_dx=qfield[0][2];
        NekDouble dU4_dx=qfield[0][3];
        NekDouble dU5_dx=qfield[0][4];
        NekDouble dU1_dy=qfield[1][0];
        NekDouble dU2_dy=qfield[1][1];
        NekDouble dU3_dy=qfield[1][2];
        NekDouble dU4_dy=qfield[1][3];
        NekDouble dU5_dy=qfield[1][4];
        NekDouble dU1_dz=qfield[2][0];
        NekDouble dU2_dz=qfield[2][1];
        NekDouble dU3_dz=qfield[2][2];
        NekDouble dU4_dz=qfield[2][3];
        NekDouble dU5_dz=qfield[2][4];
        NekDouble gamma=m_gamma;
        NekDouble Cv=m_Cv;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;

        NekDouble orho1,orho2,orho3,orho4;
        NekDouble oCv=1./Cv;
        orho1=1.0/U1;
        orho2=orho1*orho1;
        orho3=orho1*orho2;
        orho4=orho2*orho2;

        //Assume Fn=mu*Sn
        //Sn=Sx*nx+Sy*ny+Sz*nz
        NekDouble TwoThrid=2./3.;
        NekDouble FourThird=2.0*TwoThrid;
        NekDouble tmp2=gamma*mu*oPr;
        NekDouble u=U2*orho1;
        NekDouble v=U3*orho1;
        NekDouble w=U4*orho1;
        NekDouble du_dx=orho1*(dU2_dx-u*dU1_dx);
        NekDouble dv_dx=orho1*(dU3_dx-v*dU1_dx);
        NekDouble dw_dx=orho1*(dU4_dx-w*dU1_dx);
        NekDouble du_dy=orho1*(dU2_dy-u*dU1_dy);
        NekDouble dv_dy=orho1*(dU3_dy-v*dU1_dy);
        NekDouble dw_dy=orho1*(dU4_dy-w*dU1_dy);
        NekDouble du_dz=orho1*(dU2_dz-u*dU1_dz);
        NekDouble dv_dz=orho1*(dU3_dz-v*dU1_dz);
        NekDouble dw_dz=orho1*(dU4_dz-w*dU1_dz);
        NekDouble s12=FourThird*du_dx-TwoThrid*dv_dy-TwoThrid*dw_dz;
        NekDouble s13=du_dy+dv_dx;
        NekDouble s14=dw_dx+du_dz;
        NekDouble s22=s13;
        NekDouble s23=FourThird*dv_dy-TwoThrid*du_dx-TwoThrid*dw_dz;
        NekDouble s24=dv_dz+dw_dy;
        NekDouble s32=s14;
        NekDouble s33=s24;
        NekDouble s34=FourThird*dw_dz-TwoThrid*du_dx-TwoThrid*dv_dy;
        NekDouble snx=s12*nx+s22*ny+s32*nz;
        NekDouble sny=s13*nx+s23*ny+s33*nz;
        NekDouble snz=s14*nz+s24*ny+s34*nz;
        NekDouble snv=snx*u+sny*v+snz*w;
        NekDouble qx=-tmp2*(orho1*dU5_dx-U5*orho2*dU1_dx-u*(orho1*dU2_dx
                    -U2*orho2*dU1_dx)-v*(orho1*dU3_dx-U3*orho2*dU1_dx)
                    -w*(orho1*dU4_dx-U4*orho2*dU1_dx));
        NekDouble qy=-tmp2*(orho1*dU5_dy-U5*orho2*dU1_dy-u*(orho1*dU2_dy
                    -U2*orho2*dU1_dy)-v*(orho1*dU3_dy-U3*orho2*dU1_dy)
                    -w*(orho1*dU4_dy-U4*orho2*dU1_dy));
        NekDouble qz=-tmp2*(orho1*dU5_dz-U5*orho2*dU1_dz-u*(orho1*dU2_dz
                    -U2*orho2*dU1_dz)-v*(orho1*dU3_dz-U3*orho2*dU1_dz)
                    -w*(orho1*dU4_dz-U4*orho2*dU1_dz));
        NekDouble qn=qx*nx+qy*ny+qz*nz;

        //Term1 mu's derivative with U: dmu_dU*Sn
        Array<OneD,NekDouble> tmp(4,0.0);
        tmp[0]=snx;
        tmp[1]=sny;
        tmp[2]=snz;
        tmp[3]=snv-qn/mu;
        Array<OneD,NekDouble> dT_dU (5,0.0);
        dT_dU[0]=oCv*(-orho2*U5+orho3*U2*U2+orho3*U3*U3+orho3*U4*U4);
        dT_dU[1]=-oCv*orho2*U2;
        dT_dU[2]=-oCv*orho2*U3;
        dT_dU[3]=-oCv*orho2*U4;
        dT_dU[4]=oCv*orho1;
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<5;j++)
            {
                tmpArray[i+j*nrow]=dmu_dT*dT_dU[j]*tmp[i];
            }
        }

        //Term 2 +mu*dSn_dU
        NekDouble du_dx_dU1,du_dx_dU2;
        NekDouble du_dy_dU1,du_dy_dU2;
        NekDouble du_dz_dU1,du_dz_dU2;
        NekDouble dv_dx_dU1,dv_dx_dU3;
        NekDouble dv_dy_dU1,dv_dy_dU3;
        NekDouble dv_dz_dU1,dv_dz_dU3;
        NekDouble dw_dx_dU1,dw_dx_dU4;
        NekDouble dw_dy_dU1,dw_dy_dU4;
        NekDouble dw_dz_dU1,dw_dz_dU4;
        NekDouble ds12_dU1,ds12_dU2,ds12_dU3,ds12_dU4;
        NekDouble ds13_dU1,ds13_dU2,ds13_dU3;
        NekDouble ds14_dU1,ds14_dU2,ds14_dU4;
        NekDouble ds22_dU1,ds22_dU2,ds22_dU3;
        NekDouble ds23_dU1,ds23_dU2,ds23_dU3,ds23_dU4;
        NekDouble ds24_dU1,ds24_dU3,ds24_dU4;
        NekDouble ds32_dU1,ds32_dU2,ds32_dU4;
        NekDouble ds33_dU1,ds33_dU3,ds33_dU4;
        NekDouble ds34_dU1,ds34_dU2,ds34_dU3,ds34_dU4;
        NekDouble dsnx_dU1,dsnx_dU2,dsnx_dU3,dsnx_dU4;
        NekDouble dsny_dU1,dsny_dU2,dsny_dU3,dsny_dU4;
        NekDouble dsnz_dU1,dsnz_dU2,dsnz_dU3,dsnz_dU4;
        NekDouble dsnv_dU1,dsnv_dU2,dsnv_dU3,dsnv_dU4;

        du_dx_dU1=-orho2*dU2_dx+2*orho3*U2*dU1_dx;
        du_dx_dU2=-orho2*dU1_dx;
        du_dy_dU1=-orho2*dU2_dy+2*orho3*U2*dU1_dy;
        du_dy_dU2=-orho2*dU1_dy;
        du_dz_dU1=-orho2*dU2_dz+2*orho3*U2*dU1_dz;
        du_dz_dU2=-orho2*dU1_dz;
        dv_dx_dU1=-orho2*dU3_dx+2*orho3*U3*dU1_dx;
        dv_dx_dU3=-orho2*dU1_dx;
        dv_dy_dU1=-orho2*dU3_dy+2*orho3*U3*dU1_dy;
        dv_dy_dU3=-orho2*dU1_dy;
        dv_dz_dU1=-orho2*dU3_dz+2*orho3*U3*dU1_dz;
        dv_dz_dU3=-orho2*dU1_dz;
        dw_dx_dU1=-orho2*dU4_dx+2*orho3*U4*dU1_dx;
        dw_dx_dU4=-orho2*dU1_dx;
        dw_dy_dU1=-orho2*dU4_dy+2*orho3*U4*dU1_dy;
        dw_dy_dU4=-orho2*dU1_dy;
        dw_dz_dU1=-orho2*dU4_dz+2*orho3*U4*dU1_dz;
        dw_dz_dU4=-orho2*dU1_dz;
        ds12_dU1=FourThird*du_dx_dU1-TwoThrid*dv_dy_dU1-TwoThrid*dw_dz_dU1;
        ds12_dU2=FourThird*du_dx_dU2;
        ds12_dU3=-TwoThrid*dv_dy_dU3;
        ds12_dU4=-TwoThrid*dw_dz_dU4;
        ds13_dU1=du_dy_dU1+dv_dx_dU1;
        ds13_dU2=du_dy_dU2;
        ds13_dU3=dv_dx_dU3;
        ds14_dU1=dw_dx_dU1+du_dz_dU1;
        ds14_dU2=du_dz_dU2;
        ds14_dU4=dw_dx_dU4;
        ds22_dU1=du_dy_dU1+dv_dx_dU1;
        ds22_dU2=du_dy_dU2;
        ds22_dU3=dv_dx_dU3;
        ds23_dU1=FourThird*dv_dy_dU1-TwoThrid*du_dx_dU1-TwoThrid*dw_dz_dU1;
        ds23_dU2=-TwoThrid*du_dx_dU2;
        ds23_dU3=FourThird*dv_dy_dU3;
        ds23_dU4=-TwoThrid*dw_dz_dU4;
        ds24_dU1=dv_dz_dU1+dw_dy_dU1;
        ds24_dU3=dv_dz_dU3;
        ds24_dU4=dw_dy_dU4;
        ds32_dU1=dw_dx_dU1+du_dz_dU1;
        ds32_dU2=du_dz_dU2;
        ds32_dU4=dw_dx_dU4;
        ds33_dU1=dv_dz_dU1+dw_dy_dU1;
        ds33_dU3=dv_dz_dU3;
        ds33_dU4=dw_dy_dU4;
        ds34_dU1=FourThird*dw_dz_dU1-TwoThrid*du_dx_dU1-TwoThrid*dv_dy_dU1;
        ds34_dU2=-TwoThrid*du_dx_dU2;
        ds34_dU3=-TwoThrid*dv_dy_dU3;
        ds34_dU4=FourThird*dw_dz_dU4;
        dsnx_dU1=ds12_dU1*nx+ds22_dU1*ny+ds32_dU1*nz;
        dsnx_dU2=ds12_dU2*nx+ds22_dU2*ny+ds32_dU2*nz;
        dsnx_dU3=ds12_dU3*nx+ds22_dU3*ny;
        dsnx_dU4=ds12_dU4*nx+ds32_dU4*nz;
        dsny_dU1=ds13_dU1*nx+ds23_dU1*ny+ds33_dU1*nz;
        dsny_dU2=ds13_dU2*nx+ds23_dU2*ny;
        dsny_dU3=ds13_dU3*nx+ds23_dU3*ny+ds33_dU3*nz;
        dsny_dU4=ds23_dU4*ny+ds33_dU4*nz;
        dsnz_dU1=ds14_dU1*nx+ds24_dU1*ny+ds34_dU1*nz;
        dsnz_dU2=ds14_dU2*nx+ds34_dU2*nz;
        dsnz_dU3=ds24_dU3*ny+ds34_dU3*nz;
        //? why there is value if 2D
        dsnz_dU4=ds14_dU4*nx+ds24_dU4*ny+ds34_dU4*nz;
        dsnv_dU1=u*dsnx_dU1+v*dsny_dU1+w*dsnz_dU1-orho2*U2*snx
                    -orho2*U3*sny-orho2*U4*snz;
        dsnv_dU2=u*dsnx_dU2+v*dsny_dU2+w*dsnz_dU2+orho1*snx;
        dsnv_dU3=u*dsnx_dU3+v*dsny_dU3+w*dsnz_dU3+orho1*sny;
        dsnv_dU4=u*dsnx_dU4+v*dsny_dU4+w*dsnz_dU4+orho1*snz;
        tmpArray[0+0*nrow]=tmpArray[0+0*nrow]+mu*dsnx_dU1;
        tmpArray[0+1*nrow]=tmpArray[0+1*nrow]+mu*dsnx_dU2;
        tmpArray[0+2*nrow]=tmpArray[0+2*nrow]+mu*dsnx_dU3;
        tmpArray[0+3*nrow]=tmpArray[0+3*nrow]+mu*dsnx_dU4;
        tmpArray[1+0*nrow]=tmpArray[1+0*nrow]+mu*dsny_dU1;
        tmpArray[1+1*nrow]=tmpArray[1+1*nrow]+mu*dsny_dU2;
        tmpArray[1+2*nrow]=tmpArray[1+2*nrow]+mu*dsny_dU3;
        tmpArray[1+3*nrow]=tmpArray[1+3*nrow]+mu*dsny_dU4;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]+mu*dsnz_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]+mu*dsnz_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]+mu*dsnz_dU3;
        tmpArray[2+3*nrow]=tmpArray[2+3*nrow]+mu*dsnz_dU4;
        tmpArray[3+0*nrow]=tmpArray[3+0*nrow]+mu*dsnv_dU1;
        tmpArray[3+1*nrow]=tmpArray[3+1*nrow]+mu*dsnv_dU2;
        tmpArray[3+2*nrow]=tmpArray[3+2*nrow]+mu*dsnv_dU3;
        tmpArray[3+3*nrow]=tmpArray[3+3*nrow]+mu*dsnv_dU4;

        //Consider heat flux qn's effect (does not include mu's effect)
        NekDouble dqx_dU1,dqx_dU2,dqx_dU3,dqx_dU4,dqx_dU5;
        NekDouble dqy_dU1,dqy_dU2,dqy_dU3,dqy_dU4,dqy_dU5;
        NekDouble dqz_dU1,dqz_dU2,dqz_dU3,dqz_dU4,dqz_dU5;
        NekDouble tmpx=-nx*tmp2;
        dqx_dU1=tmpx*(-orho2*dU5_dx+2*orho3*U5*dU1_dx+2*orho3*U2*dU2_dx
                -3*orho4*U2*U2*dU1_dx+2*orho3*U3*dU3_dx-3*orho4*U3*U3*dU1_dx
                +2*orho3*U4*dU4_dx-3*orho4*U4*U4*dU1_dx);
        dqx_dU2=tmpx*(-orho2*dU2_dx+2*orho3*U2*dU1_dx);
        dqx_dU3=tmpx*(-orho2*dU3_dx+2*orho3*U3*dU1_dx);
        dqx_dU4=tmpx*(-orho2*dU4_dx+2*orho3*U4*dU1_dx);
        dqx_dU5=-tmpx*orho2*dU1_dx;
        NekDouble tmpy=-ny*tmp2;
        dqy_dU1=tmpy*(-orho2*dU5_dy+2*orho3*U5*dU1_dy+2*orho3*U2*dU2_dy
                -3*orho4*U2*U2*dU1_dy+2*orho3*U3*dU3_dy-3*orho4*U3*U3*dU1_dy
                +2*orho3*U4*dU4_dy-3*orho4*U4*U4*dU1_dy);
        dqy_dU2=tmpy*(-orho2*dU2_dy+2*orho3*U2*dU1_dy);
        dqy_dU3=tmpy*(-orho2*dU3_dy+2*orho3*U3*dU1_dy);
        dqy_dU4=tmpy*(-orho2*dU4_dy+2*orho3*U4*dU1_dy);
        dqy_dU5=-tmpy*orho2*dU1_dy;
        NekDouble tmpz=-nz*tmp2;
        dqz_dU1=tmpz*(-orho2*dU5_dz+2*orho3*U5*dU1_dz+2*orho3*U2*dU2_dz
                -3*orho4*U2*U2*dU1_dz+2*orho3*U3*dU3_dz-3*orho4*U3*U3*dU1_dz
                +2*orho3*U4*dU4_dz-3*orho4*U4*U4*dU1_dz);
        dqz_dU2=tmpz*(-orho2*dU2_dz+2*orho3*U2*dU1_dz);
        dqz_dU3=tmpz*(-orho2*dU3_dz+2*orho3*U3*dU1_dz);
        dqz_dU4=tmpz*(-orho2*dU4_dz+2*orho3*U4*dU1_dz);
        dqz_dU5=-tmpz*orho2*dU1_dz;
        tmpArray[3+0*nrow]=tmpArray[3+0*nrow]-dqx_dU1-dqy_dU1-dqz_dU1;
        tmpArray[3+1*nrow]=tmpArray[3+1*nrow]-dqx_dU2-dqy_dU2-dqz_dU2;
        tmpArray[3+2*nrow]=tmpArray[3+2*nrow]-dqx_dU3-dqy_dU3-dqz_dU3;
        tmpArray[3+3*nrow]=tmpArray[3+3*nrow]-dqx_dU4-dqy_dU4-dqz_dU4;
        tmpArray[3+4*nrow]=tmpArray[3+4*nrow]-dqx_dU5-dqy_dU5-dqz_dU5;
    }

    void NavierStokesImplicitCFE::v_MinusDiffusionFluxJacPoint(
        const int                                       nConvectiveFields,
        const int                                       nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const TensorOfArray3D<NekDouble>                &locDerv,
        const Array<OneD, NekDouble>                    &locmu,
        const Array<OneD, NekDouble>                    &locDmuDT,
        const Array<OneD, NekDouble>                    &normals,
              DNekMatSharedPtr                          &wspMat,
              Array<OneD,       Array<OneD, NekDouble>> &PntJacArray)
    {
        int nSpaceDim           = m_graph->GetSpaceDimension();  

        NekDouble pointmu       = 0.0;
        NekDouble pointDmuDT    = 0.0;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > pointDerv(nSpaceDim);
        for(int j = 0; j < nSpaceDim; j++)
        {   
            pointDerv[j] = Array<OneD, NekDouble>(nConvectiveFields,0.0);
        }

        Array<OneD, NekDouble > wspMatData = wspMat->GetPtr();
        Array<OneD, NekDouble > tmp1;
        Array<OneD, NekDouble > tmp2;
        Array<OneD, NekDouble > tmp3;

        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }
            for(int j = 0; j < nSpaceDim; j++)
            {   
                for(int k = 0; k < nConvectiveFields; k++)
                {
                    pointDerv[j][k] = locDerv[j][k][npnt];
                }
            }

            pointmu     = locmu[npnt];
            pointDmuDT  = locDmuDT[npnt];

            v_GetDiffusionFluxJacPoint(pointVar,pointDerv,pointmu,pointDmuDT,
                normals,wspMat);
            for (int j =0; j < nConvectiveFields; j++)
            {
                int noffset = j*nConvectiveFields;

                Vmath::Vsub(nConvectiveFields-1, 
                    tmp1 = PntJacArray[npnt] + (noffset+1),1,
                    tmp2 = wspMatData + (noffset-j),1,
                    tmp3 = PntJacArray[npnt] + (noffset+1),1);
            }
        }
    }

    void NavierStokesImplicitCFE::v_GetDiffusionFluxJacPoint(
        const Array<OneD, NekDouble>                    &conservVar, 
        const Array<OneD, const Array<OneD, NekDouble>> &conseDeriv, 
        const NekDouble                                 mu,
        const NekDouble                                 DmuDT,
        const Array<OneD, NekDouble>                    &normals,
              DNekMatSharedPtr                          &fluxJac)
    {
        switch (m_spacedim)
        {
        case 2:
            GetdFlux_dU_2D(normals,mu,DmuDT,conservVar,conseDeriv,fluxJac);
            break;

        case 3:
            GetdFlux_dU_3D(normals,mu,DmuDT,conservVar,conseDeriv,fluxJac);
            break;

        default:
            NEKERROR(ErrorUtil::efatal, "v_GetDiffusionFluxJacPoint not coded");
            break;
        }
    }

    void NavierStokesImplicitCFE::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr            &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              TensorOfArray5D<NekDouble>                &ElmtJacArray,
        const int                                       nFluxDir)
    {
        int nConvectiveFields   = inarray.size();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            explist->GetExp();
        int nTotElmt            = (*expvect).size();
        int nPts                = explist->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu[0], mu, 1);
        }

        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);
        Array<OneD, Array<OneD, NekDouble> > locnormal(nSpaceDim);

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);
        Array<OneD, NekDouble > PointFJac_data = PointFJac->GetPtr();

        for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = explist->GetPhys_Offset(nelmt);

            for(int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+noffest;
            }

            for(int j = 0; j < nSpaceDim; j++)
            {
                locnormal[j] = normals[j]+noffest;
            }

            locmu       =   mu      + noffest;
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for(int j = 0; j < nSpaceDim; j++)
                {
                    pointnormals[j] = locnormal[j][npnt];
                }

                pointmu     = locmu[npnt];

                m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,
                    pointVar,PointFJac);
                for (int j =0; j < nConvectiveFields; j++)
                {
                    ElmtJacArray[0][j][nFluxDir][nelmt][npnt] = 0.0;
                }
                for (int j =0; j < nConvectiveFields; j++)
                {
                    int noffset = j*(nConvectiveFields-1);
                    for (int i =0; i < nConvectiveFields-1; i++)
                    {
                        ElmtJacArray[i+1][j][nFluxDir][nelmt][npnt] = 
                            PointFJac_data[noffset+i];
                    }
                }
            }
        }
    }

    void NavierStokesImplicitCFE::v_GetFluxDerivJacDirctnElmt(
        const int                                       nConvectiveFields,
        const int                                       nElmtPnt,
        const int                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble>                    &locmu,
        const Array<OneD, const Array<OneD, NekDouble>> &locnormal,
              DNekMatSharedPtr                          &wspMat,
              Array<OneD,       Array<OneD, NekDouble>> &PntJacArray)
    {
        int nSpaceDim           = m_graph->GetSpaceDimension();  
        
        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);

        Array<OneD, NekDouble > wspMatData = wspMat->GetPtr();

        Array<OneD, NekDouble > tmp1;
        Array<OneD, NekDouble > tmp2;
                
        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }
            for(int j = 0; j < nSpaceDim; j++)
            {   
                pointnormals[j] = locnormal[j][npnt];
            }

            pointmu     = locmu[npnt];

            m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,
                pointVar,wspMat);
            Vmath::Zero(nConvectiveFields,PntJacArray[npnt],nConvectiveFields);
            for (int j =0; j < nConvectiveFields; j++)
            {
                int noffset = j*(nConvectiveFields-1);
                Vmath::Vcopy((nConvectiveFields-1), 
                    tmp1 = wspMatData + noffset,1, 
                    tmp2 = PntJacArray[npnt] + noffset+j+1,1);
            }
        }
    }
    
    void NavierStokesImplicitCFE::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                   &explist,
        const Array<OneD, const Array<OneD, NekDouble>>        &normals,
        const int                                              nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>>        &inarray,
              Array<OneD,       Array<OneD, DNekMatSharedPtr>> &ElmtJac)
    {
        int nConvectiveFields   = inarray.size();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =
            explist->GetExp();
        int nTotElmt            = (*expvect).size();
        int nPts                = explist->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();

        //Debug
        if(!ElmtJac.size())
        {
            ElmtJac = Array<OneD, Array<OneD, DNekMatSharedPtr> > (nTotElmt);
            for(int  nelmt = 0; nelmt < nTotElmt; nelmt++)
            {
                int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
                ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nElmtPnt);
                for(int npnt = 0; npnt < nElmtPnt; npnt++)
                {
                    ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nConvectiveFields, 
                        nConvectiveFields);
                }
            }
        }
        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, m_mu[0]);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            Array<OneD, NekDouble> muav;
            if (m_fields[0]->GetTrace()->GetTotPoints()==nPts)
            {
                muav = m_muavTrace;
            }
            else
            {
                muav = m_muav;
            }
            Vmath::Vadd(nPts, mu, 1, muav, 1, mu, 1);
        }

        // What about thermal conductivity?

        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);
        Array<OneD, Array<OneD, NekDouble> > locnormal(nSpaceDim);

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);
        Array<OneD, NekDouble > tmpMatinnData, tmpMatoutData;
        Array<OneD, NekDouble > tmp1, tmp2;
        
        for (int nelmt = 0; nelmt < nTotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = explist->GetPhys_Offset(nelmt);
            
            for (int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+noffest;
            }
            
            for (int j = 0; j < nSpaceDim; j++)
            {
                locnormal[j] = normals[j]+noffest;
            }
            
            locmu       =   mu      + noffest;
            for (int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for (int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for (int j = 0; j < nSpaceDim; j++)
                {
                    pointnormals[j] = locnormal[j][npnt];
                }
                
                pointmu     = locmu[npnt];
                
                m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,
                    pointVar,PointFJac);
                tmpMatinnData = PointFJac->GetPtr();
                tmpMatoutData = ElmtJac[nelmt][npnt]->GetPtr();
                
                Vmath::Fill(nConvectiveFields,0.0,tmpMatoutData,
                            nConvectiveFields);
                for (int j =0; j < nConvectiveFields; j++)
                {
                    Vmath::Vcopy(nConvectiveFields-1,
                        tmp1 = tmpMatinnData + (j*(nConvectiveFields-1)),1,
                        tmp2 = tmpMatoutData + (1+j*nConvectiveFields),1);
                }
            }
        }
    }
    
    void NavierStokesImplicitCFE::v_CalcPhysDeriv(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              TensorOfArray3D<NekDouble>                &qfield)
    {
        int nConvectiveFields = m_fields.size();
        int npoints           = GetTotPoints();
        const Array<OneD, Array<OneD, NekDouble>> pFwd;
        const Array<OneD, Array<OneD, NekDouble>> pBwd;
        if(!qfield.size())
        {
            qfield  = TensorOfArray3D<NekDouble> (m_spacedim);
            for(int i = 0; i< m_spacedim; i++)
            {
                qfield[i] = Array<OneD, Array<OneD, NekDouble>>(nConvectiveFields);
                for(int j = 0; j< nConvectiveFields; j++)
                {
                    qfield[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
        }
        m_diffusion->DiffuseCalcDerivative(m_fields,inarray,qfield,
            pFwd,pBwd);
    }

    void NavierStokesImplicitCFE::v_CalcMuDmuDT(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD, NekDouble>                    &mu,
              Array<OneD, NekDouble>                    &DmuDT)
    {
        int npoints = mu.size();
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature (npoints, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
            if(DmuDT.size()>0)
            {
                m_varConv->GetDmuDT(temperature,mu,DmuDT);
            }
        }
        else
        {
            Vmath::Vcopy(npoints, m_mu, 1, mu, 1);
            if(DmuDT.size()>0)
            {
                Vmath::Zero(npoints, DmuDT, 1);
            }
        }
    }    
}
