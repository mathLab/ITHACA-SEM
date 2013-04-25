///////////////////////////////////////////////////////////////////////////////
//
// File ADR2DManifold.cpp
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
// Description: ADR2Dm class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <ADR2DManifoldSolver/ADR2DManifold.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#if defined(__INTEL_COMPILER)
    #include <mathimf.h>
#else
    #include <cmath>
#endif

namespace Nektar
{
    /**
     * Basic construnctor
     */
    //    ADR2DManifold::ADR2DManifold(void):
    //        EquationSystem(),
    //        m_infosteps(100),
    //        m_explicitAdvection(true),
    //        m_explicitDiffusion(true),
    //        m_explicitReaction(true)
    //    {
    //    }

    /**
     * Constructor. Creates ...
     */
    ADR2DManifold::ADR2DManifold(
            LibUtilities::SessionReaderSharedPtr& pSession):
            EquationSystem(pSession),
            m_infosteps(10),
            m_explicitDiffusion(true),
            m_explicitReaction(true)
    {
        int i;

        m_epsilon = Array<OneD, NekDouble>(2,0.0);

        m_session->LoadParameter("epsilon0",       m_epsilon[0], 0.0);
        m_session->LoadParameter("epsilon1",       m_epsilon[1], 0.0);
        m_session->LoadParameter("Beta",           m_beta, 0.0);
        m_session->LoadParameter("UseDirDeriv",    m_UseDirDeriv, 0);

        m_session->LoadParameter("k",              mK,   0.0);
        m_session->LoadParameter("a",              mA,   0.0);
        m_session->LoadParameter("b",              mB,   0.0);
        m_session->LoadParameter("c",              mC,   0.0);
        m_session->LoadParameter("eps",            mEps, 0.0);
        m_session->LoadParameter("mu1",            mMu1, 0.0);
        m_session->LoadParameter("mu2",            mMu2, 0.0);

        m_session->LoadParameter("diam",           mDiam,1.0);

        // -2: Advection problem on sphere
        // -1: Linear Morphogenesis problem for diffusion-reaction test
        // 0: Planar propagation from -x direction
        // 1: Planar propagation from -y direction
        // 2: For Planar propagation from -x direction in xy plane
        // 10: For Planar propagation from +z direction
        m_session->LoadParameter("initialwavetype",m_initialwavetype, 0);
        m_session->LoadParameter("duration",       m_duration, 3.5);
        m_session->LoadParameter("x0c",            m_x0c, 0.0);
        m_session->LoadParameter("x1c",            m_x1c, 0.0);
        m_session->LoadParameter("x2c",            m_x2c, 0.0);
        m_session->LoadParameter("Anisotropy",     m_Anisotropy, 0);
        m_session->LoadParameter("Angularfreq",    m_Angularfreq,
                3.14159265358979323846/6.0);
        m_session->LoadParameter("Angleofaxis",    m_Angleofaxis, 0.0);

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_session->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(NoCaseStringCompare(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i;
                break;
            }
        }

        ASSERTL0(i != (int) eEquationTypeSize, "Invalid expansion type.");

        // Equation specific Setups
        int nqtot = GetNpoints();
        switch(m_equationType)
        {
        case eUnsteadyAdvection:
        {
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (nqtot);
            }

            EvaluateAdvectionVelocity();

            // m_velmagnitude = || v ||---------------------------------
            m_velmagnitude = Array<OneD, Array<OneD, NekDouble> >(1);

            m_velmagnitude[0] = Array<OneD, NekDouble> (nqtot, 0.0);
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nqtot, &m_velocity[i][0], 1, &m_velocity[i][0], 1,
                        &m_velmagnitude[0][0], 1, &m_velmagnitude[0][0], 1);
            }
            Vmath::Vsqrt(nqtot, &m_velmagnitude[0][0], 1, &m_velmagnitude[0][0], 1);
            // ------------------------------------------------------------------

            // m_vellc = \nabla \cdot m_velocity ---------------------------------
            NekDouble vellcavg;
            m_vellc = Array<OneD, Array<OneD, NekDouble> >(1);
            m_vellc[0] = Array<OneD, NekDouble>(nqtot,0.0);
            Array<OneD, NekDouble> temp0(nqtot), temp1(nqtot), temp2(nqtot);

            m_fields[0]->PhysDeriv(m_velocity[0], temp0, temp1, temp2);
            Vmath::Vadd(nqtot, &temp0[0], 1, &m_vellc[0][0], 1, &m_vellc[0][0], 1);

            m_fields[0]->PhysDeriv(m_velocity[1], temp0, temp1, temp2);
            Vmath::Vadd(nqtot, &temp1[0], 1, &m_vellc[0][0], 1, &m_vellc[0][0], 1);

            m_fields[0]->PhysDeriv(m_velocity[2], temp0, temp1, temp2);
            Vmath::Vadd(nqtot, &temp2[0], 1, &m_vellc[0][0], 1, &m_vellc[0][0], 1);

            vellcavg = Vmath::Vsum(nqtot,m_vellc[0],1);
            vellcavg = vellcavg/nqtot;

            cout << "Average gradient of velocity = " << vellcavg << endl;
            // --------------------------------------------------------------------
            m_dirForcing = Array<OneD, Array<OneD, NekDouble> >(2);
            m_dirForcing[0] = Array<OneD, NekDouble>(3*nqtot,0.0);
            m_dirForcing[1] = Array<OneD, NekDouble>(3*nqtot,0.0);

            NekDouble mag ;
            for (int i=0; i < m_spacedim; ++i)
            {
                Vmath::Vcopy(nqtot, &m_velocity[i][0], 1, &m_dirForcing[0][i*nqtot], 1);
                for(int k=0; k < nqtot; ++k)
                {
                    mag = m_velmagnitude[0][k];
                    if(fabs(mag)>0.000001)
                    {
                        m_dirForcing[0][i*nqtot+k] = m_dirForcing[0][i*nqtot+k]/mag;
                    }
                }
            }

            int nTraceNumPoints = GetTraceNpoints();
            Array<OneD, NekDouble> Fwd(nTraceNumPoints);
            Array<OneD, NekDouble> temp(nqtot);
            m_Vn = Array<OneD, NekDouble> (nTraceNumPoints,0.0);

            // m_Vn = m_dirForcing \cdot m_traceNormals when using direction derivative
                    //      = m_velocity \cdot m_traceNormals otherwise
            for(int i = 0; i < m_spacedim; ++i)
            {
                if(m_UseDirDeriv)
                {
                    Vmath::Vcopy(nqtot, &m_dirForcing[0][i*nqtot], 1, &temp[0], 1);
                }

                else
                {
                    Vmath::Vcopy(nqtot, &m_velocity[i][0], 1, &temp[0], 1);
                }

                m_fields[0]->ExtractTracePhys(temp, Fwd);
                Vmath::Vvtvp(nTraceNumPoints,&m_traceNormals[i][0],1,&Fwd[0],1,&m_Vn[0],1,&m_Vn[0],1);
            }

            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;
            goto UnsteadySetup;
            break;
        }
        case eUnsteadyDiffusion:
        {
            // Evaluate two tangential vectors on 2D manifold
            if(m_spacedim==3)
            {
                SetUpTangentialVectors();
            }

            if(GetExplicitDiffusion())
            {
                m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;
            }

            else
            {
                m_timeIntMethod = LibUtilities::eDIRKOrder3;
            }
            goto UnsteadySetup;
            break;
        }
        case eUnsteadyDiffusionReaction:
        case eNonlinearMorphogenensis:
        case eFHNtesttype1:
        case eFHNMONO:
        case eAlievPanfilov:
        {
            // Evaluate two tangential vectors on 2D manifold
            if(m_spacedim==3)
            {
                SetUpTangentialVectors();
            }

            m_timeIntMethod = LibUtilities::eIMEXdirk_3_4_3;

            goto UnsteadySetup;
            break;
        }
        UnsteadySetup:
    {
        std::string Implicit = "Implicit";
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);

        // check that any user defined boundary condition is indeed implemented
        for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Time Dependent Boundary Condition (if no use defined then this is empty)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eNoUserDefined)
            {
                if  ( (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eTimeDependent &&
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eMG   ) &&
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eWALL  )

                {
                    ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                }
            }
        }

        // Check for definition of Implicit/Explicit terms in solverinfo
        if(m_session->DefinesSolverInfo("ADVECTIONADVANCEMENT"))
        {
            std::string AdvStr = m_session->GetSolverInfo("ADVECTIONADVANCEMENT");

            if(NoCaseStringCompare(AdvStr,Implicit) == 0)
            {
                m_explicitAdvection = false;
            }
            else
            {
                m_explicitAdvection = true;
            }
        }
        else
        {
            m_explicitAdvection = true;
        }

        if(m_session->DefinesSolverInfo("DIFFUSIONADVANCEMENT"))
        {
            std::string AdvStr = m_session->GetSolverInfo("DIFFUSIONADVANCEMENT");

            if(NoCaseStringCompare(AdvStr,Implicit) == 0 )
            {
                m_explicitDiffusion = false;
                // Reset default for implicit diffusion
                if(m_equationType == eUnsteadyDiffusion)
                {
                    m_timeIntMethod = LibUtilities::eDIRKOrder3;
                }
            }
            else
            {
                m_explicitDiffusion = true;
            }
        }
        else
        {
            m_explicitDiffusion = true;
        }

        if(m_session->DefinesSolverInfo("REACTIONADVANCEMENT"))
        {
            std::string AdvStr = m_session->GetSolverInfo("REACTIONADVANCEMENT");

            if(NoCaseStringCompare(AdvStr,Implicit) == 0)
            {
                m_explicitReaction = false;
            }
            else
            {
                m_explicitReaction = true;
            }
        }
        else
        {
            m_explicitReaction = true;
        }

        // check to see if time stepping has been reset
        if(m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHE"))
        {
            std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
            int i;
            for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
            {
                if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
                {
                    m_timeIntMethod = (LibUtilities::TimeIntegrationMethod)i;
                    break;
                }
            }

            ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
        }

        break;
    }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }

        // Plot Tangential vector map ==============================
        if(m_fields.num_elements()==6)
        {
            PlotTangentialVectorMap();
        }
    }

    void ADR2DManifold::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        std::string velStr[3] = {"Vx","Vy","Vz"};

        // Get the coordinates (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        if(m_session->DefinesFunction("AdvectionVelocity"))
        {
            for(int i = 0 ; i < m_velocity.num_elements(); i++)
            {
                LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("AdvectionVelocity", velStr[i]);
                ifunc->Evaluate(x0,x1,x2,m_velocity[i]);
            }
        }

        else
        {
            NekDouble vel_phi, vel_theta, radius, sin_phi, cos_phi, sin_theta, cos_theta;
            NekDouble Tol = 0.00001;

            // theta = a*sin(z/r),  phi = a*tan(y/x);
            for (int j = 0; j < nq; j++)
            {
                radius = sqrt( x0[j]*x0[j] + x1[j]*x1[j] + x2[j]*x2[j] );

                // At North and South poles, (ax,ay,ax) = (0, Omega_f*sin(alpha),0)
                if( fabs(fabs(x2[j]) - radius) < Tol )
                {
                    sin_theta = x2[j]/radius;

                    m_velocity[0][j] = 0.0;
                    m_velocity[1][j] = m_Angularfreq*sin_theta*sin(m_Angleofaxis);
                    m_velocity[2][j] = 0.0;
                }
                else
                {
                    sin_phi = x1[j]/sqrt( x0[j]*x0[j] + x1[j]*x1[j] );
                    cos_phi = x0[j]/sqrt( x0[j]*x0[j] + x1[j]*x1[j] );

                    sin_theta = x2[j]/radius;
                    cos_theta = sqrt( x0[j]*x0[j] + x1[j]*x1[j] )/radius;

                    vel_phi = m_Angularfreq*(cos_theta*cos(m_Angleofaxis) + sin_theta*cos_phi*sin(m_Angleofaxis) );
                    vel_theta = -1.0*m_Angularfreq*sin_phi*sin(m_Angleofaxis);

                    m_velocity[0][j] = -1.0*vel_phi*sin_phi - vel_theta*sin_theta*cos_phi ;
                    m_velocity[1][j] = vel_phi*cos_phi - vel_theta*sin_theta*sin_phi ;
                    m_velocity[2][j] = vel_theta*cos_theta;
                }
            }
        }
    }

    void ADR2DManifold::ODErhs(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble t)
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
        {
            switch(m_equationType)
            {
            case eUnsteadyAdvection:
            {
                Array<OneD, Array<OneD, NekDouble> > temp(nvariables);

                SetBoundaryConditions(inarray,t);

                if(m_UseDirDeriv)
                {
                    WeakDGAdvectionDir(inarray, m_dirForcing, outarray);

                    for (i = 0; i < nvariables; ++i)
                    {
                        temp[i] = Array<OneD, NekDouble> (ncoeffs,0.0);

                        StdRegions::VarCoeffMap varcoeffs;
                        varcoeffs[StdRegions::eVarCoeffMass] = m_vellc[0];
                        MultiRegions::GlobalMatrixKey wkey(StdRegions::eMass,
                                MultiRegions::NullAssemblyMapSharedPtr,
                                StdRegions::NullConstFactorMap,
                                varcoeffs);
                        m_fields[0]->MultiRegions::ExpList::GeneralMatrixOp(wkey, inarray[i], temp[i]);
                    }
                }
                else
                {
                    WeakDGAdvection(inarray, outarray);
                }

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);

                    if(m_UseDirDeriv)
                    {
                        StdRegions::VarCoeffMap varcoeffs;
                        varcoeffs[StdRegions::eVarCoeffMass] = m_velmagnitude[0];
                        MultiRegions::GlobalMatrixKey key(StdRegions::eMass,
                                MultiRegions::NullAssemblyMapSharedPtr,
                                StdRegions::NullConstFactorMap,
                                varcoeffs);
                        m_fields[0]->MultiRegions::ExpList::GeneralMatrixOp(key, outarray[i], outarray[i]);

                        Vmath::Vadd(ncoeffs, temp[i], 1, outarray[i], 1, outarray[i], 1);

                        m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
                    }
                    Vmath::Neg(ncoeffs,outarray[i],1);
                }
                break;
            }
            case eUnsteadyDiffusion:
            {
                // Apply laplacian operator
                WeakDGDiffusion(inarray,outarray);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
                }
                break;
            }
            case eUnsteadyDiffusionReaction:
            case eNonlinearMorphogenensis:
            case eFHNtesttype1:
            case eFHNMONO:
            case eAlievPanfilov:
            {
                // Apply laplacian operator
                WeakDGDiffusion(inarray,outarray);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
                }

                Array<OneD, NekDouble> Reaction1(ncoeffs,0.0);
                Array<OneD, NekDouble> Reaction2(ncoeffs,0.0);

                // Mophogenesis reaction diffusion equation
                Vmath::Svtvp(ncoeffs,m_a,&inarray[0][0],1,&Reaction1[0],1,&Reaction1[0],1);
                Vmath::Svtvp(ncoeffs,m_b,&inarray[1][0],1,&Reaction1[0],1,&Reaction1[0],1);

                Vmath::Svtvp(ncoeffs,m_epsilon[0],&outarray[0][0],1,&Reaction1[0],1,&outarray[0][0],1);

                Vmath::Svtvp(ncoeffs,m_c,&inarray[0][0],1,&Reaction2[0],1,&Reaction2[0],1);
                Vmath::Svtvp(ncoeffs,m_d,&inarray[1][0],1,&Reaction2[0],1,&Reaction2[0],1);

                Vmath::Svtvp(ncoeffs,m_epsilon[1],&outarray[1][0],1,&Reaction2[0],1,&outarray[1][0],1);
                break;
            }
            default:
                ASSERTL0(false, "Unknown equation type.");
                break;
            }
            break;
        }
        default:
        {
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
        }
    }

    void ADR2DManifold::ODEeLinearMGReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        int ncoeffs    = inarray[0].num_elements();

        Array<OneD, NekDouble> temp(ncoeffs);

        temp = Array<OneD, NekDouble>(ncoeffs,0.0);
        Vmath::Svtvp(ncoeffs,m_a,&inarray[0][0],1,&temp[0],1,&temp[0],1);
        Vmath::Svtvp(ncoeffs,m_b,&inarray[1][0],1,&temp[0],1,&outarray[0][0],1);

        temp = Array<OneD, NekDouble>(ncoeffs,0.0);
        Vmath::Svtvp(ncoeffs,m_c,&inarray[0][0],1,&temp[0],1,&temp[0],1);
        Vmath::Svtvp(ncoeffs,m_d,&inarray[1][0],1,&temp[0],1,&outarray[1][0],1);
    }

    void ADR2DManifold::ODEeReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        const NekDouble coeff = 3.14159265*3.14159265;

        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Smul(ncoeffs, coeff, inarray[i], 1, outarray[i], 1);
        }
    }

    void ADR2DManifold::ODEeNonlinearMorphoReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        int npoints    = m_fields[0]->GetNpoints();

        NekDouble A = 2.0;
        NekDouble B = 5.0;

        Array<OneD, NekDouble> phys0(npoints);
        Array<OneD, NekDouble> phys1(npoints);
        Array<OneD, NekDouble> phys0out(npoints);
        Array<OneD, NekDouble> phys1out(npoints);
        Array<OneD, NekDouble> cubterm(npoints,0.0);

        // Back to Physicsfield
        m_fields[0]->BwdTrans(inarray[0],phys0);
        m_fields[0]->SetPhysState(true);

        m_fields[1]->BwdTrans(inarray[1],phys1);
        m_fields[1]->SetPhysState(true);

        // temp = phys0*phys0*phy1
        Vmath::Vmul(npoints,&phys0[0],1,&phys0[0],1,&cubterm[0],1);
        Vmath::Vmul(npoints,&phys1[0],1,&cubterm[0],1,&cubterm[0],1);
        
                // outarray[0] = A - B*phy0 + phy0*phy0*phy1 - phy0
        NekDouble coeff = -1.0*(B+1.0);
        Vmath::Svtvp(npoints,coeff,&phys0[0],1,&cubterm[0],1,&phys0out[0],1);
        Vmath::Sadd(npoints,A,&phys0out[0],1,&phys0out[0],1);
        
        // outarray[1] = B*phys1 - phy0*phy0*phy1
        Vmath::Svtsvtp(npoints,B,&phys1[0],1,-1.0,&cubterm[0],1,&phys1out[0],1);
        
        m_fields[0]->FwdTrans(phys0out, outarray[0]);
        m_fields[0]->SetPhysState(false);
        
        m_fields[1]->FwdTrans(phys1out, outarray[1]);
        m_fields[1]->SetPhysState(false);
    }

    void ADR2DManifold::ODEeReactionFHNtest1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        int nvariables = inarray.num_elements();
        int npoints = m_fields[0]->GetNpoints();

        NekDouble epsilon = 1.0/32.0;

        const NekDouble coeff = 2.0/epsilon;

        Array<OneD, NekDouble> physfield(npoints);
        Array<OneD, NekDouble> utemp2(npoints,0.0);
        Array<OneD, NekDouble> utemp3(npoints,0.0);

        for (int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->BwdTrans(inarray[i],physfield);
            m_fields[i]->SetPhysState(true);

            // temp3 = (2/epsilon)*(u*u - u*u*u)
            Vmath::Vmul(npoints, &physfield[0], 1, &physfield[0], 1, &utemp2[0], 1);
            Vmath::Vmul(npoints, &physfield[0], 1, &utemp2[0], 1, &utemp3[0], 1);
            Vmath::Vsub(npoints, &utemp2[0], 1, &utemp3[0], 1, &physfield[0], 1);
            Vmath::Smul(npoints, coeff, &physfield[0], 1, &physfield[0], 1);

            m_fields[i]->FwdTrans(physfield,outarray[i]);
            m_fields[i]->SetPhysState(false);
        }
    }

    void ADR2DManifold::ODEeReactionFHNmono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        NekDouble m_gamma = 0.5;

        int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> physfieldu(npoints);
        Array<OneD, NekDouble> physfieldv(npoints);

        Array<OneD, NekDouble> Ru(npoints,0.0);
        Array<OneD, NekDouble> Rv(npoints, 0.0);
        Array<OneD, NekDouble> u3(npoints,0.0);

        m_fields[0]->BwdTrans(inarray[0],physfieldu);
        m_fields[0]->SetPhysState(true);

        m_fields[1]->BwdTrans(inarray[1],physfieldv);
        m_fields[1]->SetPhysState(true);

        // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
        // physfield = u - (1.0/3.0)*u*u*u
        Vmath::Vmul(npoints, &physfieldu[0], 1, &physfieldu[0], 1, &Ru[0], 1);
        Vmath::Vmul(npoints, &physfieldu[0], 1, &Ru[0], 1, &u3[0], 1);
        Vmath::Svtvp(npoints, (-1.0/3.0), &u3[0], 1, &physfieldu[0], 1, &Ru[0], 1);

        Vmath::Vsub(npoints, &physfieldv[0], 1, &Ru[0], 1, &Ru[0], 1);
        Vmath::Smul(npoints, -1.0/m_epsilon[0], &Ru[0], 1, &Ru[0], 1);

        m_fields[0]->FwdTrans(Ru,outarray[0]);
        m_fields[0]->SetPhysState(false);

        // For v: m_epsilon*( u + m_beta - m_gamma*v )
        Vmath::Svtvp(npoints, -1.0*m_gamma, &physfieldv[0], 1, &physfieldu[0], 1, &Rv[0], 1);
        Vmath::Sadd(npoints, m_beta, &Rv[0], 1, &Rv[0], 1);
        Vmath::Smul(npoints, m_epsilon[0], &Rv[0], 1, &Rv[0], 1);

        m_fields[1]->FwdTrans(Rv,outarray[1]);
        m_fields[1]->SetPhysState(false);
    }

    void ADR2DManifold::ODEeReactionAP(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> physfieldu(npoints);
        Array<OneD, NekDouble> physfieldv(npoints);

        Array<OneD, NekDouble> Ru(npoints,0.0);
        Array<OneD, NekDouble> Rv(npoints, 0.0);
        Array<OneD, NekDouble> u2(npoints,0.0);
        Array<OneD, NekDouble> u3(npoints,0.0);

        m_fields[0]->BwdTrans(inarray[0],physfieldu);
        m_fields[0]->SetPhysState(true);

        m_fields[1]->BwdTrans(inarray[1],physfieldv);
        m_fields[1]->SetPhysState(true);

        // u2 = u*u
                Vmath::Vmul(npoints, &physfieldu[0], 1, &physfieldu[0], 1, &u2[0], 1);
                // u3 = u*u*u
                Vmath::Vmul(npoints, &physfieldu[0], 1, &u2[0], 1, &u3[0], 1);

                // Ru = au
                Vmath::Smul(npoints, mA, &physfieldu[0], 1, &Ru[0], 1);
                // Ru = (-1-a)u*u + au
                Vmath::Svtvp(npoints, (-1.0-mA), &u2[0], 1, &Ru[0], 1, &Ru[0], 1);
                // Ru = u*u*u - (1+a)u*u + au
                Vmath::Vadd(npoints, &u3[0], 1, &Ru[0], 1, &Ru[0], 1);
                // Ru = k(u*u*u - (1+a)u*u + au)
                Vmath::Smul(npoints, mK, &Ru[0], 1, &Ru[0], 1);
                // Ru = k(u*u*u - (1+a)u*u + au) + uv
                Vmath::Vvtvp(npoints, &physfieldu[0], 1, &physfieldv[0], 1, &Ru[0], 1, &Ru[0], 1);
                // Ru = -k(u*u*u - (1+a)u*u + au) - uv
                Vmath::Neg(npoints, &Ru[0], 1);

                m_fields[0]->FwdTrans(Ru,outarray[0]);
                m_fields[0]->SetPhysState(false);

                // For v:
                // Rv = mu2 + u
                Vmath::Sadd(npoints, mMu2, &physfieldu[0], 1, &Rv[0], 1);
                // Rv = v/(mu2 + u)
                Vmath::Vdiv(npoints, &physfieldv[0], 1, &Rv[0], 1, &Rv[0], 1);
                // Rv = mu1*v/(mu2 + u)
                Vmath::Smul(npoints, mMu1, &Rv[0], 1, &Rv[0], 1);
                // Rv = Eps + mu1*v/(mu2+u)
                Vmath::Sadd(npoints, mEps, &Rv[0], 1, &Rv[0], 1);

                // Ru = (-a-1) + u
                Vmath::Sadd(npoints, (-mA-1), &physfieldu[0], 1, &Ru[0], 1);
                // Ru = k(u-a-1)
                Vmath::Smul(npoints, mK, &Ru[0], 1, &Ru[0], 1);
                // Ru = ku(u-a-1) + v
                Vmath::Vvtvp(npoints, &physfieldu[0], 1, &Ru[0], 1, &physfieldv[0], 1, &Ru[0], 1);
                // Ru = -ku(u-a-1)-v
                Vmath::Neg(npoints, &Ru[0], 1);

                Vmath::Vmul(npoints, &Ru[0], 1, &Rv[0], 1, &Rv[0], 1);

                m_fields[1]->FwdTrans(Rv,outarray[1]);
                m_fields[1]->SetPhysState(false);
    }

    void ADR2DManifold::ODEeReactionBarkley(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time)

    {
        int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> physfieldu(npoints);
        Array<OneD, NekDouble> physfieldv(npoints);

        Array<OneD, NekDouble> Ru(npoints,0.0);
        Array<OneD, NekDouble> Rv(npoints, 0.0);
        Array<OneD, NekDouble> u2(npoints,0.0);
        Array<OneD, NekDouble> u3(npoints,0.0);

        m_fields[0]->BwdTrans(inarray[0],physfieldu);
        m_fields[0]->SetPhysState(true);

        m_fields[1]->BwdTrans(inarray[1],physfieldv);
        m_fields[1]->SetPhysState(true);

        // u2 = u*u
                Vmath::Vmul(npoints, &physfieldu[0], 1, &physfieldu[0], 1, &u2[0], 1);
                // u3 = u*u*u
                Vmath::Vmul(npoints, &physfieldu[0], 1, &u2[0], 1, &u3[0], 1);

                // Ru = v + b - a
                Vmath::Sadd(npoints, mB - mA, &physfieldv[0], 1, &Ru[0], 1);
                // Ru = u(v + b - a) - v
                Vmath::Vvtvm(npoints, &Ru[0], 1, &physfieldu[0], 1, &physfieldv[0], 1, &Ru[0], 1);
                // Ru = u(v + b - a) - v - b
                Vmath::Sadd(npoints, -mB, &Ru[0], 1, &Ru[0], 1);
                // Ru = 1/a(u(v + b - a) - v - b)
                Vmath::Smul(npoints, 1/mA, &Ru[0], 1, &Ru[0], 1);
                // Ru = u*(1/a)*(u(v + b - a) - v - b) - u*u*u
                Vmath::Vvtvm(npoints, &physfieldu[0], 1, &Ru[0], 1, &u3[0], 1, &Ru[0], 1);
                // Ru = k[u*(1/a)*(u(v + b - a) - v - b) - u*u*u]
                Vmath::Smul(npoints, mK, &Ru[0], 1, &Ru[0], 1);

                m_fields[0]->FwdTrans(Ru,outarray[0]);
                m_fields[0]->SetPhysState(false);

                // Rv = u - v
                Vmath::Vsub(npoints, &physfieldu[0], 1, &physfieldv[0], 1, &Rv[0], 1);

                m_fields[1]->FwdTrans(Rv,outarray[1]);
                m_fields[1]->SetPhysState(false);
    }

    void ADR2DManifold::ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        NekDouble kappa;

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        SetBoundaryConditions(inarray,time);
        for (int i = 0; i < nvariables; ++i)
        {
            kappa = 1.0/(lambda*m_epsilon[i]);

            // Multiply rhs[i] with -1.0/gamma/timestep
            Vmath::Smul(ncoeffs, -1.0*kappa, &inarray[i][0], 1, &outarray[i][0], 1);

            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];

            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());

            SolveHelmholtz(i,kappa,m_UseDirDeriv);

            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetCoeffs();
        }
    }

    void ADR2DManifold::ODEhelmSolveFHNtest1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        NekDouble kappa = 1.0/(lambda*m_epsilon[0]);

        for (int i = 0; i < nvariables; ++i)
        {
            // For v: ==============================
            // Multiply rhs[i] with -1.0/gamma/timestep
            Vmath::Smul(ncoeffs, -1.0*kappa, &inarray[i][0], 1, &outarray[i][0], 1);

            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];

            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());

            SolveHelmholtz(i,kappa,m_UseDirDeriv);

            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetCoeffs();
        }
    }

    void ADR2DManifold::ODEhelmSolveFHNmono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda)
    {
        int ncoeffs    = inarray[0].num_elements();

        NekDouble kappa = 1.0/lambda;

        SetBoundaryConditions(inarray,time);

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        // Multiply rhs[i] with -1.0/gamma/timestep
        Vmath::Smul(ncoeffs, -1.0*kappa, &inarray[0][0], 1, &outarray[0][0], 1);

        // Update coeffs to m_fields
        m_fields[0]->UpdateCoeffs() = outarray[0];

        // Backward Transformation to nodal coefficients
        m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());

        SolveHelmholtz(0,kappa,m_UseDirDeriv);

        // The solution is Y[i]
        outarray[0] = m_fields[0]->GetCoeffs();

        // For q: No helmholtz solver is needed=============================
        Vmath::Vcopy(ncoeffs, &inarray[1][0], 1, &outarray[1][0], 1);
    }

    void ADR2DManifold:: SolveHelmholtz(const int indx, const NekDouble kappa,const int m_UseDirDeriv)
    {
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = kappa;
        if(m_UseDirDeriv)
        {
            //  m_fields[indx]->HelmSolve(m_fields[indx]->GetPhys(),m_fields[indx]->UpdateCoeffs(),kappa);

            m_fields[indx]->HelmSolve(m_fields[indx]->GetPhys(),m_fields[indx]->UpdateCoeffs(),NullFlagList,factors,StdRegions::NullVarCoeffMap,m_dirForcing[indx]);
        }

        else
        {
            m_fields[indx]->HelmSolve(m_fields[indx]->GetPhys(),m_fields[indx]->UpdateCoeffs(),NullFlagList,factors);
        }

        m_fields[indx]->SetPhysState(false);
    }

    // For Continuous Galerkin projections with time-dependent dirichlet boundary conditions,
    // the time integration can be done as follows:
    // The ODE resulting from the PDE can be formulated as:
    //
    // M du/dt = F(u)  or du/dt = M^(-1) F(u)
    //
    // Now suppose that M does not depend of time, the ODE can than be written as:
    //
    // d(Mu)/dt = F(u)
    //
    // Introducing the variable u* = Mu, this yields
    //
    // du*/dt = F( M^(-1) u*  ) = F*(u*)
    //
    // So rather than solving the initial ODE, it is advised to solve this new ODE for u*
    // as this allows for an easier treatment of the dirichlet boundary conditions.
    // However, note that at the end of every time step, the actual solution u can
    // be calculated as:
    //
    // u = M^(-1) u*;
    //
    // This can be viewed as projecting the solution u* onto the known boundary conditions.
    // Note that this step is also done inside the ODE rhs function F*.
    //
    // In order for all of this to work appropriately, make sure that the operator M^(-1)
    // does include the enforcment of the dirichlet boundary conditionst

    void ADR2DManifold::GeneralTimeIntegration(int nsteps,
            LibUtilities::TimeIntegrationMethod IntMethod,
            LibUtilities::TimeIntegrationSchemeOperators ode)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->SetPhysState(false);
            fields[i]  = m_fields[i]->UpdateCoeffs();
        }

        if((m_projectionType==MultiRegions::eGalerkin)||
           (m_projectionType==MultiRegions::eMixed_CG_Discontinuous))
        {
            // calculate the variable u* = Mu
            // we are going to TimeIntegrate this new variable u*
            MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
            for(int i = 0; i < nvariables; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(ncoeffs);
                m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);

                m_fields[i]->SetPhysState(false);
            }
        }

        // Declare an array of TimeIntegrationSchemes
        // For multi-stage methods, this array will have just one entry containing
        // the actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this is an initialization scheme)
        //  - the second scheme will used for the first timestep (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps (this will be the actual scheme)
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(IntMethod)
        {
        case LibUtilities::eIMEXdirk_3_4_3:
        case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
        {
            numMultiSteps = 1;

            IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

            LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
            IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

            u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
        }
        break;
        case LibUtilities::eAdamsBashforthOrder2:
        {
            numMultiSteps = 2;

            IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

            // Used in the first time step to initalize the scheme
            LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);

            // Used for all other time steps
            LibUtilities::TimeIntegrationSchemeKey IntKey1(IntMethod);
            IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
            IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

            // Initialise the scheme for the actual time integration scheme
            u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,ode);
        }
        break;
        default:
        {
            ASSERTL0(false,"populate switch statement for integration scheme");
        }
        }

        std::string outname = m_sessionName + ".his";
        std::ofstream hisFile (outname.c_str());

        for(n = 0; n < nsteps; ++n)
        {
            if (m_time < m_duration)
            {
                SetUSERDEFINEDInitialConditions(m_initialwavetype, false);
            }

            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
            }

            m_time += m_timestep;

            if((m_projectionType==MultiRegions::eGalerkin)||
               (m_projectionType==MultiRegions::eMixed_CG_Discontinuous))
            {
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                // SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i]);
                    fields[i] = tmp[i];
                }
            }

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%m_infosteps))
            {
                cout << "Steps: " << n+1 << "\t Time: " << m_time << endl;
            }

            if(n&&(!((n+1)%m_checksteps)))
            {
                for(i = 0; i < nvariables; ++i)
                {
                    (m_fields[i]->UpdateCoeffs()) = fields[i];
                }
                Checkpoint_Output(nchk++);
            }
        }


        for(i = 0; i < nvariables; ++i)
        {
            (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }

    //-------------------------------------------------------------
    // Calculate weak DG advection in the form
    //  <\phi, \hat{F}\cdot n> - (\grad \phi \cdot F)
    // -------------------------------------------------------------
    void ADR2DManifold::WeakDGAdvectionDir(const Array<OneD, Array<OneD, NekDouble> >& InField,
            const Array<OneD, Array<OneD, NekDouble> >& dirForcing,
            Array<OneD, Array<OneD, NekDouble> >& OutField,
            bool NumericalFluxIncludesNormal, bool InFieldIsInPhysSpace, int nvariables)
    {
        int i;
        int nPointsTot      = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();

        if (!nvariables)
        {
            nvariables      = m_fields.num_elements();
        }

        Array<OneD, Array<OneD, NekDouble> > physfield (nvariables);
        // Get the variables in physical space already in physical space
        if(InFieldIsInPhysSpace == true)
        {
            for(i = 0; i < nvariables; ++i)
            {
                physfield[i] = InField[i];
            }
        }

        // otherwise do a backward transformation
        else
        {
            for(i = 0; i < nvariables; ++i)
            {
                // Could make this point to m_fields[i]->UpdatePhys();
                physfield[i] = Array<OneD, NekDouble>(nPointsTot);
                m_fields[i]->BwdTrans(InField[i],physfield[i]);
            }
        }

        // Get the directional derivative along two tangential vectors t_eta and t_xi
        // Du = v_eta*D(t_eta*u) + v_xi*D(t_xi*u)
        for(i = 0; i < nvariables; ++i)
        {
            WeakDirectionalDerivative(physfield[i], dirForcing, OutField[i]);
        }

        // Get the numerical flux and add to the modal coeffs
        // if the NumericalFlux function already includes the
        // normal in the output
        if (NumericalFluxIncludesNormal == true)
        {
            Array<OneD, Array<OneD, NekDouble> > numflux   (nvariables);
            for(i = 0; i < nvariables; ++i)
            {
                numflux[i]   = Array<OneD, NekDouble>(nTracePointsTot);
            }

            // Evaluate numerical flux in physical space which may in
            // general couple all component of vectors

            NumericalFluxdir(physfield, dirForcing, numflux);

            // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,OutField[i],1);
                m_fields[i]->AddTraceIntegral(numflux[i],OutField[i]);
                m_fields[i]->SetPhysState(false);
            }
        }
        // if the NumericalFlux function does not include the
        // normal in the output
        else
        {
            Array<OneD, Array<OneD, NekDouble> > numfluxX   (nvariables);
            Array<OneD, Array<OneD, NekDouble> > numfluxY   (nvariables);

            for(i = 0; i < nvariables; ++i)
            {
                numfluxX[i]   = Array<OneD, NekDouble>(nTracePointsTot);
                numfluxY[i]   = Array<OneD, NekDouble>(nTracePointsTot);
            }

            // Evaluate numerical flux in physical space which may in
            // general couple all component of vectors
            NumericalFlux(physfield, numfluxX, numfluxY);

            // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,OutField[i],1);
                m_fields[i]->AddTraceIntegral(numfluxX[i],numfluxY[i],OutField[i]);
                m_fields[i]->SetPhysState(false);
            }
        }
    }

    // Compute Directional derivative along a tangential vector (dirForcing)
    // input: u in physfield
    // output: v_eta ( D_{t_eta} + M_{t_eta} ) +  v_xi ( D_{t_xi} + M_{t_xi} )
    void ADR2DManifold::WeakDirectionalDerivative(const Array<OneD, NekDouble> &physfield,
            const Array<OneD, Array<OneD, NekDouble> > &dirForcing,
            Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int i;
        int nCoeffs = m_fields[0]->GetNcoeffs();
        int nq = GetNpoints();

        Array<OneD, NekDouble> iprod(nCoeffs);
        Array<OneD, Array<OneD, NekDouble> > fluxvector(m_spacedim);

        Vmath::Zero(nCoeffs,&outarray[0],1);
        // Multiply u by tangential vectors which give a vector
        for (i = 0; i < m_spacedim; ++i)
        {
            // Generate t_eta*u and t_xi*u
            fluxvector[i] = Array<OneD, NekDouble>(nq,0.0);
            Vmath::Vmul(nq, &physfield[0], 1, &dirForcing[0][i*nq], 1, &fluxvector[i][0],1);

            // Differentiate t_xi*u
            m_fields[0]->IProductWRTDerivBase(i,fluxvector[i],iprod);

            // Add each Eulerian coordiate components
            Vmath::Vadd(nCoeffs,&iprod[0],1,&outarray[0],1,&outarray[0],1);
        }
    }


    void ADR2DManifold::NumericalFluxdir(Array<OneD, Array<OneD, NekDouble> > &physfield,
            const Array<OneD, Array<OneD, NekDouble> > &dirForcing,
            Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;
        int nq         = GetTotPoints();
        int nTraceNumPoints = GetTraceNpoints();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
        Array<OneD, NekDouble > temp(nq,0.0);

        switch(m_equationType)
        {
        case eUnsteadyAdvection:
        {
            for(i = 0; i < numflux.num_elements(); ++i)
            {
                m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);

                // evaulate upwinded m_fields[i]
                                              m_fields[i]->GetTrace()->Upwind(m_Vn,Fwd,Bwd,numflux[i]);

                                              // calculate m_fields[i]*Vn
                                              Vmath::Vmul(nTraceNumPoints,&numflux[i][0],1,&m_Vn[0],1,&numflux[i][0],1);
            }
        }
        break;

        default:
            ASSERTL0(false,"unknown equationType");
        }
    }



    //----------------------------------------------------
    // inarray in coeff fields
    //----------------------------------------------------
    void ADR2DManifold::SetBoundaryConditions(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            NekDouble time)
    {
        int i,n,cnt=0;
        int nvariables = m_fields.num_elements();

        for (n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Read Boundary Condition from a function in this cpp file
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eMG)
            {
                MGBoundaryCondtions(n,cnt,time);
            }

            else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWALL)
            {
                WallBoundary2D(n,cnt,inarray);
            }

            else
            {
                // Otherwise read it from xml file
                for (i=0; i < nvariables; ++i)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }

            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }

    void ADR2DManifold::MGBoundaryCondtions(const int bcRegion, const int cnt, const NekDouble time)
    {

        int i,j,indx;
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = m_fields.num_elements();

        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);

        // Extract coordinates of grid points on trace
        Array<OneD,NekDouble> x0(nTraceNumPoints,0.0);
        Array<OneD,NekDouble> x1(nTraceNumPoints,0.0);
        Array<OneD,NekDouble> x2(nTraceNumPoints,0.0);

        m_fields[0]->GetTrace()->GetCoords(x0,x1,x2);

        int e, id1, id2, npts;
        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

            for(i = 0 ; i < nvariables; i++)
            {
                Fwd[i] = Array<OneD, NekDouble>(npts,0.0);

                for (j =0; j < npts; ++j)
                {
                    indx = id2+j;
                    Fwd[i][j] = Morphogenesis(i, x0[indx], x1[indx], x2[indx], time);
                }

                Vmath::Vcopy(npts,&Fwd[i][0], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            }
        }

        for(i =0; i < nvariables; i++)
        {
            m_fields[i]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[i]->GetBndCondExpansions()[bcRegion]->GetPhys(),
                    m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
        }
    }

    void ADR2DManifold::WallBoundary2D(const int bcRegion, const int cnt,
            const Array<OneD, const Array<OneD, NekDouble> >&inarray)
    {

        int i;
        int nq         = GetTotPoints();
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = inarray.num_elements();

        Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            physarray[i] = Array<OneD, NekDouble>(nq);
            m_fields[i]->BwdTrans(inarray[i],physarray[i]);
        }

        // get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }

        // Adjust the physical values of the trace to take
        // user defined boundaries into account
        int e, id1, id2, npts;

        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Smul(npts, -1.0, &Fwd[i][id2], 1,
                        &(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1], 1);
            }
        }

        for(i =0; i < nvariables; i++)
        {
            m_fields[i]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[i]->GetBndCondExpansions()[bcRegion]->GetPhys(),
                    m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
        }
    }


    // Evaulate flux = m_fields*ivel for i th component of Vu
    // alt
    //          flux = 0.5(m_fields*m_fields) for the
    void ADR2DManifold::GetFluxVector(const int i,
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &flux)
    {

        switch(m_equationType)
        {
        case eUnsteadyAdvection:
        case eUnsteadyDiffusion:
        case eUnsteadyDiffusionReaction:
        case eNonlinearMorphogenensis:
        case eFHNtesttype1:
        case eFHNMONO:
        case eAlievPanfilov:
        {
            ASSERTL1(flux.num_elements() == m_velocity.num_elements(),
                    "Dimension of flux array and velocity array do not match");

            for(int j = 0; j < flux.num_elements(); ++j)
            {
                Vmath::Vmul(GetNpoints(),physfield[i],1,
                        m_velocity[j],1,flux[j],1);
            }
        }
        break;
        default:
            ASSERTL0(false,"unknown equationType");
        }
    }

    // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
    void ADR2DManifold::GetFluxVector(const int i, const int j,
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        switch(m_equationType)
        {
    #ifdef OLDENUM
        case eAdvection:
    #else
        case eUnsteadyAdvection:
        case eUnsteadyDiffusion:
        case eUnsteadyDiffusionReaction:
        case eNonlinearMorphogenensis:
        case eFHNtesttype1:
        case eFHNMONO:
        case eAlievPanfilov:
    #endif
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),
                "Dimension of flux array and velocity array do not match");

        for(int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(),flux[k],1);
        }
        Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }
    break;
        default:
            ASSERTL0(false,"unknown equationType");
        }
    }


    void ADR2DManifold::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        switch(m_equationType)
        {
        case eUnsteadyAdvection:
        {
            for(i = 0; i < numflux.num_elements(); ++i)
            {
                m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);

                //evaulate upwinded m_fields[i]
                m_fields[i]->GetTrace()->Upwind(m_Vn,Fwd,Bwd,numflux[i]);

                // calculate m_fields[i]*Vn
                Vmath::Vmul(nTraceNumPoints,numflux[i],1,m_Vn,1,numflux[i],1);
            }
        }
        break;
        default:
            ASSERTL0(false,"unknown equationType");
        }
    }

    void ADR2DManifold::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
            Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int nTraceNumPoints = GetTraceNpoints();
        
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > tmp(nTraceNumPoints,0.0);

        Array<OneD, Array<OneD, NekDouble > > traceVelocity(2);

        traceVelocity[0] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);
        traceVelocity[1] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);

        switch(m_equationType)
        {
        case eUnsteadyAdvection:
        {

            // Get Edge Velocity - Could be stored if time independent
            m_fields[0]->ExtractTracePhys(m_velocity[0], traceVelocity[0]);
            m_fields[0]->ExtractTracePhys(m_velocity[1], traceVelocity[1]);

            m_fields[0]->GetFwdBwdTracePhys(physfield[0],Fwd,Bwd);

            m_fields[0]->GetTrace()->Upwind(traceVelocity,Fwd,Bwd,tmp);

            Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[0],1,numfluxX[0],1);
            Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[1],1,numfluxY[0],1);
        }
        break;
        default:
            ASSERTL0(false,"unknown equationType");
        }
    }

    // Compute the fluxes of q from the scalar functin u.
    // Input:   ufield (1 by nTraceNumPoints) - Should be in physical field
    // Output:  ufluxFwd  (2 by nTraceNumPoints) - Flux values for forward edges
    //          ufluxBwd  (2 by nTraceNumPoints) - Flux values for backward edges
    void ADR2DManifold::NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
    {
        int i,j;
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = uflux.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
        Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);

        // Get the sign of (v \cdot n), v = an arbitrary vector

        // Vn = V \cdot n, where n is tracenormal for eForward edges. Set V = (1,0)
        // Vmath::Vcopy(nTraceNumPoints,m_traceNormals_tbasis[0],1,Vn,1);

        //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
        for (j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables ; ++i)
            {
                //  Compute Forward and Backward value of ufield of i direction
                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                m_fields[i]->GetTrace()->Upwind(m_traceNormals_tbasis[j],Fwd,Bwd,fluxtemp);

                // Imposing weak boundary condition with flux
                // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd
                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforScalar(i,ufield[i],fluxtemp);
                }

                // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                Vmath::Vmul(nTraceNumPoints,m_traceNormals_tbasis[j],1,fluxtemp,1,uflux[j][i],1);
            }
        }
    }

    // Compute the fluxes of q and u vector fields for discontinuous diffusion term
    // Input:   qfield : 2 by # of total trace points
    // Output:  qflux  : 2 by # of total trace points
    void ADR2DManifold::NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, NekDouble> >  &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = qfield.num_elements();

        NekDouble C11 = 1.0;
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
        Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

        Array<OneD, NekDouble > uterm(nTraceNumPoints);

        // Get the sign of (v \cdot n), v = an arbitrary vector
        // Vn = V \cdot n, where n is tracenormal for eForward edges
        // Vmath::Vcopy(nTraceNumPoints,m_traceNormals[0],1,Vn,1);

        // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
        for(int i = 0; i < nvariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
            for(int j = 0; j < nqvar; ++j)
            {
                //  Compute Forward and Backward value of ufield of jth direction
                m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd =q+
                m_fields[i]->GetTrace()->Upwind(m_traceNormals_tbasis[j],qBwd,qFwd,qfluxtemp);
                Vmath::Vmul(nTraceNumPoints,m_traceNormals_tbasis[j],1,qfluxtemp,1,qfluxtemp,1);

                // Generate Stability term = - C11 ( u- - u+ )
                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
                Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);
                Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

                //  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
                Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);

                // Imposing weak boundary condition with flux
                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
                }

                // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                // n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
                // n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
            }
        }
    }


    // Diffusion: Imposing weak boundary condition for u with flux
    //  uflux = g_D  on Dirichlet boundary condition
    //  uflux = u_Fwd  on Neumann boundary condition
    void ADR2DManifold::WeakPenaltyforScalar(const int var,
            const Array<OneD, const NekDouble> &physfield,
            Array<OneD, NekDouble> &penaltyflux,
            NekDouble time)
    {
        unsigned int i, e, npoints, id1, id2;
        int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
        int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
        int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
        int nTraceNumPoints = GetTraceNpoints();

        Array<OneD, NekDouble > uplus(nTraceNumPoints);

        m_fields[var]->ExtractTracePhys(physfield,uplus);
        for(i = 0; i < nbnd; ++i)
        {
            // Evaluate boundary values g_D or g_N from input files
            LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("InitialConditions", i);
            npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();

            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            ifunc->Evaluate(x0,x1,x2,time,BDphysics);

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
                id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));

                // For Dirichlet boundary condition: uflux = g_D
                if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
                }

                // For Neumann boundary condition: uflux = u+
                else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                }
            }
        }
    }

    // Diffusion: Imposing weak boundary condition for q with flux
    //  uflux = g_D  on Dirichlet boundary condition
    //  uflux = u_Fwd  on Neumann boundary condition
    void ADR2DManifold::WeakPenaltyforVector(const int var,
            const int dir,
            const Array<OneD, const NekDouble> &physfield,
            Array<OneD, NekDouble> &penaltyflux,
            NekDouble C11,
            NekDouble time)
    {
        unsigned int i, e, npoints, id1, id2;
        int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
        int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
        int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
        int nTraceNumPoints = GetTraceNpoints();
        Array<OneD, NekDouble > uterm(nTraceNumPoints);
        Array<OneD, NekDouble > qtemp(nTraceNumPoints);

        m_fields[var]->ExtractTracePhys(physfield,qtemp);

        for(i = 0; i < nbnd; ++i)
        {
            // Evaluate boundary values g_D or g_N from input files
            LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("InitialConditions",i);
            npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();

            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            ifunc->Evaluate(x0,x1,x2,time,BDphysics);

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
                id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));

                // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ - g_D) (nx, ny)
                if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals_tbasis[dir][id2],1,&qtemp[id2],1,&penaltyflux[id2],1);

                    // Vmath::Vsub(Nfps,&Fwd[id2],1,&BDphysics[id1],1,&uterm[id2],1);
                    //    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&uterm[id2],1,&uterm[id2],1);
                    // Vmath::Svtvp(Nfps,-1.0*C11,&uterm[id2],1,&qFwd[id2],1,&penaltyflux[id2],1);
                }

                // For Neumann boundary condition: qflux = g_N
                else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals_tbasis[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
                }
            }
        }
    }


    void ADR2DManifold::SetUpTangentialVectors()
    {
        int i, k;
        int nvariables = m_fields.num_elements();
        int nq = m_fields[0]->GetNpoints();

        // Trace Normal vectors on the plane of tangential basis
        m_dirForcing = Array<OneD, Array<OneD, NekDouble> >(2);
        m_gradtan =  Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
        m_tanbasis  =  Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);

        for (i = 0; i < 2; ++i)
        {
            m_dirForcing[i] = Array<OneD, NekDouble>(m_spacedim*nq);
            m_tanbasis[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            for (k = 0; k < m_spacedim; ++k)
            {
                m_tanbasis[i][k] = Array<OneD, NekDouble>(nq, 0.0);
            }

            m_gradtan[i] = Array<OneD, Array<OneD, NekDouble> >(nvariables);
            for (k = 0; k < nvariables; ++k)
            {
                m_gradtan[i][k] = Array<OneD, NekDouble>(nq, 0.0);
            }
        }

        m_fields[0]->GetTangents(m_tanbasis);

        // When there is an anisotropy layers, the second tangential vector is
        // killed.
        if(m_Anisotropy>0)
        {
            for(k = 0; k < m_spacedim; ++k)
            {
                m_tanbasis[1][k] = Array<OneD, NekDouble>(nq,0.0);
            }
        }

        int nelem ,npts_e, offset0, offset1;
        // Using this tangential basis as dirForcing

        nelem = GetExpSize();
        for (int i=0; i<m_spacedim; ++i)
        {
            for (int k=0; k<2; ++k)
            {
                offset1 = 0;
                for (int j=0; j<nelem; ++j)
                {
                    npts_e = GetTotPoints(j);
                    offset0 =  GetPhys_Offset(j);

                    Vmath::Vcopy(npts_e,&m_tanbasis[k][i][offset0],1,&m_dirForcing[k][offset1+i*npts_e],1);

                    offset1 += npts_e*m_spacedim;
                }
            }
        }

        Vmath::Neg(offset1, &m_dirForcing[1][0], 1);

        Array<OneD, NekDouble> temp0(nq), temp1(nq), temp2(nq);
        // Generate \nabla \cdot tangentvector for weighted mass matrix
        for(int i = 0; i < nvariables; ++i)
        {
            for (int k = 0; k < 2; ++k)
            {
                m_fields[0]->PhysDeriv(m_tanbasis[k][0], temp0, temp1, temp2);
                Vmath::Vadd(nq, temp0, 1, m_gradtan[k][i], 1, m_gradtan[k][i], 1);

                m_fields[0]->PhysDeriv(m_tanbasis[k][1], temp0, temp1, temp2);
                Vmath::Vadd(nq, temp1, 1, m_gradtan[k][i], 1, m_gradtan[k][i], 1);

                m_fields[0]->PhysDeriv(m_tanbasis[k][2], temp0, temp1, temp2);
                Vmath::Vadd(nq, temp2, 1, m_gradtan[k][i], 1, m_gradtan[k][i], 1);
            }
        }
    }

    void ADR2DManifold::SetUSERDEFINEDInitialConditions(const int initialwavetype, bool SetRestingState, NekDouble initialtime)
    {
        int nq = m_fields[0]->GetNpoints();
        int nvariables = m_fields.num_elements();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        Array<OneD, NekDouble> rstates(2);
        if (SetRestingState)
        {
            Getrestingstate(m_epsilon[0], m_beta, rstates);
        }

        // get the coordinates (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        switch(initialwavetype)
        {
        // -2: Advection problem on sphere
        case(-2):
                    {
            for(int i = 0 ; i <nvariables ; i++)
            {
                for(int j = 0; j < nq; j++)
                {
                    (m_fields[i]->UpdatePhys())[j] = AdvectionSphere(x0[j], x1[j], x2[j], initialtime);
                }
    
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;
    
        // -1:Morphogenesis problem for diffusion-reaction test
        case(-1):
                    {
            for(int i = 0 ; i < nvariables ; i++)
            {
                for(int j = 0; j < nq; j++)
                {
                    (m_fields[i]->UpdatePhys())[j] = Morphogenesis(i, x0[j], x1[j], x2[j], initialtime);
                }
    
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;
    
        // 0: Planar propagation from -x direction
        case(0):
                    {
            // Set the left side as the initial excitation
            NekDouble xmin = Vmath::Vmin(nq,x0,1);
            for(int j = 0; j < nq; j++)
            {
                // Set U
                if( x0[j] <= (xmin+m_duration) )
                {
                    (m_fields[0]->UpdatePhys())[j] = 2.0;
                }
                else if (SetRestingState)
                {
                    (m_fields[0]->UpdatePhys())[j] = rstates[0];
                }
                // Set V
                if (SetRestingState)
                {
                    (m_fields[1]->UpdatePhys())[j] = rstates[1];
                }
            }
    
            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;

        // 1: Planar propagation from -y direction
        case(1):
                    {
            // Set the left side as the initial excitation
            NekDouble ymin = Vmath::Vmin(nq,x1,1);
            for(int j = 0; j < nq; j++)
            {
                if( x2[j] > (ymin+m_duration) )
                {
                    (m_fields[0]->UpdatePhys())[j] = 2.0;
                }
                else if (SetRestingState)
                {
                    (m_fields[0]->UpdatePhys())[j] = rstates[0];
                }
                if (SetRestingState)
                {
                    (m_fields[1]->UpdatePhys())[j] = rstates[1];
                }
            }
    
            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;

        // 2: For Planar propagation from -x direction in xy plane
        case(2):
                    {
            // Set the left side as the initial excitation
            NekDouble xmin = Vmath::Vmin(nq,x0,1);
            for(int j = 0; j < nq; j++)
            {
                if( ( x0[j] <= (xmin+m_duration) ) && ( fabs(x2[j]) <= 0.01 ) )
                {
                    (m_fields[0]->UpdatePhys())[j] = 2.0;
                }
                else if (SetRestingState)
                {
                    (m_fields[0]->UpdatePhys())[j] = rstates[0];
                }
                if (SetRestingState)
                {
                    (m_fields[1]->UpdatePhys())[j] = rstates[1];
                }
            }
    
            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;
    
        case(3):
                    {
            NekDouble rad;
            for(int j = 0; j < nq; j++)
            {
                rad = sqrt( (x0[j]-m_x0c)*(x0[j]-m_x0c) + (x1[j]-m_x1c)*(x1[j]-m_x1c) + (x2[j]-m_x2c)*(x2[j]-m_x2c) );

                if (SetRestingState)
                {
                    (m_fields[1]->UpdatePhys())[j] = rstates[1];
                }
                if( rad <= mDiam )
                {
                    (m_fields[0]->UpdatePhys())[j] = 1.0;
                }
                else if (SetRestingState)
                {
                    (m_fields[0]->UpdatePhys())[j] = rstates[0];
                }
            }
    
            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
                    }
        break;

        // 10: For Planar propagation from +z direction
        case(10):
        {
            // Set the left side as the initial excitation
            NekDouble zmax = Vmath::Vmax(nq,x2,1);
            for(int j = 0; j < nq; j++)
            {
                if (SetRestingState)
                {
                    (m_fields[1]->UpdatePhys())[j] = rstates[1];
                }
                if( x2[j] > (zmax-m_duration) )
                {
                    (m_fields[0]->UpdatePhys())[j] = 2.0;
                }
                else if (SetRestingState)
                {
                    (m_fields[0]->UpdatePhys())[j] = rstates[0];
                }
            }
    
            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
        }
        break;

        default:
        {
            if (m_time == 0.0)
            {
                SetInitialConditions(initialtime);
            }
        }
        }

        // dump initial conditions to file
        if (m_time == 0.0)
        {
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
        }
    }

    void ADR2DManifold::EvaluateUSERDEFINEDExactSolution(unsigned int field, Array<OneD, NekDouble> &outfield,
            const NekDouble time, const int initialwavetype)
    {
        int nq = m_fields[field]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates of the quad points
        m_fields[field]->GetCoords(x0,x1,x2);

        switch(initialwavetype)
        {
        case(-2):
                    {
            for(int j = 0; j < nq; ++j)
            {
                outfield[j] = AdvectionSphere(x0[j], x1[j], x2[j], time);
            }
                    }
        break;

        case(-1):
                    {
            for(int j = 0; j < nq; ++j)
            {
                outfield[j] = Morphogenesis(field, x0[j], x1[j], x2[j], time);
            }
                    }
        break;

        case(0):
        case(1):
        case(2):
        case(3):
        case(10):
        {
            for(int j = 0; j < nq; ++j)
            {
                Vmath::Zero(nq,outfield,1);
            }
        }
        break;

        default:
        {
            EvaluateFunction(m_session->GetVariable(field), outfield, "ExactSolution", time);
        }
        }
    }


    NekDouble ADR2DManifold::AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
            const NekDouble x2j, const NekDouble time)
    {

        NekDouble dist, radius, cosdiff, sin_theta, cos_theta, sin_varphi, cos_varphi;
        NekDouble pi = 3.14159265358979323846;
        NekDouble newvarphi_c, newtheta_c, returnval;
        NekDouble m_theta_c, m_varphi_c, m_radius_limit, m_c0;

        Array<OneD, NekDouble> Exactsolution;

        // Sets of parameters
        m_theta_c = 0.0;
        m_varphi_c = 3.0*pi/2.0;
        m_radius_limit = 7.0*pi/64.0;
        m_c0 = 2.0;

        newvarphi_c = m_varphi_c ;
        newtheta_c = m_theta_c ;

        radius = sqrt( x0j*x0j + x1j*x1j + x2j*x2j );

        sin_varphi = x1j/sqrt( x0j*x0j + x1j*x1j );
        cos_varphi = x0j/sqrt( x0j*x0j + x1j*x1j );

        sin_theta = x2j/radius;
        cos_theta = sqrt( x0j*x0j + x1j*x1j )/radius;

        cosdiff =  cos_varphi*cos(m_varphi_c) + sin_varphi*sin(m_varphi_c);
        dist = radius*acos( sin(m_theta_c)*sin_theta + cos(m_theta_c)*cos_theta*cosdiff );

        if(dist < m_radius_limit)
        {
            returnval = 0.5*( 1.0 + cos(pi*dist/m_radius_limit) ) + m_c0;
        }
        else
        {
            returnval = m_c0;
        }

        return returnval;
    }

    NekDouble ADR2DManifold::Morphogenesis(const int field, const NekDouble x0j, const NekDouble x1j,
            const NekDouble x2j, const NekDouble time)

    {
        NekDouble pi = 3.14159265358979323846;

        int i, m, n, ind;
        NekDouble a_n, d_n, gamma_n;
        NekDouble A_mn, C_mn, theta, phi,radius;

        std::complex<NekDouble> Spericharmonic, delta_n, varphi0, varphi1, temp;
        std::complex<NekDouble> B_mn, D_mn;

        // Set some parameter values
        int Maxn = 6;
        int Maxm = 2*Maxn-1;

        NekDouble A = 2.0;
        NekDouble B = 5.0;

        NekDouble m_mu = m_epsilon[0];
        NekDouble m_nu = m_epsilon[1];

        m_a = B-1.0;
        m_b = A*A;
        m_c = -1.0*B;
        m_d = -1.0*A*A;

        Array<OneD, Array<OneD, NekDouble> > Ainit(Maxn);
        Array<OneD, Array<OneD, NekDouble> > Binit(Maxn);

        for (i = 0; i < Maxn; ++i)
        {
            Ainit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
            Binit[i] = Array<OneD, NekDouble>(Maxm, 0.0);
        }

        Ainit[5][0] = -0.5839;
        Ainit[5][1] = -0.8436;
        Ainit[5][2] = -0.4764;
        Ainit[5][3] = 0.6475;
        Ainit[5][4] = 0.1886;
        Ainit[5][5] = 0.8709;
        Ainit[5][6] = -0.8338;
        Ainit[5][7] = 0.1795;
        Ainit[5][8] = -0.7873;
        Ainit[5][9] = 0.8842;
        Ainit[5][10] = 0.2943;

        Binit[5][0] = -0.6263;
        Binit[5][1] = 0.9803;
        Binit[5][2] = 0.7222;
        Binit[5][3] = 0.5945;
        Binit[5][4] = 0.6026;
        Binit[5][5] = -0.2076;
        Binit[5][6] = 0.4556;
        Binit[5][7] = 0.6024;
        Binit[5][8] = 0.9695;
        Binit[5][9] = -0.4936;
        Binit[5][10] = 0.1098;

        radius = sqrt(x0j*x0j + x1j*x1j + x2j*x2j) ;

        // theta is in [0, pi]
        theta = asin( x2j/radius ) + 0.5*pi;

        // phi is in [0, 2*pi]
        phi = atan2( x1j, x0j ) + pi;

        varphi0 = std::complex<NekDouble>(0.0,0.0);
        varphi1 = std::complex<NekDouble>(0.0,0.0);
        for (n = 0; n < Maxn; ++n)
        {
            // Set up parameters
            a_n = m_a - m_mu*( n*(n+1)/radius/radius );
            d_n = m_d - m_nu*( n*(n+1)/radius/radius );

            gamma_n = 0.5*( a_n + d_n );

#ifdef _MSC_VER
            temp.real( ( a_n + d_n )*( a_n + d_n ) - 4.0*( a_n*d_n - m_b*m_c ));
#else
            temp = std::complex<NekDouble>( (a_n + d_n )*( a_n + d_n ) - 4.0*( a_n*d_n - m_b*m_c ), temp.imag() ) ;
#endif

            delta_n = 0.5*sqrt( temp );

            for (m = -n; m <=n; ++m)
            {
                ind = m + n;
                A_mn = Ainit[n][ind];
                C_mn = Binit[n][ind];

                B_mn = ( (a_n - gamma_n)*Ainit[n][ind] + m_b*Binit[n][ind])/delta_n;
                D_mn = ( m_c*Ainit[n][ind] + (d_n - gamma_n)*Binit[n][ind])/delta_n;

                Spericharmonic = boost::math::spherical_harmonic(n, m, theta, phi);
                varphi0 += exp(gamma_n*time)*(A_mn*cosh(delta_n*time) + B_mn*sinh(delta_n*time))*Spericharmonic;
                varphi1 += exp(gamma_n*time)*(C_mn*cosh(delta_n*time) + D_mn*sinh(delta_n*time))*Spericharmonic;
            }
        }

        if (field==0)
        {
            return varphi0.real();
        }

        else if (field==1)
        {
            return varphi1.real();
        }

        else
        {
            ASSERTL0(false, "Invalid field value.");
            return 0.0;
        }
    }


    // Plot tangential vector map
    void ADR2DManifold::PlotTangentialVectorMap()
    {
        int i,j,k,nqtot = GetNpoints();

        for(k = 0; k < 2; ++k)
        {
            for(i = 0 ; i < 3; i++)
            {
                for(j = 0; j < nqtot; j++)
                {
                    (m_fields[i+3*k]->UpdatePhys())[j] = 0.0;
                }
            }
        }

        for(k = 0; k < 2; ++k)
        {
            for(i = 0 ; i < 3; i++)
            {
                for(j = 0; j < nqtot; j++)
                {
                    (m_fields[i+3*k]->UpdatePhys())[j] = m_dirForcing[k][i*nqtot+j];
                }
            }
        }

        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
        }

        // dump initial conditions to file
        std::string outname = m_sessionName + "_TanVecMap.fld";
        WriteFld(outname);
    }

    void ADR2DManifold::Getrestingstate(const NekDouble epsilon, const NekDouble beta,
            Array<OneD, NekDouble> rstates)
    {
        NekDouble unew=100.0, uinit=0.0, vinit, f, fd, fval, gval,Tol=0.00000001;

        // Get the initial constant u and v
        for(int j = 0; j<1000; j++)
        {
            // f = fvalue(uinit);
            f = uinit*uinit*uinit + 3.0*uinit + 6.0*m_beta;

            // fd = fderiv(uinit);
            fd = 3.0*uinit*uinit + 3.0;

            unew = uinit - f/fd;

            if(abs(unew-uinit) < Tol)
            {
                break;
            }
            uinit = unew;
        }
        vinit = uinit - (1.0/3.0)*uinit*uinit*uinit;

        fval = uinit - (1.0/3.0)*uinit*uinit*uinit - vinit;
        gval = uinit + m_beta - 0.5*vinit;

        if ( (fabs(fval)<Tol) && (fabs(gval)<Tol) )
        {
            rstates[0] = uinit;
            rstates[1] = vinit;
        }
    }


    NekDouble ADR2DManifold::USERDEFINEDError(int field, const int type, const int initialwavetype,
            Array<OneD, NekDouble> &exactsoln)
    {

        NekDouble error=0;

        if(m_fields[field]->GetPhysState() == false)
        {
            m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                    m_fields[field]->UpdatePhys());
        }

        exactsoln = Array<OneD, NekDouble>(m_fields[field]->GetNpoints());
        EvaluateUSERDEFINEDExactSolution(field,exactsoln,m_time,initialwavetype);

        switch(type)
        {
        case(0):
                    {
            error = m_fields[field]->Linf(exactsoln);
                    }
        break;

        case(1):
                    {
            error = m_fields[field]->H1(exactsoln);
                    }
        break;

        case(2):
                    {
            error = m_fields[field]->L2(exactsoln);
                    }
        break;

        default:
        {
            error = m_fields[field]->L2(exactsoln);
            break;
        }
        }

        return error;
    }

    void ADR2DManifold::AdditionalSessionSummary(std::ostream &out)
    {

        out << "\tinitialwavetype    : " << m_initialwavetype << endl;
        out << "\tinitialcenter = ( " << m_x0c << "," << m_x1c << "," << m_x2c << " )" << endl;
        out << "\tUseDirDeriv   : " << m_UseDirDeriv << endl;
        out << "\tTangentDir    : ";
        if (m_graph->CheckForGeomInfo("TangentDir"))
        {
            out << m_graph->GetGeomInfo("TangentDir") << endl;
            if (m_graph->GetGeomInfo("TangentDir") == "TangentCircular")
            {
                out << "\tTangentCentre : ";
                if (m_graph->CheckForGeomInfo("TangentCentreX")
                        && m_graph->CheckForGeomInfo("TangentCentreY"))
                {
                    out << "("  << m_graph->GetGeomInfo("TangentCentreX")
                                    << ", " << m_graph->GetGeomInfo("TangentCentreY")
                                    << ")"  << endl;
                }
                else
                {
                    out << "Not defined." << endl;
                }
            }
        }
        else
        {
            out << "Not defined." << endl;
        }

        out << "\tAnisotropy    : " << m_Anisotropy << endl;
        switch(m_equationType)
        {
        case eAlievPanfilov:
            out << "\tmK            : " << mK << endl;
            out << "\tmA            : " << mA << endl;
            out << "\tmEps          : " << mEps << endl;
            out << "\tmMu1          : " << mMu1 << endl;
            out << "\tmMu2          : " << mMu2 << endl;
            break;
        default:
            out << "\tm_mu          : " << m_epsilon[0] << endl;
            out << "\tm_nu          : " << m_epsilon[1] << endl;
            out << "\tm_beta        : " << m_beta << endl;
        }
    }

    void ADR2DManifold::Summary(std::ostream &out)
    {
        cout << "=======================================================================" << endl;
        cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
        EquationSystem::SessionSummary(out);
        AdditionalSessionSummary(out);
        switch(m_equationType)
        {
        case eUnsteadyAdvection:
            if(m_explicitAdvection)
            {
                out << "\tAdvection Advancement   : Explicit" <<endl;
            }
            else
            {
                out << "\tAdvection Advancement   : Implicit" <<endl;
            }
            out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;

            EquationSystem::TimeParamSummary(out);
            break;

        case eUnsteadyDiffusion:
            if(m_explicitDiffusion)
            {
                out << "\tDiffusion Advancement   : Explicit" <<endl;
            }
            else
            {
                out << "\tDiffusion Advancement   : Implicit" <<endl;
            }
            out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
            EquationSystem::TimeParamSummary(out);
            break;

        case eUnsteadyDiffusionReaction:
        case eNonlinearMorphogenensis:
        case eFHNtesttype1:
        case eFHNMONO:
        case eAlievPanfilov:
            if(m_explicitDiffusion)
            {
                out << "\tDiffusion Advancement   : Explicit" <<endl;
            }
            else
            {
                out << "\tDiffusion Advancement   : Implicit" <<endl;
            }
            if(m_explicitReaction)
            {
                out << "\tReaction Advancement    : Explicit" <<endl;
            }
            else
            {
                out << "\tReaction Advancement    : Implicit" <<endl;
            }
            out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
            EquationSystem::TimeParamSummary(out);
            break;
        default:
            ASSERTL0(false, "Unsupported equation type.");
            break;
        }
        cout << "=======================================================================" << endl;
    }

} //end of namespace
