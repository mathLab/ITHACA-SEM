///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesAdvection.cpp
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
// Description: Evaluation of the Navier Stokes advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.h>

namespace Nektar
{
    string NavierStokesAdvection::className  = SolverUtils::GetAdvectionFactory().RegisterCreatorFunction("Convective", NavierStokesAdvection::create);
    string NavierStokesAdvection::className2 = SolverUtils::GetAdvectionFactory().RegisterCreatorFunction("NonConservative", NavierStokesAdvection::create);
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    NavierStokesAdvection::NavierStokesAdvection():
        Advection()
	
    {
        
    }
    
    NavierStokesAdvection::~NavierStokesAdvection()
    {
    }
    

    void NavierStokesAdvection::v_InitObject(
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
    {
        m_CoeffState = MultiRegions::eLocal;
        m_homogen_dealiasing = pSession->DefinesSolverInfo("dealiasing");

        pSession->MatchSolverInfo("SPECTRALHPDEALIASING","True",m_specHP_dealiasing,false);
        if(m_specHP_dealiasing == false)
        {
            pSession->MatchSolverInfo("SPECTRALHPDEALIASING","On",m_specHP_dealiasing,false);
        }
        pSession->MatchSolverInfo("ModeType","SingleMode",m_SingleMode,false);
        pSession->MatchSolverInfo("ModeType","HalfMode",m_HalfMode,false);

        Advection::v_InitObject(pSession, pFields);
    }

    //Advection function
    
    void NavierStokesAdvection::v_Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray,
        const NekDouble                                   &time)
    {
        int i;
        int nqtot            = fields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> > velocity(nConvectiveFields);

        ASSERTL1(nConvectiveFields == inarray.num_elements(),"Number of convective fields and Inarray are not compatible");
cout << "Advection: " << m_SingleMode << ", " << m_HalfMode << endl;
//        for(i = 0; i < nConvectiveFields; ++i)
//        {
//            if(fields[i]->GetWaveSpace() && !m_SingleMode && !m_HalfMode)
//            {
//                velocity[i] = Array<OneD, NekDouble>(nqtot,0.0);
//                fields[i]->HomogeneousBwdTrans(inarray[i],velocity[i]);
//            }
//            else
//            {
//                velocity[i] = inarray[i];
//            }
//        }
//////////////////
        Array<OneD, NekDouble > Deriv;
        Deriv = Array<OneD, NekDouble> (nqtot*nConvectiveFields);

        for(i=0; i< nConvectiveFields; ++i)
        {
            v_ComputeAdvectionTerm(fields,advVel,inarray[i],outarray[i],i,time,Deriv);
            Vmath::Neg(nqtot,outarray[i],1);
        }

    }

    //Evaluation of the advective terms
    void NavierStokesAdvection::v_ComputeAdvectionTerm(
            const Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pV,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
            NekDouble m_time,
            Array<OneD, NekDouble> &pWk)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = pV.num_elements();
        Array<OneD, Array<OneD, NekDouble> > AdvVel   (pV.num_elements());
        Array<OneD, NekDouble> Outarray;
        
	
        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2,wkSp;
		
        NekDouble OneDptscale = 1.5; // factor to rescale 1d points in dealiasing 

        if(m_specHP_dealiasing)
        {
            // Get number of points to dealias a quadratic non-linearity
            nPointsTot = pFields[0]->Get1DScaledTotPoints(OneDptscale);
        }

        grad0 = Array<OneD, NekDouble> (nPointsTot);

        // interpolate Advection velocity
        int nadv = pV.num_elements();
        if(m_specHP_dealiasing) // interpolate advection field to higher space. 
        {
            AdvVel[0] = Array<OneD, NekDouble> (nPointsTot*(nadv+1));
            for(int i = 0; i < nadv; ++i)
            {
                if(i)
                {
                    AdvVel[i] = AdvVel[i-1]+nPointsTot;
                }
                // interpolate infield to 3/2 dimension
                pFields[0]->PhysInterp1DScaled(OneDptscale,pV[i],AdvVel[i]);
            }
            
            Outarray = AdvVel[nadv-1] + nPointsTot;
        }
        else
        {
            for(int i = 0; i < nadv; ++i)
            {
                AdvVel[i] = pV[i];
            }

            Outarray = pOutarray;
        }

        wkSp = Array<OneD, NekDouble> (nPointsTot);


        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            pFields[0]->PhysDeriv(pU,grad0);
            Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            break;
        case 2:
            {
                grad1 = Array<OneD, NekDouble> (nPointsTot);
                pFields[0]->PhysDeriv(pU,grad0,grad1);

                if(m_specHP_dealiasing)  // interpolate gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                    Vmath::Vcopy(nPointsTot,wkSp,1,grad0,1);
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                    Vmath::Vcopy(nPointsTot,wkSp,1,grad1,1);
                }
                
                Vmath::Vmul (nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,Outarray,1,Outarray,1);

                if(m_specHP_dealiasing) // Galerkin project solution back to origianl space 
                {
                    pFields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,pOutarray); 
                }
                
            }
            break;	 
        case 3:
            grad1 = Array<OneD, NekDouble> (pFields[0]->GetNpoints());
            grad2 = Array<OneD, NekDouble> (pFields[0]->GetNpoints());
            
            if(pFields[0]->GetWaveSpace() == false && m_homogen_dealiasing == true )
            {
                ASSERTL0(m_specHP_dealiasing == false,"Spectral/hp element dealaising is not set up for this option");

                pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);

                pFields[0]->DealiasedProd(pV[0],grad0,grad0,m_CoeffState);
                pFields[0]->DealiasedProd(pV[1],grad1,grad1,m_CoeffState);
                pFields[0]->DealiasedProd(pV[2],grad2,grad2,m_CoeffState);
                Vmath::Vadd(nPointsTot,grad0,1,grad1,1,pOutarray,1);
                Vmath::Vadd(nPointsTot,grad2,1,pOutarray,1,pOutarray,1);
            }
            else if(pFields[0]->GetWaveSpace() == true && m_homogen_dealiasing == false)
            {
                // take d/dx, d/dy  gradients in physical Fourier space
                pFields[0]->PhysDeriv(pV[pVelocityComponent],grad0,grad1);
                
                // Take d/dz derivative using wave space field 
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],pU,
                                      pOutarray);
                pFields[0]->HomogeneousBwdTrans(pOutarray,grad2);
                
                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                    Vmath::Vmul(nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                }
                else
                {
                    Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                }
		
                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                    Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[1],1,Outarray,1,
                                 Outarray,1);
                }
                else
                {
                    Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,Outarray,1,
                                 Outarray,1);
                }
		
                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad2,wkSp);
                    Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[2],1,Outarray,1,Outarray,1);
                    pFields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,grad2); 
                    pFields[0]->HomogeneousFwdTrans(grad2,pOutarray);
                }
                else
                {
                    Vmath::Vvtvp(nPointsTot,grad2,1,AdvVel[2],1,Outarray,1,grad0,1);
                    pFields[0]->HomogeneousFwdTrans(grad0,pOutarray);
                }
            }
            else if(pFields[0]->GetWaveSpace() == false && m_homogen_dealiasing == false) 
            {
                
                pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);
                
                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                    Vmath::Vmul(nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                }
                else
                {
                    Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                }


                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                    Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[1],1,Outarray,1,
                                 Outarray,1);
                }
                else
                {
                    Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,Outarray,1,
                                 Outarray,1);
                }
		
                if(m_specHP_dealiasing) //interpolate spectral/hp gradient field 
                {
                    pFields[0]->PhysInterp1DScaled(OneDptscale,grad2,wkSp);
                    Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[2],1,Outarray,1,Outarray,1);
                    pFields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,pOutarray); 
                }
                else
                {
                    Vmath::Vvtvp(nPointsTot,grad2,1,AdvVel[2],1,Outarray,1,pOutarray,1);
                }
            }
            else if(pFields[0]->GetWaveSpace() == true && m_homogen_dealiasing == true) 
            {
                ASSERTL0(m_specHP_dealiasing == false,"Spectral/hp element dealaising is not set up for this option");

                pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);

                pFields[0]->HomogeneousBwdTrans(grad0, pOutarray);
                pFields[0]->DealiasedProd(pV[0], pOutarray, grad0, 
                                          m_CoeffState);

                pFields[0]->HomogeneousBwdTrans(grad1,pOutarray);
                pFields[0]->DealiasedProd(pV[1], pOutarray, grad1,
                                          m_CoeffState);

                pFields[0]->HomogeneousBwdTrans(grad2,pOutarray);
                pFields[0]->DealiasedProd(pV[2], pOutarray, grad2,
                                          m_CoeffState);

                Vmath::Vadd(nPointsTot, grad0, 1, grad1, 1, grad0, 1);
                Vmath::Vadd(nPointsTot, grad0, 1, grad2, 1, grad0, 1);

                pFields[0]->HomogeneousFwdTrans(grad0,pOutarray);
            }
            else 
            {
                ASSERTL0(false, "Advection term calculation not implented or "
                                "possible with the current problem set up");
            }
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }

} //end of namespace

