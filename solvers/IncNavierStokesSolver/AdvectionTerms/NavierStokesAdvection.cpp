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

using namespace std;

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


    void NavierStokesAdvection::v_Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray,
        const NekDouble                                   &time)
    {
        int nqtot            = fields[0]->GetTotPoints();
        ASSERTL1(nConvectiveFields == inarray.num_elements(),"Number of convective fields and Inarray are not compatible");

        Array<OneD, NekDouble > Deriv(nqtot*nConvectiveFields);

        for(int n = 0; n < nConvectiveFields; ++n)
        {
            // use dimension of Velocity vector to dictate dimension of operation
            int ndim       = advVel.num_elements();
            Array<OneD, Array<OneD, NekDouble> > AdvVel   (advVel.num_elements());
            Array<OneD, NekDouble> Outarray;


            int nPointsTot = fields[0]->GetNpoints();
            Array<OneD, NekDouble> grad0,grad1,grad2,wkSp;

            NekDouble OneDptscale = 1.5; // factor to rescale 1d points in dealiasing

            if(m_specHP_dealiasing)
            {
                // Get number of points to dealias a quadratic non-linearity
                nPointsTot = fields[0]->Get1DScaledTotPoints(OneDptscale);
            }

            grad0 = Array<OneD, NekDouble> (nPointsTot);

            // interpolate Advection velocity
            int nadv = advVel.num_elements();
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
                    fields[0]->PhysInterp1DScaled(OneDptscale,advVel[i],AdvVel[i]);
                }

                Outarray = AdvVel[nadv-1] + nPointsTot;
            }
            else
            {
                for(int i = 0; i < nadv; ++i)
                {
                    AdvVel[i] = advVel[i];
                }

                Outarray = outarray[n];
            }

            wkSp = Array<OneD, NekDouble> (nPointsTot);


            // Evaluate V\cdot Grad(u)
            switch(ndim)
            {
            case 1:
                fields[0]->PhysDeriv(inarray[n],grad0);
                Vmath::Vmul(nPointsTot,grad0,1,advVel[0],1,outarray[n],1);
                break;
            case 2:
                {
                    grad1 = Array<OneD, NekDouble> (nPointsTot);
                    fields[0]->PhysDeriv(inarray[n],grad0,grad1);

                    if(m_specHP_dealiasing)  // interpolate gradient field
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                        Vmath::Vcopy(nPointsTot,wkSp,1,grad0,1);
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                        Vmath::Vcopy(nPointsTot,wkSp,1,grad1,1);
                    }

                    Vmath::Vmul (nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                    Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,Outarray,1,Outarray,1);

                    if(m_specHP_dealiasing) // Galerkin project solution back to origianl space
                    {
                        fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,outarray[n]);
                    }

                }
                break;
            case 3:
                grad1 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
                grad2 = Array<OneD, NekDouble> (fields[0]->GetNpoints());

                if(fields[0]->GetWaveSpace() == false && m_homogen_dealiasing == true )
                {
                    ASSERTL0(m_specHP_dealiasing == false,"Spectral/hp element dealaising is not set up for this option");

                    fields[0]->PhysDeriv(inarray[n],grad0,grad1,grad2);

                    fields[0]->DealiasedProd(advVel[0],grad0,grad0,m_CoeffState);
                    fields[0]->DealiasedProd(advVel[1],grad1,grad1,m_CoeffState);
                    fields[0]->DealiasedProd(advVel[2],grad2,grad2,m_CoeffState);
                    Vmath::Vadd(nPointsTot,grad0,1,grad1,1,outarray[n],1);
                    Vmath::Vadd(nPointsTot,grad2,1,outarray[n],1,outarray[n],1);
                }
                else if(fields[0]->GetWaveSpace() == true && m_homogen_dealiasing == false)
                {
                    // take d/dx, d/dy  gradients in physical Fourier space
                    fields[0]->PhysDeriv(advVel[n],grad0,grad1);

                    // Take d/dz derivative using wave space field
                    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],inarray[n],
                                          outarray[n]);
                    fields[0]->HomogeneousBwdTrans(outarray[n],grad2);

                    if(m_specHP_dealiasing) //interpolate spectral/hp gradient field
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                        Vmath::Vmul(nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                    }
                    else
                    {
                        Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                    }

                    if(m_specHP_dealiasing) //interpolate spectral/hp gradient field
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
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
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad2,wkSp);
                        Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[2],1,Outarray,1,Outarray,1);
                        fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,grad2);
                        fields[0]->HomogeneousFwdTrans(grad2,outarray[n]);
                    }
                    else
                    {
                        Vmath::Vvtvp(nPointsTot,grad2,1,AdvVel[2],1,Outarray,1,grad0,1);
                        fields[0]->HomogeneousFwdTrans(grad0,outarray[n]);
                    }
                }
                else if(fields[0]->GetWaveSpace() == false && m_homogen_dealiasing == false)
                {

                    fields[0]->PhysDeriv(inarray[n],grad0,grad1,grad2);

                    if(m_specHP_dealiasing) //interpolate spectral/hp gradient field
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                        Vmath::Vmul(nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                    }
                    else
                    {
                        Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,Outarray,1);
                    }


                    if(m_specHP_dealiasing) //interpolate spectral/hp gradient field
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
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
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad2,wkSp);
                        Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[2],1,Outarray,1,Outarray,1);
                        fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,outarray[n]);
                    }
                    else
                    {
                        Vmath::Vvtvp(nPointsTot,grad2,1,AdvVel[2],1,Outarray,1,outarray[n],1);
                    }
                }
                else if(fields[0]->GetWaveSpace() == true && m_homogen_dealiasing == true)
                {
                    ASSERTL0(m_specHP_dealiasing == false,"Spectral/hp element dealaising is not set up for this option");

                    fields[0]->PhysDeriv(inarray[n],grad0,grad1,grad2);

                    fields[0]->HomogeneousBwdTrans(grad0, outarray[n]);
                    fields[0]->DealiasedProd(advVel[0], outarray[n], grad0,
                                              m_CoeffState);

                    fields[0]->HomogeneousBwdTrans(grad1,outarray[n]);
                    fields[0]->DealiasedProd(advVel[1], outarray[n], grad1,
                                              m_CoeffState);

                    fields[0]->HomogeneousBwdTrans(grad2,outarray[n]);
                    fields[0]->DealiasedProd(advVel[2], outarray[n], grad2,
                                              m_CoeffState);

                    Vmath::Vadd(nPointsTot, grad0, 1, grad1, 1, grad0, 1);
                    Vmath::Vadd(nPointsTot, grad0, 1, grad2, 1, grad0, 1);

                    fields[0]->HomogeneousFwdTrans(grad0,outarray[n]);
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

            Vmath::Neg(nqtot,outarray[n],1);
        }

    }

} //end of namespace

