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

    std::string NavierStokesAdvection::navierStokesAdvectionTypeLookupIds[2] = {
        LibUtilities::SessionReader::RegisterEnumValue("SPECTRALHPDEALIASING",
            "True", 0),
        LibUtilities::SessionReader::RegisterEnumValue("SPECTRALHPDEALIASING",
            "False", 1)};

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
        m_homogen_dealiasing = pSession->DefinesSolverInfo("dealiasing");

        pSession->MatchSolverInfo("SPECTRALHPDEALIASING","True",m_specHP_dealiasing,false);
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
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
    {
        int nqtot            = fields[0]->GetTotPoints();
        ASSERTL1(nConvectiveFields == inarray.size(),"Number of convective fields and Inarray are not compatible");

        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = advVel.size();
        Array<OneD, Array<OneD, NekDouble> > AdvVel   (advVel.size());

        Array<OneD, Array<OneD, NekDouble> > velocity(ndim);
        for(int i = 0; i < ndim; ++i)
        {
            if(fields[i]->GetWaveSpace() && !m_SingleMode && !m_HalfMode &&
                !m_homogen_dealiasing)
            {
                velocity[i] = Array<OneD, NekDouble>(nqtot,0.0);
                fields[i]->HomogeneousBwdTrans(advVel[i],velocity[i]);
            }
            else
            {
                velocity[i] = advVel[i];
            }
        }

        int nPointsTot = fields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2,wkSp;

        NekDouble OneDptscale = 1.5; // factor to rescale 1d points in dealiasing

        if(m_specHP_dealiasing)
        {
            // Get number of points to dealias a quadratic non-linearity
            nPointsTot = fields[0]->Get1DScaledTotPoints(OneDptscale);
        }

        // interpolate Advection velocity
        if(m_specHP_dealiasing) // interpolate advection field to higher space.
        {
            for(int i = 0; i < ndim; ++i)
            {
                AdvVel[i] = Array<OneD, NekDouble> (nPointsTot);
                // interpolate infield to 3/2 dimension
                fields[0]->PhysInterp1DScaled(OneDptscale,velocity[i],AdvVel[i]);
            }
        }
        else
        {
            for(int i = 0; i < ndim; ++i)
            {
                AdvVel[i] = velocity[i];
            }
        }

        wkSp = Array<OneD, NekDouble> (nPointsTot);

        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            grad0 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
            for(int n = 0; n < nConvectiveFields; ++n)
            {
                fields[0]->PhysDeriv(inarray[n],grad0);
                if(m_specHP_dealiasing)  // interpolate gradient field
                {
                    Array<OneD, NekDouble> Outarray(nPointsTot);
                    fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                    Vmath::Vmul (nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                    // Galerkin project solution back to origianl spac
                    fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,outarray[n]);
                }
                else
                {
                    Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,outarray[n],1);
                }
            }
            break;
        case 2:
            grad0 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
            grad1 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
            for(int n = 0; n < nConvectiveFields; ++n)
            {
                fields[0]->PhysDeriv(inarray[n],grad0,grad1);

                if(m_specHP_dealiasing)  // interpolate gradient field
                {
                    Array<OneD, NekDouble> Outarray(nPointsTot);
                    fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                    Vmath::Vmul (nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);
                    fields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                    Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[1],1,Outarray,1,Outarray,1);
                    // Galerkin project solution back to original space
                    fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,Outarray,outarray[n]);
                }
                else
                {
                    Vmath::Vmul (nPointsTot,grad0,1,AdvVel[0],1,outarray[n],1);
                    Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,outarray[n],1,outarray[n],1);
                }
            }
            break;
        case 3:
            if(m_homogen_dealiasing == true && m_specHP_dealiasing == true)
            {
                Array<OneD, Array<OneD, NekDouble> > grad (ndim);
                Array<OneD, Array<OneD, NekDouble> > gradScaled (ndim*nConvectiveFields);
                Array<OneD, Array<OneD, NekDouble> > Outarray (nConvectiveFields);
                for (int i = 0; i < ndim; i++)
                {
                    grad[i] = Array<OneD, NekDouble>(fields[0]->GetNpoints());
                }
                for (int i = 0; i < ndim*nConvectiveFields; i++)
                {
                    gradScaled[i] = Array<OneD, NekDouble>(nPointsTot);
                }
                for (int i = 0; i < nConvectiveFields; i++)
                {
                    Outarray[i] = Array<OneD, NekDouble>(nPointsTot);
                }

                for (int n = 0; n < nConvectiveFields; n++)
                {
                    fields[0]->PhysDeriv(inarray[n],grad[0],grad[1],grad[2]);
                    for (int i = 0; i < ndim; i++)
                    {
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad[i],
                                                      gradScaled[n*ndim+i]);
                    }
                }

                fields[0]->DealiasedDotProd(AdvVel,gradScaled,Outarray);

                for (int n = 0; n < nConvectiveFields; n++)
                {
                    fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,
                                    Outarray[n],outarray[n]);
                }
            }
            else if(m_homogen_dealiasing == true && m_specHP_dealiasing == false)
            {
                Array<OneD, Array<OneD, NekDouble> > grad (ndim*nConvectiveFields);
                Array<OneD, Array<OneD, NekDouble> > Outarray (nConvectiveFields);
                for (int i = 0; i < ndim*nConvectiveFields; i++)
                {
                    grad[i] = Array<OneD, NekDouble>(nPointsTot);
                }
                for (int i = 0; i < nConvectiveFields; i++)
                {
                    Outarray[i] = Array<OneD, NekDouble>(nPointsTot);
                }

                for (int n = 0; n < nConvectiveFields; n++)
                {
                    fields[0]->PhysDeriv(inarray[n],grad[n*ndim+0],
                                                    grad[n*ndim+1],
                                                    grad[n*ndim+2]);
                }

                fields[0]->DealiasedDotProd(AdvVel,grad,outarray);
            }
            else
            {
                grad0 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
                grad1 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
                grad2 = Array<OneD, NekDouble> (fields[0]->GetNpoints());
                for(int n = 0; n < nConvectiveFields; ++n)
                {
                    if (fields[0]->GetWaveSpace() == true &&
                        fields[0]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        if (n < ndim)
                        {
                            // take d/dx, d/dy  gradients in physical Fourier space
                            fields[0]->PhysDeriv(velocity[n],grad0,grad1);
                        }
                        else
                        {
                            fields[0]->HomogeneousBwdTrans(inarray[n],wkSp);
                            fields[0]->PhysDeriv(wkSp,grad0,grad1);
                        }
                        // Take d/dz derivative using wave space field
                        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                              inarray[n],
                                              outarray[n]);
                        fields[0]->HomogeneousBwdTrans(outarray[n],grad2);
                    }
                    else if (fields[0]->GetWaveSpace() == true &&
                             fields[0]->GetExpType() == MultiRegions::e3DH2D)
                    {
                        if (n < ndim)
                        {
                            // take d/dx,  gradients in physical Fourier space
                            fields[0]->PhysDeriv(velocity[n],grad0);
                        }
                        else
                        {
                            fields[0]->HomogeneousBwdTrans(inarray[n],wkSp);
                            fields[0]->PhysDeriv(wkSp,grad0);
                        }
                        // Take d/dy derivative using wave space field
                        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],inarray[n],
                                              outarray[n]);
                        fields[0]->HomogeneousBwdTrans(outarray[n],grad1);
                        // Take d/dz derivative using wave space field
                        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],inarray[n],
                                              outarray[n]);
                        fields[0]->HomogeneousBwdTrans(outarray[n],grad2);
                    }
                    else
                    {
                        fields[0]->PhysDeriv(inarray[n],grad0,grad1,grad2);
                    }
                    if(m_specHP_dealiasing) //interpolate spectral/hp gradient field
                    {
                        Array<OneD, NekDouble> Outarray(nPointsTot);
                        fields[0]->PhysInterp1DScaled(OneDptscale,grad0,wkSp);
                        Vmath::Vmul(nPointsTot,wkSp,1,AdvVel[0],1,Outarray,1);

                        fields[0]->PhysInterp1DScaled(OneDptscale,grad1,wkSp);
                        Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[1],1,Outarray,1,
                                     Outarray,1);

                        fields[0]->PhysInterp1DScaled(OneDptscale,grad2,wkSp);
                        Vmath::Vvtvp(nPointsTot,wkSp,1,AdvVel[2],1,Outarray,1,
                                     Outarray,1);
                        fields[0]->PhysGalerkinProjection1DScaled(OneDptscale,
                                     Outarray,outarray[n]);
                    }
                    else
                    {
                        Vmath::Vmul(nPointsTot,grad0,1,AdvVel[0],1,outarray[n],1);
                        Vmath::Vvtvp(nPointsTot,grad1,1,AdvVel[1],1,outarray[n],1,
                                     outarray[n],1);
                        Vmath::Vvtvp(nPointsTot,grad2,1,AdvVel[2],1,outarray[n],1,
                                     outarray[n],1);
                    }

                    if(fields[0]->GetWaveSpace() == true)
                    {
                        fields[0]->HomogeneousFwdTrans(outarray[n],outarray[n]);
                    }
                }
            }
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }

        for(int n = 0; n < nConvectiveFields; ++n)
        {
            Vmath::Neg(nqtot,outarray[n],1);
        }

    }

} //end of namespace

