///////////////////////////////////////////////////////////////////////////////
//
// File SkewSymmetricAdvection.cpp
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

#include <IncNavierStokesSolver/AdvectionTerms/SkewSymmetricAdvection.h>

using namespace std;

namespace Nektar
{
string SkewSymmetricAdvection::className
    = SolverUtils::GetAdvectionFactory().RegisterCreatorFunction(
            "SkewSymmetric",
            SkewSymmetricAdvection::create);

/**
 *
 */
SkewSymmetricAdvection::SkewSymmetricAdvection():
    Advection()

{
}


/**
 *
 */
SkewSymmetricAdvection::~SkewSymmetricAdvection()
{
}


/**
 *
 */
void SkewSymmetricAdvection::v_InitObject(
                LibUtilities::SessionReaderSharedPtr        pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    Advection::v_InitObject(pSession, pFields);

    m_homogen_dealiasing = pSession->DefinesSolverInfo("dealiasing");
    pSession->MatchSolverInfo("ModeType","SingleMode",m_SingleMode,false);
    pSession->MatchSolverInfo("ModeType","HalfMode",m_HalfMode,false);
}


/**
 *
 */
void SkewSymmetricAdvection::v_Advect(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
          Array<OneD, Array<OneD, NekDouble> >        &outarray,
    const NekDouble                                   &time,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    // use dimension of Velocity vector to dictate dimension of operation
    int ndim             = advVel.size();
    int nqtot            = fields[0]->GetTotPoints();
    ASSERTL1(nConvectiveFields == inarray.size(),"Number of convective fields and Inarray are not compatible");

    Array<OneD, Array<OneD, NekDouble> > velocity(ndim);
    for(int i = 0; i < ndim; ++i)
    {
        if(fields[i]->GetWaveSpace() && !m_SingleMode && !m_HalfMode)
        {
            velocity[i] = Array<OneD, NekDouble>(nqtot,0.0);
            fields[i]->HomogeneousBwdTrans(advVel[i],velocity[i]);
        }
        else
        {
            velocity[i] = advVel[i];
        }
    }

    for(int n = 0; n < nConvectiveFields; ++n)
    {
        // ToDo: here we should add a check that V has right dimension

        int nPointsTot = fields[0]->GetNpoints();
        Array<OneD, NekDouble> gradV0,gradV1,gradV2, tmp, Up;

        gradV0   = Array<OneD, NekDouble> (nPointsTot);
        tmp = Array<OneD, NekDouble> (nPointsTot);

        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            fields[0]->PhysDeriv(inarray[n],gradV0);
            Vmath::Vmul(nPointsTot,gradV0,1,velocity[0],1,outarray[n],1);
            Vmath::Vmul(nPointsTot,inarray[n],1,velocity[0],1,gradV0,1);
            fields[0]->PhysDeriv(gradV0,tmp);
            Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
            Vmath::Smul(nPointsTot,0.5,outarray[n],1,outarray[n],1);
            break;
        case 2:
            gradV1 = Array<OneD, NekDouble> (nPointsTot);
            fields[0]->PhysDeriv(inarray[n],gradV0,gradV1);
            Vmath::Vmul (nPointsTot,gradV0,1,velocity[0],1,outarray[n],1);
            Vmath::Vvtvp(nPointsTot,gradV1,1,velocity[1],1,outarray[n],1,outarray[n],1);
            Vmath::Vmul(nPointsTot,inarray[n],1,velocity[0],1,gradV0,1);
            Vmath::Vmul(nPointsTot,inarray[n],1,velocity[1],1,gradV1,1);
            fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
            Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
            fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
            Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
            Vmath::Smul(nPointsTot,0.5,outarray[n],1,outarray[n],1);
            break;
        case 3:
            gradV1 = Array<OneD, NekDouble> (nPointsTot);
            gradV2 = Array<OneD, NekDouble> (nPointsTot);

            fields[0]->PhysDeriv(inarray[n],gradV0,gradV1,gradV2);

            //outarray[n] = 1/2(u*du/dx + v*du/dy + w*du/dz + duu/dx + duv/dy + duw/dz)

            if(m_homogen_dealiasing == true && fields[0]->GetWaveSpace() == false)
            {
                fields[0]->DealiasedProd(velocity[0],gradV0,gradV0);
                fields[0]->DealiasedProd(velocity[1],gradV1,gradV1);
                fields[0]->DealiasedProd(velocity[2],gradV2,gradV2);
                Vmath::Vadd(nPointsTot,gradV0,1,gradV1,1,outarray[n],1);
                Vmath::Vadd(nPointsTot,gradV2,1,outarray[n],1,outarray[n],1);
                fields[0]->DealiasedProd(inarray[n],velocity[0],gradV0);
                fields[0]->DealiasedProd(inarray[n],velocity[1],gradV1);
                fields[0]->DealiasedProd(inarray[n],velocity[2],gradV2);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                Vmath::Smul(nPointsTot,0.5,outarray[n],1,outarray[n],1);
            }
            else if(fields[0]->GetWaveSpace() == true && m_homogen_dealiasing == false)
            {
                Up = Array<OneD, NekDouble> (nPointsTot);
                //vector reused to avoid even more memory requirements
                //names may be misleading
                fields[0]->HomogeneousBwdTrans(gradV0,tmp);
                Vmath::Vmul(nPointsTot,tmp,1,velocity[0],1,outarray[n],1); // + u*du/dx
                fields[0]->HomogeneousBwdTrans(gradV1,tmp);
                Vmath::Vvtvp(nPointsTot,tmp,1,velocity[1],1,outarray[n],1,outarray[n],1);// + v*du/dy
                fields[0]->HomogeneousBwdTrans(gradV2,tmp);
                Vmath::Vvtvp(nPointsTot,tmp,1,velocity[2],1,outarray[n],1,outarray[n],1);// + w*du/dz

                fields[0]->HomogeneousBwdTrans(inarray[n],Up);
                Vmath::Vmul(nPointsTot,Up,1,velocity[0],1,gradV0,1);
                Vmath::Vmul(nPointsTot,Up,1,velocity[1],1,gradV1,1);
                Vmath::Vmul(nPointsTot,Up,1,velocity[2],1,gradV2,1);

                fields[0]->SetWaveSpace(false);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);//duu/dx
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);//duv/dy
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);//duw/dz
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->SetWaveSpace(true);

                Vmath::Smul(nPointsTot,0.5,outarray[n],1,tmp,1);
                fields[0]->HomogeneousFwdTrans(tmp,outarray[n]);
            }
            else if(fields[0]->GetWaveSpace() == false && m_homogen_dealiasing == false)
            {
                Vmath::Vmul(nPointsTot,gradV0,1,velocity[0],1,outarray[n],1);
                Vmath::Vvtvp(nPointsTot,gradV1,1,velocity[1],1,outarray[n],1,outarray[n],1);
                Vmath::Vvtvp(nPointsTot,gradV2,1,velocity[2],1,outarray[n],1,outarray[n],1);
                Vmath::Vmul(nPointsTot,inarray[n],1,velocity[0],1,gradV0,1);
                Vmath::Vmul(nPointsTot,inarray[n],1,velocity[1],1,gradV1,1);
                Vmath::Vmul(nPointsTot,inarray[n],1,velocity[2],1,gradV2,1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);
                Vmath::Vadd(nPointsTot,tmp,1,outarray[n],1,outarray[n],1);
                Vmath::Smul(nPointsTot,0.5,outarray[n],1,outarray[n],1);
            }
            else
            {
                ASSERTL0(false, "Dealiasing is not allowed in combination "
                                "with the Skew-Symmetric advection form for "
                                "efficiency reasons.");
            }
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }

        Vmath::Neg(nqtot,outarray[n],1);
    }

}

} //end of namespace

