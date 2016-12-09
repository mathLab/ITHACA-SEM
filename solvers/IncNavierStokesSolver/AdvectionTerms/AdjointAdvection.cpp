///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvection.cpp
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
// Description: Evaluation of the adjoint advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.h>

using namespace std;

namespace Nektar
{

string AdjointAdvection::className = SolverUtils
        ::GetAdvectionFactory().RegisterCreatorFunction("Adjoint",
                                                    AdjointAdvection::create);

/**
 *
 */
AdjointAdvection::AdjointAdvection():
    LinearisedAdvection()
{
}

AdjointAdvection::~AdjointAdvection()
{
}

void AdjointAdvection::v_Advect(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const NekDouble                                   &time,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    ASSERTL1(nConvectiveFields == inarray.num_elements(),
             "Number of convective fields and Inarray are not compatible");

    int nPointsTot  = fields[0]->GetNpoints();
    int ndim        = advVel.num_elements();
    int nBaseDerivs = (m_halfMode || m_singleMode) ? 2 : m_spacedim;
    int nDerivs     = (m_halfMode) ? 2 : m_spacedim;

    Array<OneD, Array<OneD, NekDouble> > velocity(ndim);
    for(int i = 0; i < ndim; ++i)
    {
        if(fields[i]->GetWaveSpace() && !m_singleMode && !m_halfMode)
        {
            velocity[i] = Array<OneD, NekDouble>(nPointsTot,0.0);
            fields[i]->HomogeneousBwdTrans(advVel[i],velocity[i]);
        }
        else
        {
            velocity[i] = advVel[i];
        }
    }

    Array<OneD, Array<OneD, NekDouble> > grad (nDerivs);
    for( int i = 0; i < nDerivs; ++i)
    {
        grad[i] = Array<OneD, NekDouble> (nPointsTot);
    }

    // Evaluation of the base flow for periodic cases
    if (m_slices > 1)
    {
        for (int i = 0; i < ndim; ++i)
        {
            UpdateBase(m_slices, m_interp[i], m_baseflow[i],
                       time, m_period);
            UpdateGradBase(i, fields[i]);
        }
    }

    //Evaluate the linearised advection term
    for( int i = 0; i < ndim; ++i)
    {
        // Calculate gradient
        switch(nDerivs)
        {
            case 1:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0]);
            }
            break;
            case 2:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1]);
            }
            break;
            case 3:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1], grad[2]);
                if(m_multipleModes)
                {
                    // transform gradients into physical Fourier space
                    fields[i]->HomogeneousBwdTrans(grad[0], grad[0]);
                    fields[i]->HomogeneousBwdTrans(grad[1], grad[1]);
                    fields[i]->HomogeneousBwdTrans(grad[2], grad[2]);
                }
            }
            break;
        }

        // Calculate -U_j du'_i/dx_j
        Vmath::Vmul(nPointsTot,grad[0], 1, m_baseflow[0], 1, outarray[i], 1);
        for( int j = 1; j < nDerivs; ++j)
        {
            Vmath::Vvtvp(nPointsTot,grad[j], 1,
                                    m_baseflow[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }
        Vmath::Neg(nPointsTot,outarray[i],1);

        // Add u'_j U_j/ dx_i
        int lim = (m_halfMode) ? 2 : ndim;
        if ( (m_halfMode || m_singleMode) && i==2)
        {
            lim = 0;
        }
        for( int j = 0; j < lim; ++j)
        {
            Vmath::Vvtvp(nPointsTot,m_gradBase[j*nBaseDerivs + i], 1,
                                    velocity[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }

        if(m_multipleModes)
        {
            fields[i]->HomogeneousFwdTrans(outarray[i],outarray[i]);
        }
        Vmath::Neg(nPointsTot,outarray[i],1);
    }
}

} //end of namespace

