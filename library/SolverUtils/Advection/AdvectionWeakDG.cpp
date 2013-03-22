///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG::type = GetAdvectionFactory().
            RegisterCreatorFunction("WeakDG", AdvectionWeakDG::create);
        
        AdvectionWeakDG::AdvectionWeakDG()
        {
        }
        
        void AdvectionWeakDG::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j;
            int nVelDim         = fields[0]->GetCoordim(0);
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(nConvectiveFields);
            //Array<OneD, Array<OneD, NekDouble> > fluxvector(nVelDim);
            Array<OneD, Array<OneD, NekDouble> > tmp       (nConvectiveFields);

            ASSERTL1(m_riemann, 
                     "Riemann solver must be provided for AdvectionWeakDG.");

            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] = Array<OneD, Array<OneD, NekDouble> >(nVelDim);
                for(j = 0; j < nVelDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            m_fluxVector(inarray, fluxvector);

            // Get the advection part (without numerical flux)
            for(i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD,NekDouble>(nCoeffs, 0.0);
                
                for (j = 0; j < nVelDim; ++j)
                {
                    fields[i]->IProductWRTDerivBase(j, fluxvector[i][j], 
                                                       outarray[i]);
                    Vmath::Vadd(nCoeffs, outarray[i], 1, tmp[i], 1, tmp[i], 1);
                }
            }
            
            // Store forwards/backwards space along trace space.
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
            
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot);
                numflux[i] = Array<OneD, NekDouble>(nTracePointsTot);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
            
            m_riemann->Solve(Fwd, Bwd, numflux);
            
            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }
    }
}
