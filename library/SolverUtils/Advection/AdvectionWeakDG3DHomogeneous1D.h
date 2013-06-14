///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG3DHomogeneous1D.h
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
// Description: Weak DG advection 3DHomogeneous1D class class .
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTIONWEAKDG3DHomogeneous1D
#define NEKTAR_SOLVERUTILS_ADVECTIONWEAKDG3DHomogeneous1D

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class AdvectionWeakDG3DHomogeneous1D : public Advection
        {
        public:
            static AdvectionSharedPtr create(std::string advType)
            {
                return AdvectionSharedPtr(new AdvectionWeakDG3DHomogeneous1D());
            }
                        
            static std::string type;
            
        protected:
            AdvectionWeakDG3DHomogeneous1D();
            
            SolverUtils::AdvectionSharedPtr m_planeAdv;
            
            virtual void v_Advect(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
        };
    }
}
    
#endif
