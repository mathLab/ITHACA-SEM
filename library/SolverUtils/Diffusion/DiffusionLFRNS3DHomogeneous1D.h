///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNS3DHomogeneous1D.h
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
///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNSNS.h
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
// Description: LFRNS diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONLFRNS3DHomogeneous1D
#define NEKTAR_SOLVERUTILS_DIFFUSIONLFRNS3DHomogeneous1D

#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/Diffusion/DiffusionLFRNS.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionLFRNS3DHomogeneous1D : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionLFRNS3DHomogeneous1D(diffType));
            }
            
            static std::string                   type[];
            
        protected:
            DiffusionLFRNS3DHomogeneous1D(std::string diffType);
            
            SolverUtils::DiffusionSharedPtr              m_planeDiff;
            
            std::string m_diffType;
            
            int                                          spaceDim;
            int                                          nPointsTot;
            int                                          nCoeffs;
            int                                          num_planes;
            int                                          nPointsTot_plane;
            int                                          nCoeffs_plane;
            int                                          nVariables;
            int                                          i, j, k;
            
            Array<OneD, unsigned int>                    planes;
            Array<OneD, Array<OneD, NekDouble> >         fluxvector;
            Array<OneD, Array<OneD, NekDouble> >         derivatives01_homo;
            Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                         fields_plane;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                         inarray_plane;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                         derivatives01_plane;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                         outarray_plane;
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                         fluxvector_homo;
            Array<OneD, Array<OneD, NekDouble> >         outarray_homo;
            
            DiffusionLFRNSSharedPtr                      diffLFRNS;
            
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            
            virtual void v_Diffuse(
                const int                                     nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
            
        }; 
    }
}

#endif
