///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNS3DHomogeneous1D.cpp
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
///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFR.cpp
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
// Description: LFR diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLFRNS3DHomogeneous1D.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLFRNS3DHomogeneous1D::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                                                          "LFRDGNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                                          "LFRSDNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                                          "LFRHUNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                                          "LFRcminNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                                          "LFRcinfNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create)};
        
        /**
         * @brief DiffusionLFRNS uses the Flux Reconstruction (FR) approach to
         * compute the diffusion term. The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         *
         * \todo Extension to triangles, tetrahedra and other shapes.
         * (Long term objective)
         */
        DiffusionLFRNS3DHomogeneous1D::DiffusionLFRNS3DHomogeneous1D(std::string diffType):m_diffType(diffType)
        {
        }
        
        /**
         * @brief Initiliase DiffusionLFRNS objects and store them before starting
         * the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics,
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to
         * initialise the objects needed by DiffusionLFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFRNS3DHomogeneous1D::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {            
        }
        
                /**
         * @brief Calculate FR Diffusion for the linear problems
         * using an LDG interface flux.
         *
         */
        void DiffusionLFRNS3DHomogeneous1D::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
        }
                       
    }// close namespace SolverUtils
}// close namespace nektar++
