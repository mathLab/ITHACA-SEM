///////////////////////////////////////////////////////////////////////////////
//
// File: ArtificialDiffusion.h
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
// Description: Abstract base class for compressible solver artificial diffusion
//              used for shock capturing artificial diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_BASE
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_BASE

#include <string>

#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
// Forward declaration
class ArtificialDiffusion;

/// A shared pointer to a artificial diffusion object
typedef std::shared_ptr<ArtificialDiffusion> ArtificialDiffusionSharedPtr;

/// Declaration of the artificial diffusion factory
typedef LibUtilities::NekFactory<std::string, ArtificialDiffusion,
        const LibUtilities::SessionReaderSharedPtr&,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&,
        const int > ArtificialDiffusionFactory;

/// Declaration of the artificial diffusion factory singleton
ArtificialDiffusionFactory& GetArtificialDiffusionFactory();

/**
 * @class ArtificialDiffusion
 * @brief Encapsulates the artificial diffusion used in shock capture
 */
class ArtificialDiffusion
{
    public:
        virtual ~ArtificialDiffusion() {}

        /// Apply the artificial diffusion
        void DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray);
        
        /// Apply the artificial diffusion the outarray is in coeff space
        void DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray);

        /// Calculate the artificial viscosity
        void GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu);

        /// Set h/p scaling
        void SetElmtHP(const Array<OneD, NekDouble> &hOverP);

    protected:
        /// Session reader
        LibUtilities::SessionReaderSharedPtr        m_session;
        /// Array of fields
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
        /// Auxiliary object to convert variables
        VariableConverterSharedPtr                  m_varConv;
        /// LDG Diffusion operator
        SolverUtils::DiffusionSharedPtr             m_diffusion;
        /// Constant scaling
        NekDouble                                   m_mu0;
        /// h/p scaling
        Array<OneD, NekDouble>                      m_hOverP;

        /// Constructor
        ArtificialDiffusion(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim);

        virtual void v_DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray);
        
        virtual void v_DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> >             &outarray);

        virtual void v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)=0;

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                &viscousTensor);


};
}

#endif
