///////////////////////////////////////////////////////////////////////////////
//
// File: CFSBndCond.h
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
// Description: Abstract base class for compressible solver boundary conditions.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_BNDCOND_CFSBNDCOND
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_BNDCOND_CFSBNDCOND

#include <string>

#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
//  Forward declaration
class CFSBndCond;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<CFSBndCond> CFSBndCondSharedPtr;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<std::string, CFSBndCond,
        const LibUtilities::SessionReaderSharedPtr&,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&,
        const Array<OneD, Array<OneD, NekDouble> >&,
        const int,
        const int,
        const int> CFSBndCondFactory;

/// Declaration of the boundary condition factory singleton
CFSBndCondFactory& GetCFSBndCondFactory();

/**
 * @class CFSBndCond
 * @brief Encapsulates the user-defined boundary conditions for compressible
 *        flow solver.
 */
class CFSBndCond
{
    public:
        virtual ~CFSBndCond() {}

        /// Apply the boundary condition
        void Apply(
            Array<OneD, Array<OneD, NekDouble> >               &Fwd,
            Array<OneD, Array<OneD, NekDouble> >               &physarray,
            const NekDouble                                    &time = 0);

        /// Apply the Weight of boundary condition
        void ApplyBwdWeight()
        {
            v_ApplyBwdWeight();
        }

    protected:
        /// Session reader
        LibUtilities::SessionReaderSharedPtr m_session;
        /// Array of fields
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
        /// Trace normals
        Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
        /// Space dimension
        int m_spacedim;
        /// Auxiliary object to convert variables
        VariableConverterSharedPtr           m_varConv;
        /// Weight for average calculation of diffusion term
        NekDouble m_diffusionAveWeight;

        /// Parameters of the flow
        NekDouble m_gamma;
        NekDouble m_rhoInf;
        NekDouble m_pInf;
        NekDouble m_pOut;
        Array<OneD, NekDouble> m_velInf;

        /// Id of the boundary region
        int       m_bcRegion;
        /// Offset
        int       m_offset;

        /// Constructor
        CFSBndCond(const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const Array<OneD, Array<OneD, NekDouble> >&       pTraceNormals,
                const int pSpaceDim,
                const int bcRegion,
                const int cnt);

        virtual void v_Apply(
            Array<OneD, Array<OneD, NekDouble> >               &Fwd,
            Array<OneD, Array<OneD, NekDouble> >               &physarray,
            const NekDouble                                    &time)=0;

        virtual void v_ApplyBwdWeight();
};
}

#endif
