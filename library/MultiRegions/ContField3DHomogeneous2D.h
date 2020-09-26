///////////////////////////////////////////////////////////////////////////////
//
// File ContField3DHomogeneous2D.h
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
// Description: Field definition in three-dimensions for a continuous
// expansion with two homogeneous directions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3DHOMO2D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3DHOMO2D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3DHomogeneous2D: public DisContField3DHomogeneous2D
        {
        public:
            MULTI_REGIONS_EXPORT ContField3DHomogeneous2D();

            MULTI_REGIONS_EXPORT ContField3DHomogeneous2D(
                         const LibUtilities::SessionReaderSharedPtr &pSession,
                         const LibUtilities::BasisKey &HomoBasis_y,
                         const LibUtilities::BasisKey &HomoBasis_z,
                         const NekDouble lhom_y,
                         const NekDouble lhom_z,
                         const bool useFFT,
                         const bool dealiasing,
                         const SpatialDomains::MeshGraphSharedPtr &graph1D,
                         const std::string &variable,
                         const Collections::ImplementationType ImpType
                         = Collections::eNoImpType);
            
            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField3DHomogeneous2D(
                                 const ContField3DHomogeneous2D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ContField3DHomogeneous2D();
            
        protected:

        private:
            Array<OneD, NekDouble> m_contCoeffs;

            // virtual functions

            virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            /// Template method virtual forwarded for LocalToGlobal()
            virtual void v_LocalToGlobal(bool useComm);

            /// Template method virtual forwarded for GlobalToLocal()
            virtual void v_GlobalToLocal(void);

            /// Solves the three-dimensional Helmholtz equation, subject to the
            /// boundary conditions specified.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);

            virtual void v_ClearGlobalLinSysManager(void);

        };

        typedef std::shared_ptr<ContField3DHomogeneous2D>  
            ContField3DHomogeneous2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3DHOMO1D_H
