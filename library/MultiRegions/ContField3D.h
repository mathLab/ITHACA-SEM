///////////////////////////////////////////////////////////////////////////////
//
// File ContField3D.h
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
// Description: Field definition in three-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContField.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3D: public ContField
        {
        public:
            MULTI_REGIONS_EXPORT ContField3D();

            /// Construct a global continuous field.
            MULTI_REGIONS_EXPORT ContField3D(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const std::string &variable  = "DefaultVar",
                        const bool DeclareCoeffPhysArrays = true,
                        const bool CheckIfSingularSystem = false,
                        const Collections::ImplementationType ImpType
                        = Collections::eNoImpType);

            /// Construct a global continuous field with solution type based on
            /// another field but using a separate input mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const std::string &variable,
                        const bool CheckIfSingularSystem = false);

            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In,
                                             bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT virtual ~ContField3D();


            /// This function return the boundary conditions expansion.
            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                    &GetBndCondExp();

            MULTI_REGIONS_EXPORT int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey);


        protected:

            
        private:

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);

        };
        typedef std::shared_ptr<ContField3D>      ContField3DSharedPtr;




    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3D_H
