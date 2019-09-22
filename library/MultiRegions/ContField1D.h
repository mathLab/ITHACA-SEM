///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.h
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
// Description: Field definition in one-dimension
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H


#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContField.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSys.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Abstraction of a global continuous one-dimensional spectral/hp
        /// element expansion which approximates the solution of a set of
        /// partial differential equations.
        class ContField1D: public ContField
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ContField1D();

            /// Set up global continuous field based on an input mesh and
            /// boundary conditions.
            MULTI_REGIONS_EXPORT ContField1D(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph1D,
                        const std::string &variable  = "DefaultVar",
                        const bool DeclareCoeffPhysArrays = true,
                        const bool CheckIfSingularSystem = false,
                        const Collections::ImplementationType ImpType
                                  = Collections::eNoImpType);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D
                (const ContField1D &In,
                 bool DeclareCoeffPhysArrays = true);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D
                (const LibUtilities::SessionReaderSharedPtr &pSession,
                 const ExpList & In);

            /// Destructor
            MULTI_REGIONS_EXPORT virtual ~ContField1D();



        protected:

        private:

#if 0
            /// Impose the Dirichlet Boundary Conditions on outarray 
            MULTI_REGIONS_EXPORT virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);
#endif

        };
        typedef std::shared_ptr<ContField1D>      ContField1DSharedPtr;


    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTSOLNFIELD1D_H
