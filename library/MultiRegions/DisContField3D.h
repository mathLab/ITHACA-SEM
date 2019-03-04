///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3D.h
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
// Description: Field definition in three-dimensions for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {        
        class AssemblyMapDG;
        

        class DisContField3D : public DisContField
        {
        public:
            MULTI_REGIONS_EXPORT DisContField3D();

            MULTI_REGIONS_EXPORT DisContField3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr   &graph3D,
                const std::string                          &variable,
                const bool                                  SetUpJustDG = true,
                const bool                                  DeclareCoeffPhysArrays = true, 
                const Collections::ImplementationType ImpType
                = Collections::eNoImpType);

            MULTI_REGIONS_EXPORT DisContField3D(
                const DisContField3D                       &In,
                const SpatialDomains::MeshGraphSharedPtr   &graph3D,
                const std::string                          &variable,
                const bool                                 SetUpJustDG = false,
                const bool                                 DeclareCoeffPhysArrays = true);
            
            /// Constructs a global discontinuous field based on another
            /// discontinuous field.
            MULTI_REGIONS_EXPORT DisContField3D(const DisContField3D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField3D();


            MULTI_REGIONS_EXPORT void EvaluateHDGPostProcessing(
                Array<OneD, NekDouble> &outarray);

        protected:
            
            virtual void v_AddFwdBwdTraceIntegral(
                const Array<OneD, const NekDouble> &Fwd, 
                const Array<OneD, const NekDouble> &Bwd, 
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);

        };

        typedef std::shared_ptr<DisContField3D> DisContField3DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3D_H
