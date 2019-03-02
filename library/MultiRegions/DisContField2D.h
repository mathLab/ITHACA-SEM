///////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.h
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
// Description: Field definition in two-dimensions for a discontinuous LDG-H
// expansion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H

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
        class DisContField2D : public DisContField
        {
        public:
            MULTI_REGIONS_EXPORT DisContField2D();

            MULTI_REGIONS_EXPORT DisContField2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr   &graph2D,
                const std::string                          &variable,
                const bool SetUpJustDG            = true,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D                     &In,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const std::string                        &variable,
                const bool SetUpJustDG            = false,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D &In,
                const bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT virtual ~DisContField2D();
            

        protected:

            /**< Bases needed for the expansion */
            Array<OneD, LibUtilities::BasisSharedPtr> m_base; 

            /** \brief This function gets the shared point to basis
             *
             *  \return returns the shared pointer to the bases
             */
            inline const Array<OneD, const LibUtilities::BasisSharedPtr>&
                GetBase() const
            {
                return(m_base);
            }

            /** \brief This function returns the type of basis used in
             *  the \a dir direction
             *
             *  The different types of bases implemented in the code
             *  are defined in the LibUtilities::BasisType enumeration
             *  list. As a result, the function will return one of the
             *  types of this enumeration list.
             *
             *  \param dir the direction \return returns the type of
             *  basis used in the \a dir direction
             */
            inline  LibUtilities::BasisType GetBasisType(const int dir) const
            {
                ASSERTL1(dir < m_base.num_elements(), "dir is larger than m_numbases");
                return(m_base[dir]->GetBasisType());
            }


            Array<OneD, Array<OneD, unsigned int> > m_mapEdgeToElmn;
            Array<OneD, Array<OneD, unsigned int> > m_signEdgeToElmn;
            Array<OneD,StdRegions::Orientation>     m_edgedir;

#if 0 
            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);
#endif


        };
        
        typedef std::shared_ptr<DisContField2D>   DisContField2DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
