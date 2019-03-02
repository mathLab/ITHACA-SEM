//////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.cpp
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
// Description: Field definition for 2D domain with boundary conditions using
// LDG flux.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h> 
#include <LocalRegions/QuadExp.h>   
#include <LocalRegions/TriExp.h>    
#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField2D
         * Abstraction of a global discontinuous two-dimensional spectral/hp
         * element expansion which approximates the solution of a set of
         * partial differential equations.
         */
        
        /**
         * @brief Default constructor.
         */
        DisContField2D::DisContField2D(void)
            : DisContField            ()
        {
        }

        DisContField2D::DisContField2D(
            const DisContField2D &In, 
            const bool            DeclareCoeffPhysArrays)
            : DisContField         (In,DeclareCoeffPhysArrays)
        {
        }

        /**
         * @brief Constructs a global discontinuous field based on an input
         * mesh with boundary conditions.
         */
        DisContField2D::DisContField2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph2D,
            const std::string                          &variable,
            const bool                                  SetUpJustDG,
            const bool                                  DeclareCoeffPhysArrays,
            const Collections::ImplementationType       ImpType)
            : DisContField(pSession, graph2D, variable, SetUpJustDG,
                           DeclareCoeffPhysArrays,ImpType)
        {
        }


        /*
         * @brief Copy type constructor which declares new boundary conditions
         * and re-uses mapping info and trace space if possible
         */
        DisContField2D::DisContField2D(
            const DisContField2D                     &In,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const std::string                        &variable,
            const bool                                SetUpJustDG,
            const bool                                DeclareCoeffPhysArrays):
            DisContField(In,graph2D, variable, SetUpJustDG, DeclareCoeffPhysArrays)
        {
        }

        /**
         *
         */
        DisContField2D::~DisContField2D()
        {
        }
        



#if 0 
        
        /**
         * @brief Calculates the result of the multiplication of a global matrix
         * of type specified by @a mkey with a vector given by @a inarray.
         * 
         * @param mkey      Key representing desired matrix multiplication.
         * @param inarray   Input vector.
         * @param outarray  Resulting multiplication.
         */
        void DisContField2D::v_GeneralMatrixOp(
               const GlobalMatrixKey             &gkey,
               const Array<OneD,const NekDouble> &inarray,
               Array<OneD,      NekDouble> &outarray)
        {
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);
            const DNekScalBlkMatSharedPtr& HDGHelm = GetBlockMatrix(gkey);

            m_traceMap->GlobalToLocalBnd(inarray, loc_lambda);
            LocLambda = (*HDGHelm) * LocLambda;
            m_traceMap->AssembleBnd(loc_lambda,outarray);
        }
#endif
        
    } // end of namespace
} //end of namespace
