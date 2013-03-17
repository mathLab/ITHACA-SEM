///////////////////////////////////////////////////////////////////////////////
//
// File MultiRegsions.hpp
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
// Description: Multiregion overall header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_H
#define MULTIREGIONS_H

#include <vector>
#include <SpatialDomains/Conditions.h>


namespace Nektar
{
    namespace MultiRegions
    {

        enum CoeffState
        {
            eNotSet,   ///< No state set
            eLocal,    ///< Local  coefficients 
            eHybrid,   ///< Hybrid coefficients
            eGlobal,   ///< Global coefficients 
        };

        // Orientation of adjacent edge for use with boundary
        // constraints
        enum AdjacentTraceOrientation
        {
            eAdjacentEdgeIsForwards,
            eAdjacentEdgeIsBackwards
        };

        // Orientation of adjacent face for use with boundary
        // constraints
        enum AdjacentFaceOrientation
        {
            eAdjacentFaceDir1FwdDir1_Dir2FwdDir2,
            eAdjacentFaceDir1FwdDir1_Dir2BwdDir2,
            eAdjacentFaceDir1BwdDir1_Dir2FwdDir2,
            eAdjacentFaceDir1BwdDir1_Dir2BwdDir2,
            eAdjacentFaceDir1FwdDir2_Dir2FwdDir1,
            eAdjacentFaceDir1FwdDir2_Dir2BwdDir1,
            eAdjacentFaceDir1BwdDir2_Dir2FwdDir1,
            eAdjacentFaceDir1BwdDir2_Dir2BwdDir1
        };

        enum GlobalSysSolnType
        {
            eNoSolnType,    ///< No Solution type specified
            eDirectFullMatrix,
            eDirectStaticCond,
            eDirectMultiLevelStaticCond,
            eIterativeFull,
            eIterativeStaticCond,
            eIterativeMultiLevelStaticCond,
            eXxtFullMatrix,
            eXxtStaticCond,
            eXxtMultiLevelStaticCond,
            eSIZE_GlobalSysSolnType
        };


        const char* const GlobalSysSolnTypeMap[] =
            {
            "No Solution Type",
            "DirectFull",
            "DirectStaticCond",
            "DirectMultiLevelStaticCond",
            "IterativeFull",
            "IterativeStaticCond",
            "IterativeMultiLevelStaticCond",
            "XxtFull",
            "XxtStaticCond",
            "XxtMultiLevelStaticCond"
        };

        /// Type of Galerkin projection.
        enum ProjectionType
        {
            eGalerkin,
            eDiscontinuous,
            eMixed_CG_Discontinuous
        };

        enum PreconditionerType
        {
            eNull,    ///< No Solution type specified
            eDiagonal,
            eInverseLinear,
	    eLowEnergy,
            eLinearLowEnergy,
            eBlock,
            eLocalLowEnergy
        };

        const char* const PreconditionerTypeMap[] =
        {
            "Null",
            "Diagonal",
            "InverseLinear",
	    "LowEnergy",
	    "LinearLowEnergy",
            "Block",
            "LocalLowEnergy"
        };


        typedef std::vector<SpatialDomains::BoundaryConditionType>  BndTypesVector;
        typedef std::vector<SpatialDomains::BoundaryConditionType>::iterator BndTypesVectorIter;


        // structure to hold information about robin boundary conditions

        struct RobinBCInfo
        {
            RobinBCInfo(const int id, const Array<OneD, const NekDouble > &primCoeffs):
                m_robinID(id),
                m_robinPrimitiveCoeffs(primCoeffs)
            {
            }

            virtual ~RobinBCInfo()
            {};

            int m_robinID; /// id of which edge/face is robin condition
            Array< OneD, const NekDouble > m_robinPrimitiveCoeffs;
            boost::shared_ptr<RobinBCInfo> next;
        };

        typedef boost::shared_ptr<RobinBCInfo> RobinBCInfoSharedPtr;


    }// end of namespace
}// end of namespace

#endif
