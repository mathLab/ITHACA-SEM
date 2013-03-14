///////////////////////////////////////////////////////////////////////////////
//
// File StdRegions.hpp
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
// Description: Definition of enum lists and constants
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDREGIONS_HPP
#define STDREGIONS_HPP


#include<map>
#include <LibUtilities/BasicUtils/SharedArray.hpp>  // for Array, etc
#include <LibUtilities/BasicUtils/ShapeType.hpp> 

namespace Nektar
{

    /** \brief The namespace associated with the the StdRegions library
     * (\ref pageStdRegions "StdRegions introduction")
     */
    namespace StdRegions
    {
        enum ElementType
        {
			//eStdPointExp,
            eStdSegExp,
            eSegExp,
            eStdQuadExp,
            eStdTriExp,
            eStdNodalTriExp,
            eQuadExp,
            eTriExp,
            eNodalTriExp,
            eStdHexExp,
            eStdPrismExp,
            eStdPyrExp,
            eStdTetExp,
            eStdNodalTetExp,
            eHexExp,
            ePrismExp,
            ePyrExp,
            eTetExp,
            eNodalTetExp,
            SIZE_ElementType
        };

        const char* const ElementTypeMap[] =
        {
			//"StdPointExp",
            "StdSegExp",
            "SegExp",
            "StdQuadExp",
            "StdTriExp",
            "StdNodalTriExp",
            "QuadExp",
            "TriExp",
            "NodalTriExp",
            "StdHexExp",
            "StdPrismExp",
            "StdPyrExp",
            "StdTetExp",
            "StdNodalTetExp",
            "HexExp",
            "PrismExp",
            "PyrExp",
            "TetExp",
            "NodalTetExp",
        };

        /** @todo we need to tidy up matrix construction approach
         *  probably using a factory type approach                 */
        enum MatrixType
        {
            eMass,
            eInvMass,
            eLaplacian,
            eLaplacian00,
            eLaplacian01,
            eLaplacian02,
            eLaplacian10,
            eLaplacian11,
            eLaplacian12,
            eLaplacian20,
            eLaplacian21,
            eLaplacian22,
            eInvLaplacianWithUnityMean,
            eWeakDeriv0,
            eWeakDeriv1,
            eWeakDeriv2,
            eWeakDirectionalDeriv,
            eMassLevelCurvature,
            eLinearAdvectionReaction,
            eLinearAdvectionDiffusionReaction,
            eNBasisTrans,
            eInvNBasisTrans,
            eBwdTrans,
            eIProductWRTBase,
            eIProductWRTDerivBase0,
            eIProductWRTDerivBase1,
            eIProductWRTDerivBase2,
            eHelmholtz,
            eHybridDGHelmholtz,
            eInvHybridDGHelmholtz,
            eHybridDGHelmBndLam,
            eHybridDGLamToQ0,
            eHybridDGLamToQ1,
            eHybridDGLamToQ2,
            eHybridDGLamToU,
            eFwdTrans,
            ePreconditioner,
            SIZE_MatrixType
        };

        const char* const MatrixTypeMap[] =
        {
            "Mass",
            "InvMass",
            "Laplacian",
            "Laplacian00",
            "Laplacian01",
            "Laplacian02",
            "Laplacian10",
            "Laplacian11",
            "Laplacian12",
            "Laplacian20",
            "Laplacian21",
            "Laplacian22",
            "InvLaplacianWithUnityMean",
            "WeakDeriv0",
            "WeakDeriv1",
            "WeakDeriv2",
            "WeakDirectionalDeriv",
            "MassLevelCurvature",
            "LinearAdvectionReaction",
            "LinearAdvectionDiffusionReaction",
            "NBasisTrans",
            "InvNBasisTrans",
            "BwdTrans",
            "IProductWRTBase",
            "IProductWRTDerivBase0",
            "IProductWRTDerivBase1",
            "IProductWRTDerivBase2",
            "Helmholtz",
            "HybridDGHelmholz",
            "InvHybridDGHelmholtz",
            "HybridDGHelmBndLam",
            "HybridDGLamToQ0",
            "HybridDGLamToQ1",
            "HybridDGLamToQ2",
            "HybridDGLamToU",
            "FwdTrans",
            "Preconditioner",
        };

        enum VarCoeffType
        {
            eVarCoeffMass,
            eVarCoeffLaplacian,
            eVarCoeffWeakDeriv,
            eVarCoeffD00,
            eVarCoeffD11,
            eVarCoeffD22,
            eVarCoeffD01,
            eVarCoeffD02,
            eVarCoeffD12,
            eVarCoeffVelX,
            eVarCoeffVelY
        };

        const char* const VarCoeffTypeMap[] = {
            "VarCoeffMass",
            "VarCoeffLaplacian",
            "VarCoeffWeakDeriv",
            "VarCoeffD00",
            "VarCoeffD11",
            "VarCoeffD22",
            "VarCoeffD01",
            "VarCoeffD02",
            "VarCoeffD12",
            "VarCoeffVelX",
            "VarCoeffVelY"
        };
        typedef std::map<StdRegions::VarCoeffType, Array<OneD, NekDouble> > VarCoeffMap;
        static VarCoeffMap NullVarCoeffMap;

        enum ConstFactorType
        {
            eFactorLambda,
            eFactorTau,
            eFactorTime,
            eFactorSVVCutoffRatio,
            eFactorSVVDiffCoeff
        };

        const char* const ConstFactorTypeMap[] = {
            "FactorLambda",
            "FactorTau",
            "FactorTime",
            "FactorSVVCutoffRatio",
            "FactorSVVDiffCoeff"
        };
        typedef std::map<ConstFactorType, NekDouble> ConstFactorMap;
        static ConstFactorMap NullConstFactorMap;

        enum IndexMapType
            {
                eEdgeToElement,
                eFaceToElement,
                eEdgeInterior,
                eFaceInterior,
                eBoundary,
                eVertex
		};
        
        const char* const IndexMapTypeMap[] =
            {
                "EdgeToElement",
                "FaceToElement",
                "EdgeInterior",
                "FaceInterior",
                "Boundary",
                "Vertex"
            };
	
        enum Orientation
            {
                eNoOrientation,
                eFwd,
                eBwd,
                eForwards,
                eBackwards,
                eDir1FwdDir1_Dir2FwdDir2,
                eDir1FwdDir1_Dir2BwdDir2,
                eDir1BwdDir1_Dir2FwdDir2,
                eDir1BwdDir1_Dir2BwdDir2,
                eDir1FwdDir2_Dir2FwdDir1,
                eDir1FwdDir2_Dir2BwdDir1,
                eDir1BwdDir2_Dir2FwdDir1,
                eDir1BwdDir2_Dir2BwdDir1,
                SIZE_Orientation
            };
	
        const char* const OrientationMap[] =
            {
                "NoOrientation",
                "Fwd",
                "Bwd",
                "Forwards",
                "Backwards",
                "Dir1FwdDir1_Dir2FwdDir2",
                "Dir1FwdDir1_Dir2BwdDir2",
                "Dir1BwdDir1_Dir2FwdDir2",
                "Dir1BwdDir1_Dir2BwdDir2",
                "Dir1FwdDir2_Dir2FwdDir1",
                "Dir1FwdDir2_Dir2BwdDir1",
                "Dir1BwdDir2_Dir2FwdDir1",
                "Dir1BwdDir2_Dir2BwdDir1"
            };
        
        // Defines a "fast find"
        // Assumes that first/last define the beginning/ending of
        // a continuous range of classes, and that start is
        // an iterator between first and last

        template<class InputIterator, class EqualityComparable>
        InputIterator find(InputIterator first, InputIterator last,
                           InputIterator startingpoint,
                           const EqualityComparable& value)
        {
            InputIterator val;

            if(startingpoint == first)
            {
                val = find(first,last,value);
            }
            else
            {
                val = find(startingpoint,last,value);
                if(val == last)
                {
                    val = find(first,startingpoint,value);
                    if(val == startingpoint)
                    {
                        val = last;
                    }
                }
            }
            return val;
        }

    } // end of namespace
} // end of namespace

#endif //STDREGIONS_H

