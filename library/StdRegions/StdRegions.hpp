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


// Headers from LibUtilities needed in StdRegions
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

#include <string>
#include <loki/Factory.h>


namespace Nektar
{
    /** \brief The namespace associated with the the StdRegions library 
     * (\ref pageStdRegions "StdRegions introduction")
     */
    namespace StdRegions
    {
        /** \page pageStdRegions The StdRegions library
         *
         */

        enum ElementType
        {
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
            "HybridDGLamToU"
        };

        /** enum list of StdExpansion regions */
        enum ExpansionType
        {
            eNoExpansionType,
            eSegment,
            eTriangle,
            eQuadrilateral,
            eTetrahedron,
            ePyramid,
            ePrism,
            eHexahedron,
            SIZE_ExpansionType
        };
    

        const char* const ExpansionTypeMap[] = 
        {
            "NoExpansionType",
            "Segment",
            "Triangle",
            "Quadrilateral",
            "Tetrahedron",
            "Pyramid",
            "Prism",
            "Hexahedron"
        };
    



        // Hold the dimension of each of the types of shapes.
        const unsigned int ExpansionTypeDimMap[SIZE_ExpansionType] =
        {
            0,  // Unknown
            1,  // eSegment
            2,  // eTriangle
            2,  // eQuadrilateral
            3,  // eTetrahedron
            3,  // ePyramid
            3,  // ePrism
            3,  // eHexahedron
        };

        enum EdgeOrientation
        {
            eForwards,
            eBackwards
        };

        const char* const EdgeOrientationMap[] = 
        {
            "Forwards",
            "Backwards"
        };

        enum FaceOrientation
        {
            eDir1FwdDir1_Dir2FwdDir2,
            eDir1FwdDir1_Dir2BwdDir2,
            eDir1BwdDir1_Dir2FwdDir2,
            eDir1BwdDir1_Dir2BwdDir2,
            eDir1FwdDir2_Dir2FwdDir1,
            eDir1FwdDir2_Dir2BwdDir1,
            eDir1BwdDir2_Dir2FwdDir1,
            eDir1BwdDir2_Dir2BwdDir1
        };

        const char* const FaceOrientationMap[] = 
        {
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

/**
* $Log: StdRegions.hpp,v $
* Revision 1.37  2009/11/09 15:41:38  sehunchun
* eWeakDirectionalDerivative and eMassLevelCurvature are added.
*
* Revision 1.36  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.35  2009/04/03 14:57:34  sherwin
* Linear Advection matrices added, corrected unsigned int intialisation
*
* Revision 1.34  2009/02/16 16:06:20  pvos
* Update of TimeIntegration classes
*
* Revision 1.33  2008/12/18 14:11:35  pvos
* NekConstants Update
*
* Revision 1.32  2008/11/19 16:02:47  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.31  2008/11/12 12:12:10  pvos
* Time Integration update
*
* Revision 1.30  2008/11/05 16:08:15  pvos
* Added elemental optimisation functionality
*
* Revision 1.29  2008/09/12 11:26:39  pvos
* Updates for mappings in 3D
*
* Revision 1.28  2008/08/14 22:09:51  sherwin
* Modifications to remove HDG routines from StdRegions and removed StdExpMap
*
* Revision 1.27  2008/07/19 21:12:54  sherwin
* Removed MapTo function and made orientation convention anticlockwise in UDG routines
*
* Revision 1.26  2008/07/10 13:03:49  pvos
* Added eUnifiedDGHelmBndSysForce to the enum list
*
* Revision 1.25  2008/06/05 15:06:06  pvos
* Added documentation
*
* Revision 1.24  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.23  2008/03/12 15:25:09  pvos
* Clean up of the code
*
* Revision 1.22  2008/01/23 09:09:25  sherwin
* New matrix listing for DG stuff
*
* Revision 1.21  2007/11/20 16:29:48  sherwin
* Added enum matrixtype definitions for UDG solver
*
* Revision 1.20  2007/07/20 02:16:55  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.19  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.18  2007/07/10 19:27:58  kirby
* Update for new matrix structures
*
* Revision 1.17  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.16  2007/05/15 05:18:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.15  2007/04/18 16:09:13  pvos
* Added some new Tensor Operations routines
*
* Revision 1.14  2007/04/03 03:56:13  bnelson
* Moved Lapack.hpp, Blas.hpp, Transf77.hpp to LinearAlgebra
*
* Revision 1.13  2007/03/21 20:56:43  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.12  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.11  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.10  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.9  2007/02/23 19:26:08  jfrazier
* General bug fix and formatting.
*
* Revision 1.8  2007/02/22 22:02:28  sherwin
* Update with executing StdMatManager
*
* Revision 1.7  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.6  2007/02/17 03:40:20  jfrazier
* Couple changes to reflect additions and corrections to reflect linear algebra calls.
*
* Revision 1.5  2007/01/28 18:34:24  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.4  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.3  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.2  2006/06/01 13:43:20  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:33  kirby
* *** empty log message ***
*
* Revision 1.8  2006/03/23 21:55:02  jfrazier
* Declared NekConstants namespace to replace NekConstants class.
*
* Revision 1.7  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.6  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.5  2006/03/05 22:11:03  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.4  2006/03/04 20:26:55  bnelson
* Added comments after #endif.
*
* Revision 1.3  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.2  2006/02/19 13:26:13  sherwin
*
* Coding standard revisions so that libraries compile
*
* Revision 1.1  2006/02/15 08:06:36  sherwin
*
* Put files into coding standard (although they do not compile)
*
*
**/
