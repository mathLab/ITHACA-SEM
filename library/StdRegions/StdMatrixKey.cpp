///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrixKey.cpp
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
// Description: Definition of StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#include "StdRegions/StdExpansion.h"
#include "StdRegions/StdMatrixKey.h"

namespace Nektar
{
    namespace StdRegions
    {
    
        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const StdExpansion &stdExpansion,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(stdExpansion.GetBase()),
            m_ncoeffs(stdExpansion.GetNcoeffs()),
            m_matrixType(matrixType),
            m_nconstants(0),
            m_nodalPointsType(nodalType)
        {
        }

        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const StdExpansion &stdExpansion,
                                   const NekDouble const0,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(stdExpansion.GetBase()),
            m_ncoeffs(stdExpansion.GetNcoeffs()),
            m_matrixType(matrixType),
            m_nconstants(1),
            m_nodalPointsType(nodalType)
        {
            m_constant = Array<OneD,NekDouble>(m_nconstants);
            m_constant[0] = const0;
        }


        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const StdExpansion &stdExpansion,
                                   const NekDouble const0,
                                   const NekDouble const1,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(stdExpansion.GetBase()),
            m_ncoeffs(stdExpansion.GetNcoeffs()),
            m_matrixType(matrixType),
            m_nconstants(2),
            m_nodalPointsType(nodalType)
        {
            m_constant = Array<OneD,NekDouble>(m_nconstants);
            m_constant[0] = const0;
            m_constant[1] = const1;
        }


        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                                   const int ncoeffs,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(base),
            m_ncoeffs(ncoeffs),
            m_matrixType(matrixType),
            m_nconstants(0),
            m_nodalPointsType(nodalType)
        {
        }


        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                                   const int ncoeffs,
                                   const NekDouble const0,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(base),
            m_ncoeffs(ncoeffs),
            m_matrixType(matrixType),
            m_nconstants(1),
            m_nodalPointsType(nodalType)
        {
            m_constant = Array<OneD,NekDouble>(m_nconstants);
            m_constant[0] = const0;
        }

        StdMatrixKey::StdMatrixKey(const MatrixType matrixType, 
                                   const ShapeType shapeType,
                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                                   const int ncoeffs,
                                   const NekDouble const0,
                                   const NekDouble const1,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(base),
            m_ncoeffs(ncoeffs),
            m_matrixType(matrixType),
            m_nconstants(2),
            m_nodalPointsType(nodalType)
        {
            m_constant = Array<OneD,NekDouble>(m_nconstants);
            m_constant[0] = const0;
            m_constant[1] = const1;
        }
        
        StdMatrixKey::StdMatrixKey(const StdMatrixKey& rhs) :
            m_shapeType(rhs.m_shapeType),
            m_base(rhs.m_base),
            m_ncoeffs(rhs.m_ncoeffs),
            m_matrixType(rhs.m_matrixType),
            m_nconstants(rhs.m_nconstants),
            m_constant(rhs.m_constant),
            m_nodalPointsType(rhs.m_nodalPointsType)
        {
        }
        
        bool StdMatrixKey::opLess::operator()(const StdMatrixKey &lhs, 
                          const StdMatrixKey &rhs) const
        {        
            return (lhs.m_matrixType < rhs.m_matrixType);
        }

        bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs)
        {   
            if(lhs.m_matrixType < rhs.m_matrixType)
            {
                return true;
            }
            
            if(lhs.m_matrixType > rhs.m_matrixType)
            {
                return false;
            }
            
            if(lhs.m_ncoeffs < rhs.m_ncoeffs)
            {
                return true;
            }
            
            if(lhs.m_ncoeffs > rhs.m_ncoeffs)
            {
                return false;
            }
            
            for(unsigned int i = 0; i < ShapeTypeDimMap[lhs.m_shapeType]; ++i)
            {
                if(lhs.m_base[i].get() < rhs.m_base[i].get())
                {
                    return true;
                }
                
                if(lhs.m_base[i].get() > rhs.m_base[i].get())
                {
                    return false;
                }
            }

            if(lhs.m_nconstants != rhs.m_nconstants)
            {
                return false;
            }
            else {

                for(unsigned int i = 0; i < lhs.m_nconstants; ++i)
                {
                    if(lhs.m_constant[i] < rhs.m_constant[i])
                    {
                        return true;
                    }
                    
                    if(lhs.m_constant[i] > rhs.m_constant[i])
                    {
                        return false;
                    }
                }
            }
            
            if(lhs.m_nodalPointsType < rhs.m_nodalPointsType)
            {
                return true;
            }
            
            if(lhs.m_nodalPointsType > rhs.m_nodalPointsType)
            {
                return false;
            }
            
            return false;
        }
        
        std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs)
        {
            os << "MatrixType: " << MatrixTypeMap[rhs.GetMatrixType()] << ", ShapeType: " 
                << ShapeTypeMap[rhs.GetShapeType()] << ", Ncoeffs: " << rhs.GetNcoeffs() 
                << std::endl;

            if(rhs.GetNconstants())
            {
                os << "Constants: " << endl;
                for(unsigned int i = 0; i < rhs.GetNconstants(); ++i)
                {
                    os << "\t value " << i <<" : " << rhs.GetConstant(i) << endl;
                }
            }
            
            for(unsigned int i = 0; i < ShapeTypeDimMap[rhs.GetShapeType()]; ++i)
            {
                os << rhs.GetBase()[i]->GetBasisKey();
            }

            return os;
        }
    }
}

/**
* $Log: StdMatrixKey.cpp,v $
* Revision 1.11  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.10  2007/08/11 23:42:26  sherwin
* A few changes
*
* Revision 1.9  2007/07/26 02:39:21  bnelson
* Fixed Visual C++ compiler errors when compiling in release mode.
*
* Revision 1.8  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.7  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.6  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.5  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.4  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.3  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.2  2007/03/05 08:07:13  sherwin
* Modified so that StdMatrixKey has const calling arguments in its constructor.
*
* Revision 1.1  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.6  2007/02/28 09:53:17  sherwin
* Update including adding GetBasis call to StdExpansion
*
* Revision 1.5  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.4  2007/02/23 19:26:04  jfrazier
* General bug fix and formatting.
*
* Revision 1.3  2007/02/22 22:02:27  sherwin
* Update with executing StdMatManager
*
* Revision 1.2  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.1  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
***/

