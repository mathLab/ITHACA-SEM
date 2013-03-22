///////////////////////////////////////////////////////////////////////////////
//
// File MatrixKey.cpp
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
// Description: Definition of MatrixKey based on StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////


#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
    namespace LocalRegions
    {
        MatrixKey::MatrixKey(const StdRegions::MatrixType matrixType,
                  const LibUtilities::ShapeType shapeType,
                  const StdRegions::StdExpansion &stdExpansion,
                  const StdRegions::ConstFactorMap &factorMap,
                  const StdRegions::VarCoeffMap &varCoeffMap,
                  LibUtilities::PointsType nodalType) :
            StdMatrixKey(matrixType, shapeType, stdExpansion, factorMap, varCoeffMap, nodalType),
            m_metricinfo(stdExpansion.GetMetricInfo())
        {
        }

        MatrixKey::MatrixKey(const MatrixKey& mkey,
                      const StdRegions::MatrixType matrixType) :
            StdRegions::StdMatrixKey(mkey, matrixType),
            m_metricinfo(mkey.m_metricinfo)
        {
        }

        MatrixKey::MatrixKey(const StdRegions::StdMatrixKey &mkey) :
            StdRegions::StdMatrixKey(mkey)
        {
        }

        bool MatrixKey::opLess::operator()(const MatrixKey &lhs, const MatrixKey &rhs) const
        {        
            {
                return (lhs.GetMatrixType() < rhs.GetMatrixType());
            }
        }

        bool operator<(const MatrixKey &lhs, const MatrixKey &rhs)
        {
            if(lhs.m_metricinfo.get() < rhs.m_metricinfo.get())
            {
                return true;
            }


            if(lhs.m_metricinfo.get() > rhs.m_metricinfo.get())
            {
                return false;
            }    
            
            return (*dynamic_cast<const StdRegions::StdMatrixKey*>(&lhs)
                    < *dynamic_cast<const StdRegions::StdMatrixKey*>(&rhs));
        }

    }
}

/**
* $Log: MatrixKey.cpp,v $
* Revision 1.23  2008/12/17 16:57:04  pvos
* Performance updates
*
* Revision 1.22  2008/12/17 12:29:36  pvos
* Minor update
*
* Revision 1.21  2008/12/16 14:09:41  pvos
* Performance updates
*
* Revision 1.20  2008/11/19 16:46:43  pvos
* Fixed minor bug
*
* Revision 1.19  2008/11/19 16:01:41  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.18  2008/07/09 11:39:47  sherwin
* Removed m_scalefactor and made operator< dependent upon StdMatKey
*
* Revision 1.17  2008/06/02 23:33:46  ehan
* Fixed warning : no new line at end of file
*
* Revision 1.16  2008/05/30 00:33:48  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.15  2008/05/29 01:02:13  bnelson
* Added precompiled header support.
*
* Revision 1.14  2007/12/17 13:04:29  sherwin
* Modified GenMatrix to take a StdMatrixKey and removed m_constant from MatrixKey
*
* Revision 1.13  2007/11/20 16:28:45  sherwin
* Added terms for UDG Helmholtz solver
*
* Revision 1.12  2007/07/28 05:09:32  sherwin
* Fixed version with updated MemoryManager
*
* Revision 1.11  2007/07/26 02:39:21  bnelson
* Fixed Visual C++ compiler errors when compiling in release mode.
*
* Revision 1.10  2007/07/20 00:45:50  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/07/12 12:52:58  sherwin
* Updated to have a helmholtz matrix
*
* Revision 1.8  2007/07/10 17:17:22  sherwin
* Introduced Scaled Matrices into the MatrixManager
*
* Revision 1.7  2007/05/27 16:10:28  bnelson
* Update to new Array type.
*
* Revision 1.6  2007/04/08 03:33:30  jfrazier
* Minor reformatting and fixing SharedArray usage.
*
* Revision 1.5  2007/04/04 21:49:23  sherwin
* Update for SharedArray
*
* Revision 1.4  2007/03/20 09:13:37  kirby
* new geomfactor routines; update for metricinfo; update style
*
* Revision 1.3  2007/03/09 20:41:50  jfrazier
* *** empty log message ***
*
* Revision 1.2  2007/03/05 08:06:07  sherwin
* Updated to use MemoryManager
*
* Revision 1.1  2007/03/04 20:42:14  sherwin
* Keys for matrix managers
*
***/
