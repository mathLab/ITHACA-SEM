///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrixKeys.h
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
// Description: Headers for StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDMATRIXKEY_H
#define STDMATRIXKEY_H

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion;

        class StdMatrixKey
        {
        public:
            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const Array<OneD, NekDouble>& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          const Array<OneD, NekDouble>& varcoeffs0,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const Array<OneD,NekDouble>& constants,
                          const Array<OneD,NekDouble >& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const Array<OneD,NekDouble>& constants,
                          const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          const NekDouble const1,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          const NekDouble const1,
                          const Array<OneD, NekDouble>& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ExpansionType expansionType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          const NekDouble const0,
                          const NekDouble const1,
                          const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                          const int matrixid = NekConstants::kNekUnsetDouble,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            StdMatrixKey(const MatrixType matrixType, 
                         const ExpansionType expansionType,
                         const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                         const int ncoeffs,
                         LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);
            
            StdMatrixKey(const MatrixType matrixType, 
                         const ExpansionType expansionType,
                         const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                         const int ncoeffs,
                         const NekDouble const0,
                         LibUtilities::PointsType nodalType);


            StdMatrixKey(const MatrixType matrixType, 
                         const ExpansionType expansionType,
                         const Array<OneD, const LibUtilities::BasisSharedPtr> &base,
                         const int ncoeffs,
                         const NekDouble const0,
                         const NekDouble const1,
                         LibUtilities::PointsType nodalType);


            StdMatrixKey(const StdMatrixKey& rhs);

            virtual ~StdMatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                bool operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;
            };

            /// Used for finding value given the key in NekManager.
            friend bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            friend bool opLess::operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;

            MatrixType GetMatrixType() const
            {
                return m_matrixType;
            }

            ExpansionType GetExpansionType() const
            {
                return m_expansionType;
            }

            LibUtilities::PointsType GetNodalPointsType() const
            {
                return m_nodalPointsType;
            }
           
            int GetNcoeffs() const
            {
                return m_ncoeffs;
            }

            int GetNconstants() const
            {
                return m_nconstants;
            }

            NekDouble GetConstant(int i) const
            {
                ASSERTL1(i < m_nconstants,"requesting constant which has not been definied");
                return m_constant[i];
            }

            inline const Array<OneD,NekDouble>& GetConstants() const
            {         
                return m_constant;
            }

            inline const Array<OneD, const LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return m_base;
            }

            int GetMatrixID() const
            {
                return m_matrixid;
            }

            int GetNvariableCoefficients() const
            {
                return m_nvariablecoefficients;
            }

            int GetNvariableLaplacianCoefficients() const
            {
                return m_nvariablecoefficients;
            }

            inline const Array<OneD,NekDouble>& GetVariableCoefficient(int i) const
            {
	      // ASSERTL1(m_nvariablecoefficients > 0,"No coeffcients have been defined");   
	      if(m_nvariablecoefficients>0)
		{
		  return m_variablecoefficient[i];
		}          

	      else
		{ 
		  return NullNekDouble1DArray;  
		}
            }

            inline const Array<OneD,NekDouble>& GetVariableCoefficient() const
            {
                ASSERTL1(m_nvariablecoefficients > 0,"No coeffcients have been defined");                
                return m_variablecoefficient[0];
            }

            inline const Array<OneD,NekDouble>& GetVariableLaplacianCoefficient() const
            {
                ASSERTL1(m_nvariablecoefficients > 0,"No coeffcients have been defined");                
                return m_variablecoefficient[0];
            }

            inline const Array<OneD,NekDouble>& Get2DVariableLaplacianCoefficient(int i, int j) const
            {
                int idx = i+j;
                ASSERTL1(idx < m_nvariablecoefficients,"requesting a coefficient which has not been defined");                
                return m_variablecoefficient[idx];
            }

            inline const Array<OneD,NekDouble>& Get3DVariableLaplacianCoefficient(int i, int j) const
            {
                int idx = (i*j)+(i+j)-( (int) ((i+j)-(i+j)%3)/3 )*(i*j-1);
                ASSERTL1(idx < m_nvariablecoefficients,"requesting a coefficient which has not been defined");                
                return m_variablecoefficient[idx];
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(int dir) const
            {
                return(m_base[dir]);
            }

        protected:
            ExpansionType m_expansionType;
            Array<OneD, const LibUtilities::BasisSharedPtr> m_base;

            unsigned int m_ncoeffs;
            MatrixType   m_matrixType;
            LibUtilities::PointsType m_nodalPointsType;
            
            int                   m_nconstants;
            Array<OneD,NekDouble> m_constant;

            int m_nvariablecoefficients;
            Array<OneD, Array<OneD,NekDouble> >  m_variablecoefficient;

            int                   m_matrixid;

        private:
            StdMatrixKey();
        };

        std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs);

        typedef  boost::shared_ptr<StdMatrixKey> StdMatrixKeySharedPtr;

    } // end of namespace
} // end of namespace

#endif //STDMATRIXKEY_H

/**
* $Log: StdMatrixKey.h,v $
* Revision 1.24  2009/11/09 15:41:38  sehunchun
* eWeakDirectionalDerivative and eMassLevelCurvature are added.
*
* Revision 1.23  2009/11/07 21:09:11  sehunchun
* Add more functions with various parameters
*
* Revision 1.22  2009/07/09 21:44:49  sehunchun
* *** empty log message ***
*
* Revision 1.21  2008/11/19 16:02:47  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.20  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.19  2008/04/06 06:04:15  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.18  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.17  2007/08/11 23:42:26  sherwin
* A few changes
*
* Revision 1.16  2007/07/26 02:39:21  bnelson
* Fixed Visual C++ compiler errors when compiling in release mode.
*
* Revision 1.15  2007/07/22 23:04:26  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.14  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.13  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.12  2007/05/25 17:48:58  jfrazier
* Changed the class so the bases are not stored as an Array, but rather a ConstArray.
*
* Revision 1.11  2007/05/22 02:01:50  bnelson
* Changed Array::size to Array::num_elements.
*
* Fixed some compiler errors in assertions.
*
* Revision 1.10  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.9  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.8  2007/04/08 03:36:58  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.7  2007/04/04 21:49:25  sherwin
* Update for SharedArray
*
* Revision 1.6  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.5  2007/03/21 20:56:43  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.4  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.3  2007/03/05 08:07:14  sherwin
* Modified so that StdMatrixKey has const calling arguments in its constructor.
*
* Revision 1.2  2007/03/02 12:01:53  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.1  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
*
***/
