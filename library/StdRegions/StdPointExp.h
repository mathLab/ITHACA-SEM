///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/PointExp.h,v $
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
// Description: Definition of a Point expansion 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDPOINTEXP_H
#define NEKTAR_LIBS_STDREGIONS_STDPOINTEXP_H

#include <StdRegions/StdExpansion0D.h>
#include <StdRegions/StdRegionsDeclspec.h>


namespace Nektar
{
    namespace StdRegions
    {

		class StdPointExp: virtual public StdExpansion0D
    {
        
    public:
		
		/// Class representing a point  in reference space
		
        /// All interface of this class sits in StdExpansion class
		
		STD_REGIONS_EXPORT StdPointExp();
		
		STD_REGIONS_EXPORT StdPointExp(const LibUtilities::BasisKey &Ba);
		
		STD_REGIONS_EXPORT StdPointExp(const StdPointExp &T);
		
		STD_REGIONS_EXPORT ~StdPointExp();
		
    protected:
            Array<OneD, NekDouble > m_coeffs; //!< Array containing expansion coefficients
            Array<OneD, NekDouble > m_phys; //!< Array containing physical point which is likely to be the same as the coefficient but is defined for consistency (It is also used in Robin boundary conditions) 
            
            //----------------------------
            // Evaluations Methods
            //---------------------------
            
            STD_REGIONS_EXPORT virtual void v_GetCoords(
                                          Array<OneD, NekDouble> &coords_0,
                                          Array<OneD, NekDouble> &coords_1,
                                          Array<OneD, NekDouble> &coords_2);
            
            //----------------------------
            // Helper functions
            //---------------------------
            
            STD_REGIONS_EXPORT virtual LibUtilities::ShapeType
                v_DetShapeType() const;
            
            //-----------------------------
            // Transforms
            //-----------------------------
            
            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_FwdTrans_BndConstrained(
	  		      const Array<OneD, const NekDouble>& inarray,
                              Array<OneD, NekDouble> &outarray);
            
            //----------------------------
            // Inner product functions
            //----------------------------
            
            
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                              const Array<OneD, const NekDouble>& inarray,
                              Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                              const Array<OneD, const NekDouble>& base,
                              const Array<OneD, const NekDouble>& inarray,
                              Array<OneD, NekDouble> &outarray,
                              int coll_check);
		
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
		               const Array<OneD, const NekDouble>& inarray,
                               Array<OneD, NekDouble> &outarray);
		
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase (
				const int dir,
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray);

            //----------------------------
            // Evaluations Methods
            //---------------------------
            
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                                 const StdMatrixKey &mkey);
		
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                                 const StdMatrixKey &mkey);
		
            
    private:
            virtual int v_GetNverts() const
            {
                return 1;
            }
            
    };

    // type defines for use of PointExp in a boost vector
    typedef boost::shared_ptr<StdPointExp> StdPointExpSharedPtr;

    } //end of namespace
} //end of namespace

#endif // STDPOINTEXP_H
