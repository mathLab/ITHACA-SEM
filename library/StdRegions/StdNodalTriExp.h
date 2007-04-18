///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.h
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
// Description: Header field for Nodal triangle routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDNODALTRIEXP_H
#define STDNODALTRIEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdTriExp.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/LocalRegionsDeclarations.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        class StdNodalTriExp: public StdTriExp
        {

        public:

            //Constructors
            StdNodalTriExp(const LibUtilities::BasisKey &Ba, 
			   const LibUtilities::BasisKey &Bb,
			   const LibUtilities::PointsType Ntype);
	    
            //Copy Constructor
            StdNodalTriExp(const StdNodalTriExp &T);

            //Destructor
            ~StdNodalTriExp();

            ShapeType DetShapeType()
            {
                return eTriangle;
            }

            /////////////////////////////
            /// Nodal basis specific routines
            ///////////////////////////

            DNekMatSharedPtr GenNBasisTransMatrix();
	    

            void NodalToModal();

            void NodalToModal(NekDoubleSharedArray &in_out_array);

            void NodalToModalTranspose();
            void NodalToModalTranspose(NekDoubleSharedArray &in_out_array);

            void ModalToNodal();
            void ModalToNodal(NekDoubleSharedArray &in_out_array);

            void GetNodalPoints(ConstNekDoubleSharedArray &x, 
				ConstNekDoubleSharedArray &y)
	    {
		LibUtilities::PointsManager()[*m_nodalPointsKey]->GetPoints(x,y);
            }

            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

            void IProductWRTBase(ConstNekDoubleSharedArray inarray, 
				 NekDoubleSharedArray &outarray);

	    /** \brief Fill outarray with nodal mode \a mode of expansion
	     *   and put in m_phys
	     */
            void FillMode(const int mode, NekDoubleSharedArray &outarray);

	    void  MapTo(const int edge_ncoeffs,
			const LibUtilities::BasisType Btype, 
			const int eid, 
			const EdgeOrientation eorient,
			StdExpMap &Map);

	    void  MapTo_ModalFormat(const int edge_ncoeffs, 
				    const LibUtilities::BasisType Btype, 
				    const int eid, 
				    const EdgeOrientation eorient,
				    StdExpMap &Map);

            //-----------------------------
            // Evaluations Methods
            //-----------------------------

	    void BwdTrans(ConstNekDoubleSharedArray inarray,
			  NekDoubleSharedArray &outarray);
	    void FwdTrans(ConstNekDoubleSharedArray inarray,
			  NekDoubleSharedArray &outarray);
	    
        protected:

	    /** \brief Calculate the inner product of inarray with respect to
	     *  the basis B=base0[p]*base1[pq] and put into outarray
	     *
	     *  This function uses the StdTriExp routine and then 
	     *  calls ModalToNodal to transform to Nodal basis
	     */
	    void IProductWRTBase(ConstNekDoubleSharedArray base0, 
				 ConstNekDoubleSharedArray base1,
				 ConstNekDoubleSharedArray inarray, 
				 NekDoubleSharedArray &outarray);

        private:

            boost::shared_ptr<LibUtilities::PointsKey> m_nodalPointsKey;

            virtual ShapeType v_DetShapeType()
            {
                return DetShapeType();
            }

	    virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i)
	    {
		return GetEdgeBasisType(i);
	    }

            virtual DNekMatSharedPtr v_GenNBasisTransMatrix() 
            {
		return GenNBasisTransMatrix();
	    }


            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

	    virtual NekDouble v_Integral(ConstNekDoubleSharedArray inarray)
	    {
		return Integral(inarray);
	    }

	    virtual void v_IProductWRTBase(ConstNekDoubleSharedArray inarray,
					   NekDoubleSharedArray &outarray)
	    {
		IProductWRTBase(inarray, outarray);
	    }

	    virtual void v_FillMode(const int mode, NekDoubleSharedArray &outarray)
	    {
		FillMode(mode, outarray);
	    }
            //-----------------------------
            // Differentiation Methods
            //-----------------------------
	    virtual void v_PhysDeriv(ConstNekDoubleSharedArray inarray,
				     NekDoubleSharedArray &out_d0,
				     NekDoubleSharedArray &out_d1 = NullNekDoubleSharedArray,
				     NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray)
            {
		PhysDeriv(inarray, out_d1, out_d2);
	    }

            //-----------------------------
            // Evaluations Methods
            //-----------------------------

	    virtual void v_BwdTrans(ConstNekDoubleSharedArray inarray,
				    NekDoubleSharedArray &outarray)
            {
		BwdTrans(inarray,outarray);
	    }

	    /** \brief Virtual call to StdTriExp::FwdTrans */
	    virtual void v_FwdTrans(ConstNekDoubleSharedArray inarray,
				    NekDoubleSharedArray &outarray)
            {
		FwdTrans(inarray,outarray);
	    }

	    virtual NekDouble v_PhysEvaluate(ConstNekDoubleSharedArray coords)
            {
		return StdTriExp::PhysEvaluate(coords);
	    }
	    
	    
	    virtual void v_WriteToFile(std::ofstream &outfile)
	    {
		WriteToFile(outfile);
	    }
	    
	    virtual void v_MapTo(const int edge_ncoeffs, 
				 const LibUtilities::BasisType Btype, 
				 const int eid, 
				 const EdgeOrientation eorient,
				 StdExpMap &Map)
	    {
		MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
	    }

	    virtual void v_MapTo_ModalFormat(const int edge_ncoeffs, 
					     const LibUtilities::BasisType Btype, 
					     const int eid, 
					     const EdgeOrientation eorient,
					     StdExpMap &Map)
	    {
		MapTo_ModalFormat(edge_ncoeffs,Btype,eid,eorient,Map);
	    }
	    

	    virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
	    {
		WriteCoeffsToFile(outfile);
	    }
        };

    } //end of namespace
} //end of namespace

#endif //STDNODALTRIEXP_H

/**
* $Log: StdNodalTriExp.h,v $
* Revision 1.7  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.6  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.5  2007/01/15 21:13:46  sherwin
* Nodal stuff correction and added Field Classes
*
* Revision 1.4  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.3  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.2  2006/05/23 15:08:19  jfrazier
* Minor cleanup to correct compile warnings.
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.16  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.15  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.14  2006/03/06 17:12:46  sherwin
*
* Updated to properly execute all current StdRegions Demos.
*
* Revision 1.13  2006/03/04 20:26:55  bnelson
* Added comments after #endif.
*
* Revision 1.12  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
**/

