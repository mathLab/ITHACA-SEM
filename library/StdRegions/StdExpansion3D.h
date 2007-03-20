///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion3D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP3D_H
#define STDEXP3D_H

#include <StdRegions/StdExpansion.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion3D: public StdExpansion
        {

        public:
            StdExpansion3D();
            StdExpansion3D(int numcoeffs, const LibUtilities::BasisKey &Ba, 
			    const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc);
            StdExpansion3D(const StdExpansion3D &T);
            ~StdExpansion3D();

            // Differentiation

	    /** \brief Calculate the 3D derivative in the local 
	     *  tensor/collapsed coordinate at the physical points 
	     *	
	     *	This function is independent of the expansion basis and can
	     *	therefore be defined for all tensor product distribution of
	     *	quadrature points in a generic manner.  The key operations are:
	     *
	     *	- \f$ \frac{d}{d\eta_1} \rightarrow {\bf D^T_0 u } \f$ \n
	     *	- \f$ \frac{d}{d\eta_2} \rightarrow {\bf D_1 u } \f$
	     *	- \f$ \frac{d}{d\eta_3} \rightarrow {\bf D_2 u } \f$
	     *
	     *  \param inarray array of physical points to be differentiated
	     *  \param  outarray_d1 the resulting array of derivative in the 
	     *  \f$\eta_1\f$ direction will be stored in outarray_d1 as output
	     *  of the function
	     *  \param outarray_d2 the resulting array of derivative in the 
	     *  \f$\eta_2\f$ direction will be stored in outarray_d2 as output 
	     *  of the function
	     *  \param outarray_d3 the resulting array of derivative in the 
	     *  \f$\eta_3\f$ direction will be stored in outarray_d3 as output 
	     *  of the function
	     */
            void PhysTensorDeriv(const NekDoubleSharedArray &inarray, NekDoubleSharedArray &outarray_d1,
                NekDoubleSharedArray &outarray_d2, NekDoubleSharedArray &outarray_d3);

            /** \brief Evaluate a function at points coords which is assumed
	     *  to be in local collapsed coordinate format. The function is
	     *  assumed to be in physical space
	     */
            double PhysEvaluate(const NekDoubleSharedArray &coords);

            void PhysDeriv(const NekDoubleSharedArray &inarray, 
			   NekDoubleSharedArray &outarray_d0, 
			   NekDoubleSharedArray &outarray_d1, 
			   NekDoubleSharedArray &outarray_d2)
            {
                v_PhysDeriv(inarray, outarray_d0, outarray_d1, outarray_d2);
            }

            void StdPhysDeriv(const NekDoubleSharedArray &inarray,
			  NekDoubleSharedArray &outarray_d0,
			  NekDoubleSharedArray &outarray_d1,  
			  NekDoubleSharedArray &outarray_d2)
            {
                v_StdPhysDeriv(inarray, outarray_d0, outarray_d1, outarray_d2);
            }


        protected:

        private:

            // Virtual Functions ----------------------------------------

	    virtual int v_GetNverts() = 0;
	    virtual int v_GetNedges() = 0;
	    virtual int v_GetNfaces() = 0;

            virtual void v_GenMassMatrix(NekDoubleSharedArray & outarray)   = 0;
            virtual void v_GenLapMatrix (NekDoubleSharedArray & outarray)   = 0;
            virtual ShapeType v_DetShapeType()                = 0;

            virtual DNekMatSharedPtr v_GetMassMatrix()        = 0;
            virtual DNekMatSharedPtr v_GetLapMatrix()         = 0;

            virtual int v_get_nodalpoints(const NekDoubleSharedArray &x, const NekDoubleSharedArray &y)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return 0;
            }

            virtual void v_GenNBasisTransMatrix(NekDoubleSharedArray &outarray)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
            }

            virtual void   v_BwdTrans (const NekDoubleSharedArray &inarray, 
				       NekDoubleSharedArray &outarray)      = 0;
            virtual void   v_FwdTrans (const NekDoubleSharedArray &inarray,
				       NekDoubleSharedArray &outarray)      = 0;

            virtual NekDouble v_Integral(const NekDoubleSharedArray &inarray ) = 0;
            virtual NekDouble v_Evaluate(const NekDoubleSharedArray &coords) = 0;

            virtual void   v_PhysDeriv(const NekDoubleSharedArray &inarray, 
				       NekDoubleSharedArray &outarray_d0,
				       NekDoubleSharedArray &outarray_d1, 
				       NekDoubleSharedArray &outarray_d2) = 0;
            virtual void   v_StdPhysDeriv(const NekDoubleSharedArray &inarray, 
					  NekDoubleSharedArray &outarray_d0,
					  NekDoubleSharedArray &outarray_d1,
					  NekDoubleSharedArray &outarray_d2) = 0;

	    virtual int v_GetCoordim(void)
	    {
                return 3; 
	    }

        };

    } //end of namespace
} //end of namespace

#endif //STDEXP3D_H

/**
* $Log: StdExpansion3D.h,v $
* Revision 1.4  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.3  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.2  2006/06/13 18:05:02  sherwin
* Modifications to make MultiRegions demo ProjectLoc2D execute properly.
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.12  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.11  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.10  2006/03/06 17:12:45  sherwin
*
* Updated to properly execute all current StdRegions Demos.
*
* Revision 1.9  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.8  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.7  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
**/

