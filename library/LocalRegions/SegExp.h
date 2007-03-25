///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.h
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
// Description: Header file for SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SEGEXP_H
#define SEGEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdSegExp.h>
#include <SpatialDomains/SegGeom.h>

#include <SpatialDomains/GeomFactors.h>

#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/LinSys.hpp>

#include <fstream>

namespace Nektar
{
    namespace LocalRegions 
    {
    
	class SegExp: public StdRegions::StdSegExp
        {

        public:

            /// \brief Constructor using BasisKey class for quadrature
            /// points and order definition 
            SegExp(const LibUtilities::BasisKey &Ba, 
		   SpatialDomains::SegGeomSharedPtr &geom);

            ///Copy Constructor
            SegExp(const SegExp &S);

            ///Destructor
            ~SegExp();

            /// Return Shape of region, using  ShapeType enum list. i.e. Segment  
            StdRegions::ShapeType DetShapeType() 
            { 
                return StdRegions::eSegment;
            }    

            void GetCoords(NekDoubleSharedArray &coords_1,
			   NekDoubleSharedArray &coords_2 = NullNekDoubleSharedArray,
			   NekDoubleSharedArray &coords_3 = NullNekDoubleSharedArray);
            void GetCoord(const NekDoubleSharedArray &Lcoords, 
			  NekDoubleSharedArray &coords);


            SpatialDomains::SegGeomSharedPtr GetGeom()
            {
                return m_geom;
            }

            void WriteToFile(FILE *outfile);
            void WriteToFile(std::ofstream &outfile, const int dumpVar);


            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            double Integral(const NekDoubleSharedArray &inarray);


            /// \brief  Inner product of \a inarray over region with respect to the 
            /// expansion basis (this)->_Base[0] and return in \a outarray 
            void IProductWRTBase(const NekDoubleSharedArray &inarray, 
				 NekDoubleSharedArray &outarray);


            //-----------------------------
            // Differentiation Methods
            //-----------------------------


            /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
            physical quadrature points given by \a inarray and return in \a
            outarray. */
            void PhysDeriv(const NekDoubleSharedArray &inarray, 
			   NekDoubleSharedArray &out_d0 = NullNekDoubleSharedArray,
			   NekDoubleSharedArray &out_d1 = NullNekDoubleSharedArray,
			   NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray);

            //----------------------------
            // Evaluations Methods
            //---------------------------

            /** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->_coeffs  */
            void FwdTrans(const NekDoubleSharedArray &inarray, 
			  NekDoubleSharedArray &outarray);

            double PhysEvaluate(const NekDoubleSharedArray &coord);

        protected:

            void GenMetricInfo();    

            DNekMatSharedPtr    CreateMatrix(const MatrixKey &mkey);
            DNekLinSysSharedPtr CreateLinSys(const LinSysKey &mkey);

            SpatialDomains::SegGeomSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekMat, MatrixKey::opLess> m_matrixManager;

            LibUtilities::NekManager<LinSysKey, DNekLinSys, LinSysKey::opLess> m_linSysManager;

            /// \brief  Inner product of \a inarray over region with respect to
            /// the expansion basis \a base and return in \a outarray 
            inline void IProductWRTBase(const NekDouble *base, 
					const NekDoubleSharedArray &inarray, 
					NekDoubleSharedArray &outarray, 
					const int coll_check);

        private:

            virtual StdRegions::ShapeType v_DetShapeType() 
            {
                return DetShapeType();
            }

            virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
            {
                return m_metricinfo;
            }

            virtual void v_GetCoords(NekDoubleSharedArray &coords_0,
				     NekDoubleSharedArray &coords_1 = NullNekDoubleSharedArray,
				     NekDoubleSharedArray &coords_2 = NullNekDoubleSharedArray)
	    {
                GetCoords(coords_0, coords_1, coords_2);
            }

            virtual void v_GetCoord(const NekDoubleSharedArray &lcoord, 
				     NekDoubleSharedArray &coord)
            {
                GetCoord(lcoord, coord);
            }

            virtual int v_GetCoordim()
            {
                return m_geom->GetCoordim();
            }

            virtual void v_WriteToFile(FILE *outfile)
            {
                WriteToFile(outfile);
            }

            virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                WriteToFile(outfile,dumpVar);
            }

            virtual SpatialDomains::GeomType v_MetricInfoType()
            {
                return m_metricinfo->GetGtype();
            }

            /// \brief Virtual call to integrate the physical point list \a inarray
            /// over region (see SegExp::Integral) 
            virtual NekDouble v_Integral(const NekDoubleSharedArray &inarray )
            {
                return Integral(inarray);
            }

            /** \brief Virtual call to SegExp::IProduct_WRT_B */
            virtual void v_IProductWRTBase(const NekDoubleSharedArray &inarray,
					   NekDoubleSharedArray &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            /// Virtual call to SegExp::PhysDeriv
            virtual void v_StdPhysDeriv(const NekDoubleSharedArray & inarray, 
					NekDoubleSharedArray &outarray)
            {
                StdSegExp::PhysDeriv(inarray, outarray);
            }

            virtual void v_PhysDeriv(const NekDoubleSharedArray &inarray, 
				     NekDoubleSharedArray &out_d0 = NullNekDoubleSharedArray,
				     NekDoubleSharedArray &out_d1 = NullNekDoubleSharedArray,
				     NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray)
	    {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            /// Virtual call to SegExp::FwdTrans
            virtual void v_FwdTrans(const NekDoubleSharedArray &inarray, 
				    NekDoubleSharedArray &outarray)
            {
                FwdTrans(inarray,outarray);
            }

	    /** \brief Virtual call to SegExp::FwdTrans */
	    virtual void v_FwdTrans(const StdExpansion1D &in)
	    {
		FwdTrans(((SegExp &) in).GetPhys(), m_coeffs);
	    }
      

            /// Virtual call to SegExp::Evaluate
            virtual double v_PhysEvaluate(const NekDoubleSharedArray &coords)
            {
                return PhysEvaluate(coords);
            }

            /** \brief Virtual function to evaluate the discrete \f$ L_\infty\f$
            error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
            u_{exact}\f$ is given by the array \a sol. 

            The full function is defined in StdExpansion::Linf 

            Input: 

            - \a _phys: takes the physical value space array as
            approximate solution

            - \a sol: array of solution function  at physical quadrature points

            output: 

            - returns the \f$ L_\infty \f$ error as a double. 
            */
            virtual double v_Linf(const NekDoubleSharedArray &sol)
            {
                return Linf(sol);
            }

            /** \brief Virtual function to evaluate the \f$ L_\infty \f$ norm of
            the function defined at the physical points \a (this)->_phys. 

            The full function is defined in StdExpansion::Linf 

            Input: 

            - \a _phys: uses the physical value space array as discrete
            function to be evaulated.

            output: 

            - returns the \f$ L_\infty \f$  as a double. 
            */
            virtual double v_Linf()
            {
                return Linf();
            }

            /** \brief Virtual function to evaluate the \f$ L_2\f$, \f$ |
            \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2 dx
            \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by the
            array sol.

            The full function is defined in StdExpansion::L2 

            Input: 

            - \a _phys: takes the physical value space array as
            approximate solution
            - \a sol: array of solution function  at physical quadrature points

            output: 

            - returns the \f$ L_2 \f$ error as a double. 
            */
            virtual double v_L2(const NekDoubleSharedArray &sol)
            {
                return StdExpansion::L2(sol);
            }

            /** \brief Virtual function to evaluate the \f$ L_2\f$ norm of the
            function defined at the physical points \a (this)->_phys.  

            The full function is defined in StdExpansion::L2 

            Input: 

            - \a _phys: uses the physical value space array as discrete
            function to be evaulated.

            output: 

            - returns the \f$ L_2 \f$  as a double. 
            */
            virtual double v_L2()
            {
                return StdExpansion::L2();
            }
        };

        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<SegExp> SegExpSharedPtr;
        typedef std::vector< SegExpSharedPtr > SegExpVector;
        typedef std::vector< SegExpSharedPtr >::iterator SegExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // SEGEXP_H

//
// $Log: SegExp.h,v $
// Revision 1.8  2007/03/20 09:13:38  kirby
// new geomfactor routines; update for metricinfo; update style
//
// Revision 1.7  2007/03/14 21:24:07  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.6  2007/03/02 12:01:55  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.5  2007/01/15 21:12:26  sherwin
// First definition
//
// Revision 1.4  2006/06/13 18:05:01  sherwin
// Modifications to make MultiRegions demo ProjectLoc2D execute properly.
//
// Revision 1.3  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.2  2006/05/29 17:05:49  sherwin
// Modified to put shared_ptr around geom definitions
//
// Revision 1.1  2006/05/04 18:58:46  kirby
// *** empty log message ***
//
// Revision 1.33  2006/03/13 18:20:33  sherwin
//
// Fixed error in ResetGmat
//
// Revision 1.32  2006/03/12 21:59:48  sherwin
//
// compiling version of LocalRegions
//
//
