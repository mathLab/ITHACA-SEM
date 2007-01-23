///////////////////////////////////////////////////////////////////////////////
//
// File StdSegExp.h
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
// Description: Header file for Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H
#define NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpMap.h>

namespace Nektar
{
    namespace StdRegions
    {
      
	class StdSegExp: public StdExpansion1D
        {
	public:
	    /** \brief Default constructor */
	    StdSegExp();

	    /** \brief Constructor using BasisKey class for quadrature points and 
	     *  order definition
	     *
	     *  \param Ba BasisKey class definition containing order and quadrature 
	     *  points.
	     */
	    StdSegExp(const LibUtilities::BasisKey &Ba);
	    
	    /** \brief Constructor using BasisKey class for quadrature
	     *  points and order definition where m_coeffs and m_phys are all set.
	     *
	     *  \param Ba BasisKey definition containing order and quadrature points.
	     *  \param coeffs list of expansions coefficient to be set in m_coeffs
	     *  \param phys list of physical values to be set in m_phys
	     */
	    StdSegExp(const LibUtilities::BasisKey &Ba, double *coeffs, double *phys);

	    /** \brief Copy Constructor */
	    StdSegExp(const StdSegExp &T);
	    
	    /** \brief Destructor */
	    ~StdSegExp();
	    
	    /** \brief Return Shape of region, using  ShapeType enum list.
	     *  i.e. Segment
	     */
	    ShapeType DetShapeType()
	    {
		return eSegment;
	    };

	    //----------------------------
	    // Integration Methods
	    //----------------------------
	    
	    /** \brief Integrate the physical point list \a inarray over region 
	     *  and return the value
	     *
	     *  \param inarray definition of function to be integrated evauluated at 
	     *  quadrature point of expansion. 
	     *  \return returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
	     *  = u(\xi_{1i}) \f$
	     */
	    double Integral(const double *inarray);
	    
	    
	    /** \brief Inner product of \a inarray over region with respect to the
	     *  expansion basis (this)->m_base[0] and return in \a outarray
	     *
	     *  Wrapper call to StdSegExp::IProduct_WRT_B
	     *  \param inarray array of function values evaluated at the physical
	     *  collocation points
	     *  \param outarray  the values of the inner product with respect to 
	     *  each basis over region will be stored in the array \a outarray as
	     *  output of the function
	     */
	    void IProductWRTBase(const double * inarray, double * outarray);
	    
	    void FillMode(const int mode, double *outarray);
	    
	    //----------------------------------
	    // Local Matrix Routines
	    //----------------------------------
	    
	    /** \brief Generate local mass matrix \f$ {\bf M}[i][j] =
	     *  \int^1_{-1} \phi_i(\xi_1) \phi_j(\xi_1) d\xi_1 \f$ in standard
	     *  region and store in \a outarray  
	     *
	     *  \param outarray As output of the function, the local mass matrix is 
	     *  stored in the array \a outarray. The matrix is in row major format and 
	     *  is stored as \a outarray[j*(this)->m_ncoeffs + i]
	     */
	    DNekMatSharedPtr GenMassMatrix();
	    
	    /** \brief Generate local weak Laplacian matrix \f$ {\bf L}[i][j] =
	     *  \int^1_{-1} \frac{d \phi_i(\xi_1)}{d \xi_1}\frac{d
	     *  \phi_j(\xi_1)}{d \xi_1} d\xi_1 \f$ in standard region and store in
	     *  \a outarray 
	     *
	     *  \param outarray As output of the function, the local Laplacian matrix  
	     *  is stored in the array \a outarray. The matrix is in row major format 
	     *  and is stored as \a outarray[j*(this)->m_ncoeffs + i]
	     */
	    DNekMatSharedPtr GenLapMatrix();

	    //-----------------------------
	    // Differentiation Methods
	    //-----------------------------
	    
	    /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the physical 
	     *  quadrature opoints in the expansion (i.e. (this)->m_phys) and return 
	     *  in \a outarray.
	     *
	     *  This is a wrapper function around StdExpansion1D::TensorDeriv
	     *
	     *  This function takes the physical value space array \a m_phys as 
	     *  discrete function to be evaluated
	     *  \param outarray the resulting array of the derivative \f$
	     *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array \a outarray
	     */
	    
	    void PhysDeriv(double * outarray);
	    
	    /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the physical 
	     *  quadrature points given by \a inarray and return in \a outarray.
	     *
	     *  This is a wrapper around StdExpansion1D::Tensor_Deriv
	     *  \param inarray array of a function evaluated at the quadrature points
	     *  \param  outarray the resulting array of the derivative \f$
	     *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array \a outarra 
	     */
	    void PhysDeriv(const double *inarray, double * outarray);
	    
	    
	    //----------------------------
	    // Evaluations Methods
	    //---------------------------
	    
	    /** \brief Backward transform from coefficient space (stored in \a
	     *  (this)->m_coeffs) and evaluate at the physical quadrature points 
	     *  \a  outarray 
	     *
	     *  Operation can be evaluated as \f$ u(\xi_{1i}) =
	     *  \sum_{p=0}^{order-1} \hat{u}_p \phi_p(\xi_{1i}) \f$ or equivalently 
	     *  \f$ {\bf u} = {\bf B}^T {\bf \hat{u}} \f$ where
	     *  \f${\bf B}[i][j] = \phi_i(\xi_{1j}), \mbox{\_coeffs}[p] = {\bf
	     *  \hat{u}}[p] \f$
	     *
	     *  The function takes the coefficient space array \a m_coeffs as input
	     *  for the transformation
	     *
	     *  \param outarray the resulting array of the values of the function at 
	     *  the physical quadrature points will be stored in the array \a outarray
	     */
	    void BwdTrans(double * outarray);
	    
	    /** \brief Forward transform from physical quadrature space stored in 
	     *  \a inarray and evaluate the expansion coefficients and store in 
	     *  \a (this)->m_coeffs
	     *
	     *  Perform a forward transform using a Galerkin projection by taking the 
	     *  inner product of the physical points and multiplying by the inverse of
	     *  the mass matrix using the Solve method of the standard matrix 
	     *  container holding the local mass matrix, i.e. \f$ {\bf \hat{u}} = 
	     *  {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =  \int^1_{-1} 
	     *  \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$
	     *
	     *  This function stores the expansion coefficients calculated by the 
	     *  transformation in the coefficient space array \a m_coeffs
	     *
	     *  \param inarray array of physical quadrature points to be transformed
	     */ 
	    void FwdTrans(const double * inarray);
	    
	    /** \brief Single Point Evaluation: \f$ u(x) = \sum_p \phi_p(x) \hat{u}_p 
	     *  = \sum_p h_p(x) u(x_p)\f$
	     */
	    double Evaluate(const double *Lcoords);
	    
	    void MapTo(EdgeOrientation dir, StdExpMap& Map);
	    
	protected:
	    
	    /** \brief  Inner product of \a inarray over region with respect to 
	     *  expansion basis \a base and return in \a outarray 
	     *
	     *  Calculate \f$ I[p] = \int^{1}_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1
	     *  = \sum_{i=0}^{nq-1} \phi_p(\xi_{1i}) u(\xi_{1i}) w_i \f$ where
	     *  \f$ outarray[p] = I[p], inarray[i] = u(\xi_{1i}), base[p*nq+i] =
	     *  \phi_p(\xi_{1i}) \f$.
	     *
	     *  \param  base an array defining the local basis for the inner product 
	     *  usually passed from Basis->GetBdata() or Basis->GetDbdata()
	     *  \param inarray: physical point array of function to be integrated
	     *  \f$ u(\xi_1) \f$
	     *  \param coll_check flag to identify when a Basis->Collocation() call 
	     *  should be performed to see if this is a GLL_Lagrange basis with a 
	     *  collocation property. (should be set to 0 if taking the inner product 
	     *  with respect to the derivative of basis)
	     *  \param outarray  the values of the inner product with respect to 
	     *  each basis over region will be stored in the array \a outarray as
	     *  output of the function
	     */
	    void IProductWRTBase(const double *base, const double *inarray,
				 double *outarray, int coll_check);
	    
	private:

	    virtual int v_GetNverts()
	    {
		return 2;
	    } 
	    
	    virtual ShapeType v_DetShapeType()
		{
		    return DetShapeType();
		};
	    
	    /** \brief Virtual call to integrate the physical point list \a inarray 
	     *  over region (see StdSegExp::Integral)
	     */
	    virtual double v_Integral(const double *inarray )
	    {
		return Integral(inarray);
	    } 
	    
	    /** \brief Virtual call to StdSegExp::IProduct_WRT_B */
	    virtual void v_IProductWRTBase(const double * inarray, double * outarray)
	    {
		IProductWRTBase(inarray,outarray);
	    } 
	    
	    virtual void v_FillMode(const int mode, double *outarray)
	    {
		FillMode(mode,outarray);
	    } 
 
	    /** \brief Virtual call to GenMassMatrix */
	    virtual DNekMatSharedPtr v_GenMassMatrix()
	    {
		return GenMassMatrix();
	    }
      
	    virtual DNekMatSharedPtr v_GenLapMatrix()
	    {
		return GenLapMatrix();
	    }

	    /** \brief Virtual call to StdSegExp::Deriv */
	    virtual void v_PhysDeriv(double * outarray)
	    {
		PhysDeriv(&m_phys[0], outarray);
	    }
	    
	    /** \brief Virtual call to StdSegExp::Deriv */
	    virtual void v_StdPhysDeriv(double * outarray)
	    {
		PhysDeriv(&m_phys[0], outarray);
	    }
      
	    /** \brief Virtual call to StdSegExp::Deriv */
	    virtual void v_PhysDeriv(const double *inarray, double * outarray)
	    {
		PhysDeriv(inarray, outarray);
	    }
      
	    /** \brief Virtual call to StdSegExp::Deriv */
	    virtual void v_StdPhysDeriv(const double *inarray, double * outarray)
            {
		PhysDeriv(inarray, outarray);
	    }
	    
	    /** \brief Virtual call to StdSegExp::BwdTrans */
	    virtual void v_BwdTrans(double * outarray)
	    {
		BwdTrans(outarray);
	    }
      
	    /** \brief Virtual call to StdSegExp::FwdTrans */
	    virtual void v_FwdTrans(const double * inarray)
	    {
		FwdTrans(inarray);
	    }
      
	    /** \brief Virtual call to StdSegExp::Evaluate */
	    virtual double v_Evaluate(const double * Lcoords)
	    {
		return Evaluate(Lcoords);
	    }
      
	    virtual void v_MapTo(EdgeOrientation dir, StdExpMap &Map)
	    {
		MapTo(dir,Map);
	    }
	};
	
    } //end of namespace
} //end of namespace

#endif //STDSEGEXP_H

/**
 * $Log: StdSegExp.h,v $
 * Revision 1.5  2007/01/15 15:07:26  pvos
 * updating doxygen documentation
 *
 * Revision 1.4  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.3  2006/07/02 17:16:19  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.2  2006/06/01 14:13:37  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.41  2006/03/13 18:29:35  sherwin
 *
 * Corrected error with definition of GetCoords
 *
 * Revision 1.40  2006/03/05 22:11:03  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.39  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.38  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.37  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/
