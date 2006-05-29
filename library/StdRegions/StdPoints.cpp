///////////////////////////////////////////////////////////////////////////////
//
// File StdPoints.cpp
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
// Description: Points definition 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <StdRegions/StdPoints.h>

#include <LibUtilities/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>

namespace Nektar
{
    namespace StdRegions 
    {
	
	Points::Points(const int npts):
	    m_pointsorder(npts)
	{
	    m_zeros       = new double [npts];
	    m_weights     = new double [npts];
	    m_derivmatrix = new double [npts*npts];
	}
	
	Points::~Points()
	{
	    if(m_zeros)       delete[] m_zeros;
	    if(m_weights)     delete[] m_weights;
	    if(m_derivmatrix) delete[] m_derivmatrix;
	}
	
	
	PolyPoints::PolyPoints(const int npts, double *z):
	    Points(npts)
	{
	    Blas::Dcopy(npts,z,1,m_zeros,1);
	    CalcWeights();
	    CalcDerivMatrix();
	}
	
	
	void PolyPoints::CalcWeights()
	{
	    ErrorUtil::Error(ErrorUtil::efatal,"PolyPoints::CalcWeights",
			     "Routine needs defining");
	    // need to set up Lagrange routine
	}
	
	
	void  PolyPoints::CalcDerivMatrix()
	{
	    ErrorUtil::Error(ErrorUtil::efatal,"PolyPoints::CalcDerivMatrix",
			     "Routine needs defining");
	    
	    // need to set up Lagrange routine
	}
	
	GaussPolyPoints::GaussPolyPoints(const int npts, PointsType ptype, 
				       const double alpha, const double beta):
	    PolyPoints(npts),
	    m_gausspointstype(ptype),
	    m_alpha(alpha),
	    m_beta(beta)
	{
	    
	    switch(m_gausspointstype)
	    {
	    case eGauss:
		Polylib::zwgj(m_zeros,m_weights,m_pointsorder,m_alpha,
			      m_beta);
		Polylib::Dgj(m_derivmatrix,m_zeros,m_pointsorder,m_alpha,
			     m_beta);
		break;
	    case eLobatto:
		Polylib::zwglj(m_zeros,m_weights,m_pointsorder,m_alpha,
			       m_beta);
		
		Polylib::Dglj(m_derivmatrix,m_zeros,m_pointsorder,m_alpha,
			      m_beta);
		break;
	    case eRadauM:
		Polylib::zwgrjm(m_zeros,m_weights,m_pointsorder,m_alpha,
				m_beta);
		Polylib::Dgrjm(m_derivmatrix,m_zeros,m_pointsorder,m_alpha,
			       m_beta);
		break;
	    case eRadauP:
		Polylib::zwgrjp(m_zeros,m_weights,m_pointsorder,m_alpha,
				m_beta);
		Polylib::Dgrjp(m_derivmatrix,m_zeros,m_pointsorder,m_alpha,
			       m_beta);
		break;
	    default:
		ErrorUtil::Error(ErrorUtil::efatal,
				 "GaussPolyPoints::GaussPolyPoints",
				 "Unknown points type");
		break;
	    }      
	}
	
	
	FourierPoints::FourierPoints(const int npts):
	    Points(npts)
	{
	    
	    
	    if(m_pointsorder%2)
	    {
		ErrorUtil::Error(ErrorUtil::efatal,
				 "FourierPoints::FourierPoints",
				 "Fourier points need to be of even order");
	    }
	    
	    // define points in the region [-1,1]
	    for(int i = 0; i < m_pointsorder; ++i)
	    {
		m_zeros  [i] = -1.0 + i*2.0/(double)m_pointsorder;
		m_weights[i] =  2.0/(double)m_pointsorder;
	    }      
	    
	    ErrorUtil::Error(ErrorUtil::ewarning,
			     "FourierPoints::FourierPoints",
			     "Fourier points Derivative matrix "
			     "need defining for Fourier spacing");
	}     
	
    } // end of namespace stdregion
} // end of namespace stdregion

/** 
 * $Log: StdPoints.cpp,v $
 * Revision 1.4  2006/05/29 16:02:00  sherwin
 * Used default spacing
 *
 * Revision 1.3  2006/05/29 13:43:29  sherwin
 * U
 *
 * Revision 1.2  2006/05/29 13:40:40  sherwin
 * Tab checking
 *
 * Revision 1.1  2006/05/08 05:41:52  sherwin
 * Initial attempt at sorting out points structure
 *
 **/ 




