///////////////////////////////////////////////////////////////////////////////
//
// File StdPoints.h
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
// Description: Header file of Points definition 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDPOINTS_H
#define STDPOINTS_H

#include <math.h>

#include <StdRegions/PolyManager.h>
#include <loki/Factory.h>

#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/ErrorUtil.hpp>

namespace Nektar
{
  namespace StdRegions 
  {

    typedef boost::shared_ptr<NekMatrix> SharedNekMatrixPtr;

    struct PointsKey
    {
        int m_pointsorder; //< Number of points
        PointsType m_pointstype;
        double m_params[2];
    };

    // base class for points definition whic can store the zeros and
    // must have a definition of weight and a derivative matrix
    
    class Points
    {
    public:
        Points():
            m_pointsorder(0),
            m_zeros(NULL),
	        m_weights(NULL),
	        m_derivmatrix(NULL)
        {
            }

      Points(const int npts);

      ~Points();

      inline double *GetZ() const
      {
	return m_zeros; 
      }

      inline double *GetW() const
      {
	return m_weights; 
      }
      
      inline void GetZW(const double *&z, const double *&w) 
      {
	z = m_zeros;
	w = m_weights;
      }

      inline double *GetD() const
      {
	return m_derivmatrix;
      }

    protected:
      double *m_zeros;       //< definition of points
      double *m_weights;     //< definition of integration weights
    
        SharedNekMatrixPtr
      double *m_derivmatrix; //< defintiion of derivative matrix
      
    private:

    };


    
    // definition of points using given points through which we assume
    // a polynomial fit
    class PolyPoints: public Points
    {

    public:
      PolyPoints()
      {
      }

      PolyPoints(const int npts):
	Points(npts)
      {
      }

      PolyPoints(const int npts, double *z);

      void CalcWeights();
      
      void CalcDerivMatrix();

    protected:

    private:

    };


    // definition of points using given points through which we assume
    // a polynomial fit
    class GaussPolyPoints: public PolyPoints
    {
    public:
      GaussPolyPoints()
      {
      }

      GaussPolyPoints(const int npts, PointsType ptype, 
		      const double alpha, const double beta);
      
    protected:

    private:
      
    };


    // definition of points using discrete Fourier fit through
    // equispaced points
    class FourierPoints: public Points
    {
    public:

    protected:

    private:
      FourierPoints()
      {
      }

      FourierPoints(const int npts);
    };

  } // end of namespace
} // end of namespace 

#endif //STDPOINTS_H
