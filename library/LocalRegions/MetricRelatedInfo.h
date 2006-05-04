///////////////////////////////////////////////////////////////////////////////
//
// File MetricRelatedInfo.h
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
// Description: Metric information related to LocalRegions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef METRICINFO_H
#define METRICINFO_H

#include <SpatialDomains/GeoFac.h>
#include <StdRegions/StdMatrix.h>

namespace Nektar
{
  namespace LocalRegions 
  {
  
    class SegExp;

    class MetricRelatedInfo: public SpatialDomains::GeoFac
    {
    
    public:

      MetricRelatedInfo(const StdRegions::GeomType gtype, const int expdim, 
			const int coordim):  
	GeoFac(gtype,expdim,coordim) 
      { 
      } 

      /** \brief One dimensional geometric factors based on one, two or three
	  dimensional coordinate description 
      **/
      MetricRelatedInfo(const StdRegions::GeomType gtype, const int coordim,   
  			const StdRegions::StdExpansion1D **Coords):  
  	GeoFac(gtype,coordim,Coords)  
      {  
      }  

      ~MetricRelatedInfo() 
      { 
      } 
    
       StdRegions::StdMatContainer *GetMassMatrix(SegExp *S);
       
    protected: 
    
    private: 
       StdRegions::StdMatrix m_elmtmats; 
       
    }; 
      
  } //end of namespace
} //end of namespace

#endif // METRICINFO_H

/** 
 *    $Log: MetricRelatedInfo.h,v $
 *    Revision 1.5  2006/03/13 19:47:54  sherwin
 *
 *    Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
 *
 *    Revision 1.4  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.3  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
