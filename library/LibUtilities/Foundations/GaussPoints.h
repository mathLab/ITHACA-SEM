///////////////////////////////////////////////////////////////////////////////
//
// File GaussPoints.h
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
// Description: Header file of GaussPoints Distributions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        class GaussPoints: public Points<NekDouble>
        {
        public:
            virtual ~GaussPoints()
            {
            }            

            LIB_UTILITIES_EXPORT static boost::shared_ptr< Points<NekDouble> > Create(const PointsKey &pkey);

            LIB_UTILITIES_EXPORT boost::shared_ptr< NekMatrix<NekDouble> > CreateMatrix(const PointsKey &pkey);

            LIB_UTILITIES_EXPORT const boost::shared_ptr<NekMatrix<NekDouble> > GetI(const PointsKey &pkey);
            LIB_UTILITIES_EXPORT const boost::shared_ptr<NekMatrix<NekDouble> > GetI(const Array<OneD, const NekDouble>& x);
            LIB_UTILITIES_EXPORT const boost::shared_ptr<NekMatrix<NekDouble> > GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x);

            LIB_UTILITIES_EXPORT boost::shared_ptr< NekMatrix<NekDouble> > CreateGPMatrix(const PointsKey &pkey);
            
            LIB_UTILITIES_EXPORT const boost::shared_ptr<NekMatrix<NekDouble> > GetGalerkinProjection(const PointsKey &pkey);
            


            GaussPoints(const PointsKey &pkey):PointsBaseType(pkey)
            {
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussGaussLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauPLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussLobattoLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussGaussChebyshev),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMChebyshev),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauPChebyshev),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussLobattoChebyshev),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha1Beta0),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha2Beta0),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussKronrodLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauKronrodMLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauKronrodMAlpha1Beta0),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eGaussLobattoKronrodLegendre),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, eFourierEvenlySpaced),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
                m_InterpManager.RegisterCreator(PointsKey(0, ePolyEvenlySpaced),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));
		m_InterpManager.RegisterCreator(PointsKey(0, eBoundaryLayerPoints),
                    boost::bind(&GaussPoints::CreateMatrix, this, _1));            
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussGaussLegendre),
                   boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauPLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussLobattoLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussGaussChebyshev),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMChebyshev),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauPChebyshev),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussLobattoChebyshev),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha1Beta0),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha2Beta0),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussKronrodLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauKronrodMLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussRadauKronrodMAlpha1Beta0),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eGaussLobattoKronrodLegendre),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eFourierEvenlySpaced),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
                m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, ePolyEvenlySpaced),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));
		m_GalerkinProjectionManager.RegisterCreator(PointsKey(0, eBoundaryLayerPoints),
                    boost::bind(&GaussPoints::CreateGPMatrix, this, _1));            }


        private:
            /// These should not be called.  All creation is done
            /// using the constructor requiring the key, declared
            /// above.
            GaussPoints();
            GaussPoints(const GaussPoints &points);

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp);


            boost::shared_ptr<NekMatrix<NekDouble> > CalculateGalerkinProjectionMatrix(const PointsKey &pkey);


	    /// functions used by the Kronrod points
	    NekDouble LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts, const Array<OneD, const NekDouble>& funcvals);
            NekDouble LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts);     
            NekDouble LagrangePolyDeriv(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts);

        }; // class GaussPoints
    } // end of namespace
} // end of namespace 

#endif //GAUSSPOINTS_H
