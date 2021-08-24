///////////////////////////////////////////////////////////////////////////////
//
// File FourierPoints.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header file of 1D Evenly-Spaced Point Definitions 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FOURIERPOINTS_H
#define FOURIERPOINTS_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>  // for NekManager


namespace Nektar
{
    namespace LibUtilities 
    {
        class FourierPoints: public Points<NekDouble>
        {
            public:
                virtual ~FourierPoints()
                {
                }            

                LIB_UTILITIES_EXPORT static std::shared_ptr< PointsBaseType > Create(const PointsKey &key);
                LIB_UTILITIES_EXPORT std::shared_ptr< NekMatrix<NekDouble> > CreateMatrix(const PointsKey &pkey);

                LIB_UTILITIES_EXPORT const MatrixSharedPtrType GetI(const PointsKey &pkey);
                LIB_UTILITIES_EXPORT const MatrixSharedPtrType GetI(const Array<OneD, const NekDouble>& x);
                LIB_UTILITIES_EXPORT const MatrixSharedPtrType GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x);

                FourierPoints(const PointsKey &key):PointsBaseType(key)
                {
                    namespace pl = std::placeholders;
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussGaussLegendre),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMLegendre),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauPLegendre),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussLobattoLegendre),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussGaussChebyshev),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMChebyshev),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauPChebyshev),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussLobattoChebyshev),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha1Beta0),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eGaussRadauMAlpha2Beta0),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, ePolyEvenlySpaced),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                    m_InterpManager.RegisterCreator(PointsKey(0, eFourierEvenlySpaced),
                        std::bind(&FourierPoints::CreateMatrix, this, pl::_1));
                }

            private:
                static bool initPointsManager[];
                /// Default constructor should not be called except by Create method.
                FourierPoints();

                /// Copy constructor should not be called.
                FourierPoints(const FourierPoints &points);

                void CalculatePoints();
                void CalculateWeights();
                void CalculateDerivMatrix();

                void CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, 
                                           Array<OneD, NekDouble>& interp);
                               NekDouble PeriodicSincFunction(const NekDouble x, const NekDouble h);
        }; // class FourierPoints
    } // end of namespace
} // end of namespace 


#endif //FOURIERPOINTS_H
