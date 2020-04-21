///////////////////////////////////////////////////////////////////////////////
//
// File Expansion1D.h
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
// Description: Header file for Expansion1D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION1D_H
#define EXPANSION1D_H

#include <LocalRegions/Expansion.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <StdRegions/StdExpansion1D.h>
#include <SpatialDomains/Geometry1D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        class Expansion2D;
        typedef std::shared_ptr<Expansion2D>  Expansion2DSharedPtr;
        typedef std::weak_ptr<Expansion2D>    Expansion2DWeakPtr;
        
        class Expansion1D;
        typedef std::shared_ptr<Expansion1D>  Expansion1DSharedPtr;
        typedef std::weak_ptr<Expansion1D>    Expansion1DWeakPtr;
        typedef std::vector< Expansion1DSharedPtr > Expansion1DVector;

        class Expansion1D: virtual public Expansion,
            virtual public StdRegions::StdExpansion1D
        {
            public:
                LOCAL_REGIONS_EXPORT Expansion1D(SpatialDomains::
                                                 Geometry1DSharedPtr pGeom)
                                                 : Expansion(pGeom),
                                                   StdExpansion1D()
                {
                    m_elementEdgeLeft  = -1;
                    m_elementEdgeRight = -1;
                }

                LOCAL_REGIONS_EXPORT virtual ~Expansion1D() {}

                LOCAL_REGIONS_EXPORT void AddNormTraceInt(
                        const int dir,
                        Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray);

                inline Expansion2DSharedPtr GetLeftAdjacentElementExp() const;

                inline Expansion2DSharedPtr GetRightAdjacentElementExp() const;

                inline int GetLeftAdjacentElementEdge() const;

                inline int GetRightAdjacentElementEdge() const;

                inline void SetAdjacentElementExp(
                    int                  edge,
                    Expansion2DSharedPtr &e);

                void AddHDGHelmholtzTraceTerms(
                    const NekDouble                      tau,
                    const Array<OneD, const NekDouble>  &inarray,
                          Array<OneD, NekDouble>        &outarray);

                inline SpatialDomains::Geometry1DSharedPtr GetGeom1D() const;

            protected:
                virtual DNekMatSharedPtr v_GenMatrix(
                    const StdRegions::StdMatrixKey      &mkey);

                virtual void v_AddRobinMassMatrix(
                    const int                            vert,
                    const Array<OneD, const NekDouble > &primCoeffs,
                    DNekMatSharedPtr                    &inoutmat);

                virtual void v_AddRobinEdgeContribution(
                    const int                            vert,
                    const Array<OneD, const NekDouble > &primCoeffs,
                    const Array<OneD, NekDouble>        &incoeffs,
                    Array<OneD, NekDouble>        &coeffs);

                virtual NekDouble v_VectorFlux(
                    const Array<OneD, Array<OneD, NekDouble> > &vec);

            private:
                Expansion2DWeakPtr m_elementLeft;
                Expansion2DWeakPtr m_elementRight;
                int                m_elementEdgeLeft;
                int                m_elementEdgeRight;

        };
        
        inline Expansion2DSharedPtr Expansion1D::
            GetLeftAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(),
                     "Left adjacent element not set.");
            return m_elementLeft.lock();
        }

        inline Expansion2DSharedPtr Expansion1D::
            GetRightAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(),
                     "Right adjacent element not set.");
            
            return m_elementRight.lock();
        }

        inline int Expansion1D::GetLeftAdjacentElementEdge() const
        {
            return m_elementEdgeLeft;
        }

        inline int Expansion1D::GetRightAdjacentElementEdge() const
        {
            return m_elementEdgeRight;
        }

        inline void Expansion1D::SetAdjacentElementExp(
            int                  edge,
            Expansion2DSharedPtr &e)
        {
            if (m_elementLeft.lock().get())
            {
                ASSERTL1(!m_elementRight.lock().get(),
                         "Both adjacent elements already set.");
                
                m_elementRight     = e;
                m_elementEdgeRight = edge;
            }
            else
            {
                m_elementLeft     = e;
                m_elementEdgeLeft = edge;
            }
        }

        inline SpatialDomains::Geometry1DSharedPtr Expansion1D
            ::GetGeom1D() const
        {
            return std::dynamic_pointer_cast<SpatialDomains
                ::Geometry1D>(m_geom);
        }
    } //end of namespace
} //end of namespace

#endif

