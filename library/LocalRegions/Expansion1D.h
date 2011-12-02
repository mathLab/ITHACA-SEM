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
// Description: Header file for Expansion1D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION1D_H
#define EXPANSION1D_H

#include <LocalRegions/Expansion.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
//#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        class Expansion2D;
        typedef boost::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        
        class Expansion1D: virtual public Expansion, virtual public StdRegions::StdExpansion1D
        {
            public:
                LOCAL_REGIONS_EXPORT void AddNormTraceInt(
                        const int dir,
                        Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray);

                inline Expansion2DSharedPtr GetLeftAdjacentElementExp() const;

                inline Expansion2DSharedPtr GetRightAdjacentElementExp() const;

                inline int GetLeftAdjacentElementEdge() const;

                inline int GetRightAdjacentElementEdge() const;

                inline void SetAdjacentElementExp(
                        int edge,
                        Expansion2DSharedPtr &e);

            protected:
                virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey);

                virtual void v_AddHDGHelmholtzTraceTerms(
                        const NekDouble tau,
                        const Array<OneD,const NekDouble> &inarray,
                              Array<OneD,NekDouble> &outarray);

                virtual void v_AddRobinMassMatrix(
                        const int vert,
                        const Array<OneD, const NekDouble > &primCoeffs,
                        DNekMatSharedPtr &inoutmat);

                virtual void v_AddRobinEdgeContribution(
                        const int vert,
                        const Array<OneD, const NekDouble > &primCoeffs,
                              Array<OneD, NekDouble> &coeffs);

            private:
                Expansion2DSharedPtr m_elementLeft;
                Expansion2DSharedPtr m_elementRight;
                int m_elementEdgeLeft;
                int m_elementEdgeRight;

        };
        
        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion1D> Expansion1DSharedPtr;
        typedef std::vector< Expansion1DSharedPtr > Expansion1DVector;
        typedef std::vector< Expansion1DSharedPtr >::iterator Expansion1DVectorIter;
        
        inline Expansion2DSharedPtr Expansion1D::GetLeftAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.get(), "Left adjacent element not set.");
            return m_elementLeft;
        }

        inline Expansion2DSharedPtr Expansion1D::GetRightAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.get(), "Right adjacent element not set.");
            return m_elementRight;
        }

        inline int Expansion1D::GetLeftAdjacentElementEdge() const
        {
            return m_elementEdgeLeft;
        }

        inline int Expansion1D::GetRightAdjacentElementEdge() const
        {
            return m_elementEdgeRight;
        }

        inline void Expansion1D::SetAdjacentElementExp(int edge, Expansion2DSharedPtr &e)
        {
            if (m_elementLeft.get())
            {
                ASSERTL1(!m_elementRight.get(),
                         "Both adjacent elements already set.");
                m_elementRight = e;
                m_elementEdgeRight = edge;
            }
            else
            {
                m_elementLeft = e;
                m_elementEdgeLeft = edge;
            }
        }

    } //end of namespace
} //end of namespace

#define EXPANSION1D_H
#endif

/** 
 *    $Log: Expansion1D.h,v $
 *    Revision 1.2  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
