///////////////////////////////////////////////////////////////////////////////
//
// File Expansion2D.h
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
// Description: Header file for Expansion2D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION2D_H
#define EXPANSION2D_H

#include <LocalRegions/Expansion.h>
#include <StdRegions/StdExpansion2D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions
    {
      class Expansion3D;
      typedef boost::shared_ptr<Expansion3D> Expansion3DSharedPtr;
      typedef boost::weak_ptr<Expansion3D> Expansion3DWeakPtr;

        class Expansion2D: virtual public Expansion, virtual public StdRegions::StdExpansion2D
        {
        public:
            LOCAL_REGIONS_EXPORT Expansion2D();
            LOCAL_REGIONS_EXPORT virtual ~Expansion2D() {}
            
            LOCAL_REGIONS_EXPORT void SetTraceToGeomOrientation(Array<OneD, StdRegions::StdExpansionSharedPtr> &EdgeExp,
                                           Array<OneD, NekDouble> &inout);

            ExpansionSharedPtr GetEdgeExp(int edge, bool SetUpNormal=true);
            
            void SetEdgeExp(const int edge, ExpansionSharedPtr &e);

            inline void AddNormTraceInt(const int dir,
                                 Array<OneD,StdRegions::StdExpansionSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray);

            inline void AddNormTraceInt(const int dir,
                                 Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,StdRegions::StdExpansionSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray,
                                 const StdRegions::VarCoeffMap &varcoeffs);

            inline void AddEdgeBoundaryInt(const int edge,
                                    StdRegions::StdExpansionSharedPtr &EdgeExp,
                                    Array <OneD,NekDouble > &outarray,
                                    const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

            inline void AddHDGHelmholtzEdgeTerms(const NekDouble tau,
                                          const int edge,
                                          Array <OneD, StdRegions::StdExpansionSharedPtr > &EdgeExp,
                                          const StdRegions::VarCoeffMap &dirForcing,
                                          Array <OneD,NekDouble > &outarray);

            inline void AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,StdRegions::StdExpansionSharedPtr> &EdgeExp,
                                           const StdRegions::VarCoeffMap &dirForcing,
                                           Array<OneD,NekDouble> &outarray);

            inline Expansion3DSharedPtr GetLeftAdjacentElementExp() const;

            inline Expansion3DSharedPtr GetRightAdjacentElementExp() const;

            inline int GetLeftAdjacentElementFace() const;

            inline int GetRightAdjacentElementFace() const;

            inline void SetAdjacentElementExp(
                        int face,
                        Expansion3DSharedPtr &f);

        protected:
            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey);

            // Hybridized DG routines
            virtual void v_DGDeriv(const int dir,
                         const Array<OneD, const NekDouble>&incoeffs,
                         Array<OneD,StdRegions::StdExpansionSharedPtr> &EdgeExp,
                         Array<OneD, NekDouble> &out_d);

            virtual void v_AddEdgeNormBoundaryInt(
                const int                            edge,
                StdRegions::StdExpansionSharedPtr   &EdgeExp,
                const Array<OneD, const NekDouble>  &Fx,
                const Array<OneD, const NekDouble>  &Fy,
                      Array<OneD,       NekDouble>  &outarray);

            virtual void v_AddEdgeNormBoundaryInt(
                const int                            edge,
                StdRegions::StdExpansionSharedPtr   &EdgeExp,
                const Array<OneD, const NekDouble>  &Fn,
                      Array<OneD,       NekDouble>  &outarray);

            virtual void v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat);

            virtual void v_AddRobinEdgeContribution(const int edgeid, const Array<OneD, const NekDouble> &primCoeffs, Array<OneD, NekDouble> &coeffs);

            void GetPhysEdgeVarCoeffsFromElement(
                    const int edge,
                    StdRegions::StdExpansionSharedPtr &EdgeExp,
                    const Array<OneD, const NekDouble>  &varcoeff,
                    Array<OneD,NekDouble> &outarray);

        private:
            std::vector<ExpansionWeakPtr> m_edgeExp;

            Expansion3DWeakPtr m_elementLeft;
            Expansion3DWeakPtr m_elementRight;
            int m_elementFaceLeft;
            int m_elementFaceRight;

         };

        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        typedef boost::weak_ptr<Expansion2D> Expansion2DWeakPtr;
        typedef std::vector< Expansion2DSharedPtr > Expansion2DVector;
        typedef std::vector< Expansion2DSharedPtr >::iterator Expansion2DVectorIter;

        inline ExpansionSharedPtr Expansion2D::GetEdgeExp(int edge, bool SetUpNormal)
        {
            ASSERTL1(edge < GetNedges(), "Edge out of range.");
            return m_edgeExp[edge].lock();
        }

        inline void Expansion2D::SetEdgeExp(const int edge, ExpansionSharedPtr &e)
        {
            unsigned int nEdges = GetNedges();
            ASSERTL1(edge < nEdges, "Edge out of range.");
            if (m_edgeExp.size() < nEdges)
            {
                m_edgeExp.resize(nEdges);
            }
            m_edgeExp[edge] = e;
        }

        inline Expansion3DSharedPtr Expansion2D::GetLeftAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(), "Left adjacent element not set.");
            return m_elementLeft.lock();
        }
        
        inline Expansion3DSharedPtr Expansion2D::GetRightAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(), "Right adjacent element not set.");
            return m_elementRight.lock();
        }

        inline int Expansion2D::GetLeftAdjacentElementFace() const
        {
            return m_elementFaceLeft;
        }

        inline int Expansion2D::GetRightAdjacentElementFace() const
        {
            return m_elementFaceRight;
        }
        
        inline void Expansion2D::SetAdjacentElementExp(int face, Expansion3DSharedPtr &f)
        {
            if (m_elementLeft.lock().get())
            {
                ASSERTL1(!m_elementRight.lock().get(),
                         "Both adjacent elements already set.");
                m_elementRight = f;
                m_elementFaceRight = face;
            }
            else
            {
                m_elementLeft = f;
                m_elementFaceLeft = face;
            }
	}
    } //end of namespace
} //end of namespace

#endif
