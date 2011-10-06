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
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions
    {
        class Expansion2D: virtual public Expansion, virtual public StdRegions::StdExpansion2D
        {
        public:
            LOCAL_REGIONS_EXPORT void SetTraceToGeomOrientation(Array<OneD, StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                           Array<OneD, NekDouble> &inout);

            LOCAL_REGIONS_EXPORT void Getnormalindir(const int edge,
                                StdRegions::StdExpansion1DSharedPtr &EdgeExp_e,
                                const Array<OneD, const Array<OneD, NekDouble> > &normals,
                                const Array<OneD, const NekDouble> &directional,
                                Array<OneD, NekDouble> &outarray);

            Expansion1DSharedPtr GetEdgeExp(int edge, bool SetUpNormal=true);

            void SetEdgeExp(const int edge, Expansion1DSharedPtr &e);

            inline void AddNormTraceInt(const int dir,
                                 Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray);

            inline void AddNormTraceInt(const int dir,
                                 Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray,
                                 const Array<OneD, NekDouble> &directional = NullNekDouble1DArray );

            inline void AddEdgeBoundaryInt(const int edge,
                                    const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                    Array <OneD,NekDouble > &outarray);

            inline void AddHDGHelmholtzEdgeTerms(const NekDouble tau,
                                          const int edge,
                                          Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp,
                                          const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                          Array <OneD,NekDouble > &outarray);

            inline void AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                           const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                           Array<OneD,NekDouble> &outarray);

        protected:
            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey);

            // Hybridized DG routines
            virtual void v_AddNormTraceInt(const int dir,
                                 Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray);

            virtual void v_AddNormTraceInt(const int dir,
                                 Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                 Array<OneD,NekDouble> &outarray,
                                 const Array<OneD, NekDouble> &directional = NullNekDouble1DArray );


            virtual void v_AddHDGHelmholtzEdgeTerms(const NekDouble tau,
                                          const int edge,
                                          Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp,
                                          const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                          Array <OneD,NekDouble > &outarray);

            virtual void v_AddEdgeBoundaryInt(const int edge,
                                    const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                    Array <OneD,NekDouble > &outarray);

            virtual void v_DGDeriv(int dir,
                         const Array<OneD, const NekDouble>&incoeffs,
                         Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                         Array<OneD, NekDouble> &out_d);

            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                           const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                           Array<OneD,NekDouble> &outarray);

            virtual void v_AddEdgeNormBoundaryInt(const int edge,
                                        StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                        const Array<OneD, const NekDouble> &Fx,
                                        const Array<OneD, const NekDouble> &Fy,
                                        Array<OneD, NekDouble> &outarray);


            virtual void v_AddEdgeNormBoundaryInt(const int edge,
                                   StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                   const Array<OneD, const NekDouble> &Fn,
                                        Array<OneD, NekDouble> &outarray);

            virtual void v_AddEdgeNormBoundaryBiInt(const int edge,
                                          StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                          const Array<OneD, const NekDouble> &Fwd,
                                          const Array<OneD, const NekDouble> &Bwd,
                                          Array<OneD, NekDouble> &outarray);

            virtual void v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat);

            virtual void v_AddRobinEdgeContribution(const int edgeid, const Array<OneD, const NekDouble> &primCoeffs, Array<OneD, NekDouble> &coeffs);


        private:
            std::vector<Expansion1DSharedPtr> m_edgeExp;

         };

        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        typedef std::vector< Expansion2DSharedPtr > Expansion2DVector;
        typedef std::vector< Expansion2DSharedPtr >::iterator Expansion2DVectorIter;

        inline Expansion1DSharedPtr Expansion2D::GetEdgeExp(int edge, bool SetUpNormal)
        {
            ASSERTL1(edge < GetNedges(), "Edge out of range.");
            return m_edgeExp[edge];
        }

        inline void Expansion2D::SetEdgeExp(const int edge, Expansion1DSharedPtr &e)
        {
            int nEdges = GetNedges();
            ASSERTL1(edge < nEdges, "Edge out of range.");
            if (m_edgeExp.size() < nEdges)
            {
                m_edgeExp.resize(nEdges);
            }
            m_edgeExp[edge] = e;
        }

        inline void Expansion2D::AddNormTraceInt(const int dir,
                             Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                             Array<OneD,NekDouble> &outarray)
        {
            v_AddNormTraceInt(dir, EdgeExp, outarray);
        }

        inline void Expansion2D::AddNormTraceInt(const int dir,
                             Array<OneD, const NekDouble> &inarray,
                             Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                             Array<OneD,NekDouble> &outarray,
                             const Array<OneD, NekDouble> &directional )
        {
            v_AddNormTraceInt(dir, inarray, EdgeExp, outarray, directional);
        }

        inline void Expansion2D::AddEdgeBoundaryInt(const int edge,
                                const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                Array <OneD,NekDouble > &outarray)
        {
            v_AddEdgeBoundaryInt(edge, EdgeExp, outarray);
        }

        inline void Expansion2D::AddHDGHelmholtzEdgeTerms(const NekDouble tau,
                                      const int edge,
                                      Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp,
                                      const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                      Array <OneD,NekDouble > &outarray)
        {
            v_AddHDGHelmholtzEdgeTerms(tau, edge, EdgeExp, dirForcing, outarray);
        }

        inline void Expansion2D::AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                       const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                       Array<OneD,NekDouble> &outarray)
        {
            v_AddHDGHelmholtzTraceTerms(tau, inarray, EdgeExp, dirForcing, outarray);
        }

    } //end of namespace
} //end of namespace

#define EXPANSION2D_H
#endif

/**
 *    $Log: Expansion2D.h,v $
 *    Revision 1.11  2009/11/17 17:43:36  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.10  2009/11/09 15:43:51  sehunchun
 *    HDG2DManifold Solver with Variable coefficients
 *
 *    Revision 1.9  2009/11/06 21:43:56  sherwin
 *    DGDeriv function
 *
 *    Revision 1.8  2009/09/06 22:24:00  sherwin
 *    Updates for Navier-Stokes solver
 *
 *    Revision 1.7  2009/07/07 16:31:47  sehunchun
 *    Adding AddEdgeBoundaryBiInt to line integrate depending on Fwd and Bwd
 *
 *    Revision 1.6  2009/04/02 13:04:36  sherwin
 *    Modified Hybrid Matrix call to use matrix D M^{-1}D' formulation and removed operations based version
 *
 *    Revision 1.5  2008/11/01 22:08:29  bnelson
 *    Fixed compiler warning
 *
 *    Revision 1.4  2008/10/04 19:34:09  sherwin
 *    Added an upwind method which takes the normal flux rather than than individual components
 *
 *    Revision 1.3  2008/08/27 16:35:13  pvos
 *    Small efficiency update
 *
 *    Revision 1.2  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
