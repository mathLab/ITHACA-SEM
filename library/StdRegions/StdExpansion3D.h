///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion3D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP3D_H
#define STDEXP3D_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

	class StdExpansion3D;
	typedef std::shared_ptr<StdExpansion3D> StdExpansion3DSharedPtr;

    class StdExpansion3D: virtual public StdExpansion
        {

        public:
            STD_REGIONS_EXPORT StdExpansion3D();
            STD_REGIONS_EXPORT StdExpansion3D(int numcoeffs, const LibUtilities::BasisKey &Ba,
                           const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc);
            STD_REGIONS_EXPORT StdExpansion3D(const StdExpansion3D &T);
            STD_REGIONS_EXPORT virtual ~StdExpansion3D();

            // Differentiation

            /** \brief Calculate the 3D derivative in the local
             *  tensor/collapsed coordinate at the physical points
             *
             *    This function is independent of the expansion basis and can
             *    therefore be defined for all tensor product distribution of
             *    quadrature points in a generic manner.  The key operations are:
             *
             *    - \f$ \frac{d}{d\eta_1} \rightarrow {\bf D^T_0 u } \f$ \n
             *    - \f$ \frac{d}{d\eta_2} \rightarrow {\bf D_1 u } \f$
             *    - \f$ \frac{d}{d\eta_3} \rightarrow {\bf D_2 u } \f$
             *
             *  \param inarray array of physical points to be differentiated
             *  \param  outarray_d1 the resulting array of derivative in the
             *  \f$\eta_1\f$ direction will be stored in outarray_d1 as output
             *  of the function
             *  \param outarray_d2 the resulting array of derivative in the
             *  \f$\eta_2\f$ direction will be stored in outarray_d2 as output
             *  of the function
             *  \param outarray_d3 the resulting array of derivative in the
             *  \f$\eta_3\f$ direction will be stored in outarray_d3 as output
             *  of the function
             *
             *  Recall that:
             *  \f$
             *  \hspace{1cm} \begin{array}{llll}
             *  \mbox{Shape}    & \mbox{Cartesian coordinate range} &
             *  \mbox{Collapsed coord.}      &
             *  \mbox{Collapsed coordinate definition}\\
             *  \mbox{Hexahedral}  & -1 \leq \xi_1,\xi_2, \xi_3 \leq  1
             *  & -1 \leq \eta_1,\eta_2, \eta_3 \leq 1
             *  & \eta_1 = \xi_1, \eta_2 = \xi_2, \eta_3 = \xi_3 \\
             *  \mbox{Tetrahedral}  & -1 \leq \xi_1,\xi_2,\xi_3; \xi_1+\xi_2 +\xi_3 \leq  -1
             *  & -1 \leq \eta_1,\eta_2, \eta_3 \leq 1
             *  & \eta_1 = \frac{2(1+\xi_1)}{-\xi_2 -\xi_3}-1, \eta_2 = \frac{2(1+\xi_2)}{1 - \xi_3}-1, \eta_3 = \xi_3 \\
             *  \end{array} \f$
             */
            STD_REGIONS_EXPORT void PhysTensorDeriv(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray_d1,
                                 Array<OneD, NekDouble> &outarray_d2,
                                 Array<OneD, NekDouble> &outarray_d3);

            STD_REGIONS_EXPORT void BwdTrans_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray,
                      Array<OneD,       NekDouble>& wsp,
                bool                                doCheckCollDir0,
                bool                                doCheckCollDir1,
                bool                                doCheckCollDir2);

            STD_REGIONS_EXPORT void IProductWRTBase_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2);

        protected:

            /** \brief This function evaluates the expansion at a single
             *  (arbitrary) point of the domain
             *
             *
             *  Based on the value of the expansion at the quadrature points,
             *  this function calculates the value of the expansion at an
             *  arbitrary single points (with coordinates \f$ \mathbf{x_c}\f$
             *  given by the pointer \a coords). This operation, equivalent to
             *  \f[ u(\mathbf{x_c})  = \sum_p \phi_p(\mathbf{x_c}) \hat{u}_p \f]
             *  is evaluated using Lagrangian interpolants through the quadrature
             *  points:
             *  \f[ u(\mathbf{x_c}) = \sum_p h_p(\mathbf{x_c}) u_p\f]
             *
             *  This function requires that the physical value array
             *  \f$\mathbf{u}\f$ (implemented as the attribute #phys)
             *  is set.
             *
             *  \param coords the coordinates of the single point
             *  \return returns the value of the expansion at the single point
             */
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coords,
                    const Array<OneD, const NekDouble>& physvals);


            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, DNekMatSharedPtr >& I,
                    const Array<OneD, const NekDouble >& physvals);

            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray,
                      Array<OneD,       NekDouble>& wsp,
                bool                                doCheckCollDir0,
                bool                                doCheckCollDir1,
                bool                                doCheckCollDir2) = 0;

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble>&       outarray,
                      Array<OneD, NekDouble>&       wsp,
                bool                                doCheckCollDir0,
                bool                                doCheckCollDir1,
                bool                                doCheckCollDir2) = 0;

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp_MatFree(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp_MatFree(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
                const Array<OneD, const NekDouble>& inarray);

            STD_REGIONS_EXPORT virtual int v_GetTraceNcoeffs(const int i) const
            {
                return GetFaceNcoeffs(i);
            }

            std::map<int, NormalVector> m_faceNormals;

        private:

            virtual int v_GetShapeDimension() const
            {
                return 3;
            }

            virtual int v_GetCoordim(void)
            {
                return 3;
            }
            STD_REGIONS_EXPORT const NormalVector & v_GetSurfaceNormal(const int id) const;
            STD_REGIONS_EXPORT const NormalVector & v_GetFaceNormal(const int face) const;
            
        };

        STD_REGIONS_EXPORT LibUtilities::BasisKey EvaluateTriFaceBasisKey(
            const int                     facedir,
            const LibUtilities::BasisType faceDirBasisType,
            const int                     numpoints,
            const int                     nummodes);

        STD_REGIONS_EXPORT LibUtilities::BasisKey EvaluateQuadFaceBasisKey(
            const int                     facedir,
            const LibUtilities::BasisType faceDirBasisType,
            const int                     numpoints,
            const int                     nummodes);
    } //end of namespace
} //end of namespace

#endif //STDEXP3D_H
