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
	typedef boost::shared_ptr<StdExpansion3D> StdExpansion3DSharedPtr;

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


        protected:

            /** \brief This function evaluates the expansion at a single
             *  (arbitrary) point of the domain
             *
             *  This function is a wrapper around the virtual function
             *  \a v_PhysEvaluate()
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
             *  \f$\mathbf{u}\f$ (implemented as the attribute #m_phys)
             *  is set.
             *
             *  \param coords the coordinates of the single point
             *  \return returns the value of the expansion at the single point
             */
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coords,
                    const Array<OneD, const NekDouble>& physvals);

            STD_REGIONS_EXPORT virtual void v_NegateFaceNormal(
                const int face);
            NormalVector m_surfaceNormal;

            std::map<int, NormalVector> m_faceNormals;
            std::map<int, bool> m_negatedNormals;

        private:

            virtual int v_GetShapeDimension() const
            {
                return 3;
            }

            virtual int v_GetCoordim(void)
            {
                return 3;
            }
            STD_REGIONS_EXPORT const NormalVector & v_GetSurfaceNormal() const;
            STD_REGIONS_EXPORT const NormalVector & v_GetFaceNormal(const int face) const;
        };
    } //end of namespace
} //end of namespace

#endif //STDEXP3D_H

/**
 * $Log: StdExpansion3D.h,v $
 * Revision 1.19  2008/09/16 13:37:03  pvos
 * Restructured the LocalToGlobalMap classes
 *
 * Revision 1.18  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.17  2008/06/06 23:21:13  ehan
 * Added doxygen documentation
 *
 * Revision 1.16  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.15  2008/05/15 22:39:34  ehan
 * Clean up the codes
 *
 * Revision 1.14  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.13  2008/02/12 02:46:33  jfrazier
 * Moved typedef to the top of the file.
 *
 * Revision 1.12  2008/02/12 01:30:07  ehan
 * Added typedef StdExpansion3DSharedPtr.
 *
 * Revision 1.11  2007/11/08 14:27:53  ehan
 * Fixed PhysTensorDerivative3D matrix and improved L1 error up to 1e-15.
 *
 * Revision 1.10  2007/10/29 20:31:04  ehan
 * Fixed floating point approximation up to 1-e15 for PhysEvaluate.
 *
 * Revision 1.9  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.8  2007/05/15 05:18:23  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.7  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.6  2007/03/29 19:35:09  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.5  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.4  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.3  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.2  2006/06/13 18:05:02  sherwin
 * Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.12  2006/04/25 20:23:33  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.11  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.10  2006/03/06 17:12:45  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.9  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.8  2006/03/01 08:25:03  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.7  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 **/

