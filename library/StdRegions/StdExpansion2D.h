///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.h
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
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP2D_H
#define STDEXP2D_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

    class StdExpansion2D: virtual public StdExpansion
    {
        public:
            STD_REGIONS_EXPORT StdExpansion2D();
            STD_REGIONS_EXPORT StdExpansion2D(int numcoeffs, const LibUtilities::BasisKey &Ba,
                           const LibUtilities::BasisKey &Bb);
            STD_REGIONS_EXPORT StdExpansion2D(const StdExpansion2D &T);
            STD_REGIONS_EXPORT ~StdExpansion2D();

            // Generic operations in different element

            /** \brief Calculate the 2D derivative in the local
             *  tensor/collapsed coordinate at the physical points
             *
             *  This function is independent of the expansion basis and can
             *  therefore be defined for all tensor product distribution of
             *  quadrature points in a generic manner.  The key operations are:
             *
             *  - \f$ \frac{d}{d\eta_1} \rightarrow {\bf D^T_0 u } \f$ \n
             *  - \f$ \frac{d}{d\eta_2} \rightarrow {\bf D_1 u } \f$
             *
             *  \param inarray array of physical points to be differentiated
             *  \param  outarray_d0 the resulting array of derivative in the
             *  \f$\eta_1\f$ direction will be stored in outarray_d0 as output
             *  of the function
             *  \param outarray_d1 the resulting array of derivative in the
             *  \f$\eta_2\f$ direction will be stored in outarray_d1 as output
             *  of the function
             *
             *  Recall that:
             *  \f$
             *  \hspace{1cm} \begin{array}{llll}
             *  \mbox{Shape}    & \mbox{Cartesian coordinate range} &
             *  \mbox{Collapsed coord.}      &
             *  \mbox{Collapsed coordinate definition}\\
             *  \mbox{Quadrilateral}  & -1 \leq \xi_1,\xi_2 \leq  1
             *  & -1 \leq \eta_1,\eta_2 \leq 1
             *  & \eta_1 = \xi_1, \eta_2 = \xi_2\\
             *  \mbox{Triangle}  & -1 \leq \xi_1,\xi_2; \xi_1+\xi_2 \leq  0
             *  & -1 \leq \eta_1,\eta_2 \leq 1
             *  & \eta_1 = \frac{2(1+\xi_1)}{(1-\xi_2)}-1, \eta_2 = \xi_2 \\
             *  \end{array} \f$
             */
            STD_REGIONS_EXPORT void PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray_d0,
                                 Array<OneD, NekDouble> &outarray_d1);

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
            STD_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble> & physvals);

            STD_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble>& inarray,
                               const Array<OneD, const NekDouble>& w0,
                               const Array<OneD, const NekDouble>& w1);


        protected:
            std::map<int, NormalVector> m_edgeNormals;
            NormalVector m_surfaceNormal;

        private:

            // Virtual Functions ----------------------------------------
            virtual int v_GetShapeDimension() const
            {
                return 2;
            }

            virtual int v_GetNfaces() const
            {
                return 0;
            }

            virtual int v_GetCoordim(void)
            {
                return 2;
            }

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
            {
                return PhysEvaluate(coords);
            }

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble> & physvals)
            {
                return PhysEvaluate(coords,physvals);
            }

            STD_REGIONS_EXPORT virtual void v_SetUpPhysNormals(const int edge);

            STD_REGIONS_EXPORT const NormalVector & v_GetEdgeNormal(const int edge) const;

    };

        typedef boost::shared_ptr<StdExpansion2D> StdExpansion2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDEXP2D_H

/**
 * $Log: StdExpansion2D.h,v $
 * Revision 1.25  2008/09/16 13:37:03  pvos
 * Restructured the LocalToGlobalMap classes
 *
 * Revision 1.24  2008/08/14 22:09:50  sherwin
 * Modifications to remove HDG routines from StdRegions and removed StdExpMap
 *
 * Revision 1.23  2008/07/29 22:21:15  sherwin
 * A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
 *
 * Revision 1.22  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.21  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.20  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.19  2008/04/02 22:18:10  pvos
 * Update for 2D local to global mapping
 *
 * Revision 1.18  2007/11/08 16:55:14  pvos
 * Updates towards 2D helmholtz solver
 *
 * Revision 1.17  2007/10/15 20:37:40  ehan
 * Make changes of column major matrix
 *
 * Revision 1.16  2007/07/22 23:04:26  bnelson
 * Backed out Nektar::ptr.
 *
 * Revision 1.15  2007/07/20 02:16:53  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.14  2007/05/30 20:49:13  sherwin
 * Updates to do with LocalRegions and SpatialDomains
 *
 * Revision 1.13  2007/05/15 05:18:23  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.12  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.11  2007/04/06 08:44:43  sherwin
 * Update to make 2D regions work at StdRegions level
 *
 * Revision 1.10  2007/04/05 15:20:11  sherwin
 * Updated 2D stuff to comply with SharedArray philosophy
 *
 * Revision 1.9  2007/03/29 19:35:09  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.8  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.7  2007/03/14 21:24:09  sherwin
 * Update for working version of MultiRegions up to ExpList1D
 *
 * Revision 1.6  2007/03/05 19:26:56  bcarmo
 * StdExpansion2D.h modified according to StdExpansion1D. Compiles.
 *
 * Revision 1.6  2007/01/17 16:05:39  bcarmo
 * Version with StdExpansion2D compiling
 *
 * Revision 1.5  2007/01/17 16:05:39  pvos
 * updated doxygen documentation
 *
 * Revision 1.4  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
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
 * Revision 1.19  2006/05/02 21:21:12  sherwin
 * Corrected libraries to compile new version of spatialdomains and demo Graph1D
 *
 * Revision 1.18  2006/04/25 20:23:33  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.17  2006/03/06 17:12:45  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.16  2006/03/05 22:11:02  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.15  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.14  2006/03/01 22:59:12  sherwin
 *
 * First working version of Project1D
 *
 * Revision 1.13  2006/03/01 08:25:03  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.12  2006/02/26 23:37:29  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/



