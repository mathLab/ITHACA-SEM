///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion1D.h
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
// which are common to 1d expansion. Typically this inolves physical
// space operations.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef STDEXP1D_H
#define STDEXP1D_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion1D: virtual public StdExpansion
        {

        public:

            STD_REGIONS_EXPORT StdExpansion1D();
            STD_REGIONS_EXPORT StdExpansion1D(int numcoeffs, const LibUtilities::BasisKey &Ba);
            STD_REGIONS_EXPORT StdExpansion1D(const StdExpansion1D &T);
            STD_REGIONS_EXPORT virtual ~StdExpansion1D();

            /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
            *  physical quadrature points given by \a inarray and return in
            *  \a outarray.
            *
            *  \param inarray array of a function evaluated at the quadrature
            *  points
            *  \param outarray the resulting array of the derivative \f$
            *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array
            *  \a outarray as output of the function
            */
            STD_REGIONS_EXPORT void PhysTensorDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,       NekDouble>& outarray);

        protected:
            std::map<int, NormalVector> m_vertexNormals;
            
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coords,
                    const Array<OneD, const NekDouble>& physvals);


        private:

            // Virtual Functions ----------------------------------------

            virtual int v_GetCoordim(void)
            {
                return 1;
            }

            virtual int v_GetShapeDimension() const
            {
                return 1;
            }

            virtual int v_GetNedges() const
            {
                return 0;
            }

            virtual int v_GetNfaces() const
            {
                return 0;
            }

            virtual void v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                                  Array<OneD, const NekDouble> &inarray,
                                                  Array<OneD, NekDouble> &outarray)
            {
                ASSERTL0(false,"method only valid in local regions 1D classes");
            }

            virtual void v_SetCoeffsToOrientation(StdRegions::Orientation dir)
            {
                ASSERTL0(false,"method only valid in local regions 1D classes");
            }
			
			STD_REGIONS_EXPORT virtual void v_SetUpPhysNormals(const int vertex);
						
			STD_REGIONS_EXPORT const NormalVector & v_GetVertexNormal(const int vertex) const;

        };

        typedef boost::shared_ptr<StdExpansion1D> StdExpansion1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDEXP1D_H

/**
* $Log: StdExpansion1D.h,v $
* Revision 1.32  2008/09/16 13:37:03  pvos
* Restructured the LocalToGlobalMap classes
*
* Revision 1.31  2008/08/14 22:09:50  sherwin
* Modifications to remove HDG routines from StdRegions and removed StdExpMap
*
* Revision 1.30  2008/07/29 22:21:15  sherwin
* A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
*
* Revision 1.29  2008/07/04 10:18:40  pvos
* Some updates
*
* Revision 1.28  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.27  2008/04/06 06:04:15  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.26  2008/04/02 22:18:10  pvos
* Update for 2D local to global mapping
*
* Revision 1.25  2008/02/29 19:15:19  sherwin
* Update for UDG stuff
*
* Revision 1.24  2007/12/06 22:44:47  pvos
* 2D Helmholtz solver updates
*
* Revision 1.23  2007/11/08 16:55:14  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.22  2007/07/22 23:04:26  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.21  2007/07/20 02:16:53  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.20  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.19  2007/05/17 17:59:28  sherwin
* Modification to make Demos work after introducion of Array<>
*
* Revision 1.18  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.17  2007/04/08 03:36:58  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.16  2007/04/04 20:48:17  sherwin
* Update to handle SharedArrays
*
* Revision 1.15  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.14  2007/03/25 15:48:22  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.13  2007/03/21 20:56:43  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.12  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.11  2007/03/20 09:12:47  kirby
* update of geofac and metric info; fix style issues
*
* Revision 1.10  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.9  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.8  2007/02/07 12:51:53  sherwin
* Compiling version of Project1D
*
* Revision 1.7  2007/01/28 18:34:23  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.6  2007/01/21 02:28:07  sherwin
* Compiling under new revision
*
* Revision 1.5  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.4  2007/01/15 11:30:22  pvos
* Updating doxygen documentation
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
* Revision 1.17  2006/05/02 21:21:12  sherwin
* Corrected libraries to compile new version of spatialdomains and demo Graph1D
*
* Revision 1.16  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.15  2006/03/13 18:29:35  sherwin
*
* Corrected error with definition of GetCoords
*
* Revision 1.14  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.13  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.12  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.11  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.10  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/



