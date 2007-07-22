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

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion1D: public StdExpansion
        {

        public:

            StdExpansion1D();
            StdExpansion1D(int numcoeffs, const LibUtilities::BasisKey &Ba);
            StdExpansion1D(const StdExpansion1D &T);
            ~StdExpansion1D();


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
            void PhysTensorDeriv(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble>& outarray);

            void PhysDeriv (const ConstArray<OneD, NekDouble>& inarray,
                            Array<OneD, NekDouble>& out_d1 = NullNekDouble1DArray,
                            Array<OneD, NekDouble>& out_d2 = NullNekDouble1DArray,
                            Array<OneD, NekDouble>& out_d3 = NullNekDouble1DArray)
            {
                v_PhysDeriv (inarray, out_d1, out_d2, out_d3);
            }

            void StdPhysDeriv (const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble>& outarray)
            {
                v_StdPhysDeriv (inarray,outarray);
            }

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
            NekDouble PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
            {
                return v_PhysEvaluate(coords);
            }



            /** \brief Evaluate a function at points coords which is assumed
            *  to be in local collapsed coordinate format. The function is
            *  assumed to be in physical space
            */
            NekDouble PhysEvaluate1D(const ConstArray<OneD, NekDouble>& coords);

            /** \brief wrapper around virtual call */
            void FwdTrans(const StdExpansion1D &in)
            {
                v_FwdTrans(in);
            }

            void BwdTrans(const StdExpansion1D &in)
            {
                v_BwdTrans(in);
            }

            void  BwdTrans (const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray)
            {
                v_BwdTrans (inarray, outarray);
            }




        protected:

        private:

            // Virtual Functions ----------------------------------------

            virtual void v_FwdTrans(const StdExpansion1D &in) = 0;
            virtual void v_BwdTrans(const StdExpansion1D &in) = 0;

            virtual void v_BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray) = 0;


            virtual int v_GetNverts() const = 0;
            virtual int v_GetNedges() const
            {
                ASSERTL0(false,"This function is only valid for 2 and 3D expansions");
                return 0;
            }
            virtual int v_GetNfaces() const
            {
                ASSERTL0(false,"This function is only valid for 2 and 3D expansions");
                return 0;
            }

            virtual ShapeType v_DetShapeType() const = 0;

            virtual int v_get_nodalpoints(const ConstArray<OneD, NekDouble>& x, 
                Array<OneD, ConstArray<OneD, NekDouble> >& y)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return 0;
            }

            virtual void v_GenNBasisTransMatrix(Array<OneD, NekDouble> & outarray)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
            }


            virtual int v_GetCoordim(void)
            {
                return 1; 
            }

            virtual void   v_PhysDeriv (const ConstArray<OneD, NekDouble>& inarray,
                                        Array<OneD, NekDouble> &out_d0,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2) = 0;
            
            virtual void   v_StdPhysDeriv (const ConstArray<OneD, NekDouble>& inarray, 
                                           Array<OneD, NekDouble> &outarray) = 0;

            virtual void   v_StdPhysDeriv (const int dir, 
                                           const ConstArray<OneD, NekDouble>& inarray, 
                                           Array<OneD, NekDouble> &outarray) = 0;

            virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "local expansions");
                return 0;
            }
        };

        typedef boost::shared_ptr<StdExpansion1D> StdExpansion1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDEXP1D_H

/**
* $Log: StdExpansion1D.h,v $
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



