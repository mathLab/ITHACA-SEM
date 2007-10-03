///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.h
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
// Description: Header field for Quadrilateral routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDQUADEXP_H
#define NEKTAR_LIB_STDREGIONS_STDQUADEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdExpMap.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdQuadExp: public StdExpansion2D
        {

        public:

            StdQuadExp();

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition 
            */
            StdQuadExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb);

            /** \brief Copy Constructor */
            StdQuadExp(const StdQuadExp &T);

            /** \brief Destructor */
            ~StdQuadExp();

            /** \brief Return Shape of region, using ShapeType enum list. 
            *  i.e. Quadrilateral
            */
            ShapeType DetShapeType() const
            {
                return eQuadrilateral;
            };


            /** \brief Fill outarray with mode \a mode of expansion
            *
            *    Note for quadrilateral expansions _base[0] (i.e. p)  modes run 
            *  fastest
            */
            void FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            //////////////////////////////
            // Integration Methods
            //////////////////////////////

            NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray);

            void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray);

            /** \brief Calculate the inner product of inarray with respect to
            *  the basis B=base0*base1 and put into outarray
            *
            *  \f$ 
            *  \begin{array}{rcl}
            *  I_{pq} = (\phi_q \phi_q, u) & = & \sum_{i=0}^{nq_0}
            *  \sum_{j=0}^{nq_1}
            *  \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} 
            *  \xi_{1,j}) \\
            *  & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
            *  \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} 
            *  \end{array}
            *  \f$ 
            *
            *  where
            *
            *  \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
            *
            *  which can be implemented as
            *
            *  \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) 
            *  \tilde{u}_{i,j} = {\bf B_1 U}  \f$
            *  \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
            *  {\bf B_0 F}  \f$
            */
            void IProductWRTBase(const ConstArray<OneD, NekDouble>& base0, 
                const ConstArray<OneD, NekDouble>& base1,
                const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray,
                int coll_check);

            //----------------------------------
            // Generate Matrix Routine
            //----------------------------------

            DNekMatSharedPtr GenMatrix(MatrixType mtype);

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the derivative of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
            void PhysDeriv(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &out_d0, 
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);


            //----------------------------
            // Evaluations Methods
            //-----------------------------

            void BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);
            void FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);

            NekDouble PhysEvaluate(ConstArray<OneD, NekDouble>& coords);

            void MapTo(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient, 
                StdExpMap &Map);

            void MapTo_ModalFormat(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient, 
                StdExpMap &Map);

            const int GetEdgeNcoeffs(const int i) const
            {
                ASSERTL2((i > 0)&&(i < 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisNumModes(0);
                }
                else
                {
                    return  GetBasisNumModes(1); 
                }

            }

            const LibUtilities::BasisType GetEdgeBasisType(const int i) const
            {
                ASSERTL2((i > 0)&&(i < 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisType(0);
                }
                else
                {
                    return  GetBasisType(1);
                }

            }

            void GetCoords(Array<OneD, NekDouble> &coords_0, 
                Array<OneD, NekDouble> &coords_1);

        protected:


        private:

            virtual int v_GetNverts() const
            {
                return 4;
            }

            virtual int v_GetNedges() const
            {
                return 4;
            }

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual int v_NumBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                return 4 + 2*(GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
            } 

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            virtual ShapeType v_DetShapeType() const
            {
                return DetShapeType();
            }

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
            {
                GetCoords(coords_0,coords_1);
            }


            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &array)
            {
                FillMode(mode,array);
            }

            virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual DNekMatSharedPtr v_GenMatrix(MatrixType mtype)
            {
                return GenMatrix(mtype);
            }

            virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray,out_d0, out_d1);
            }

            virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1)
            {
                PhysDeriv(inarray, out_d0,  out_d1);
            }

            virtual void v_BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray, outarray);
            }

            virtual void v_FwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray, outarray);
            }

            virtual NekDouble v_PhysEvaluate(ConstArray<OneD, NekDouble>& Lcoords)
            {
                return PhysEvaluate(Lcoords);
            }

            virtual void v_MapTo(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient,
                StdExpMap &Map)
            {
                MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
            }

            virtual void v_MapTo_ModalFormat(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient,
                StdExpMap &Map)
            {
                MapTo_ModalFormat(edge_ncoeffs,Btype,eid,eorient,Map);
            }

        };

    } //end of namespace
} //end of namespace

#endif //STDQUADEXP_H

/**
* $Log: StdQuadExp.h,v $
* Revision 1.19  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.18  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.17  2007/07/10 20:41:51  kirby
* more fixes
*
* Revision 1.16  2007/07/10 19:27:57  kirby
* Update for new matrix structures
*
* Revision 1.15  2007/05/31 19:13:12  pvos
* Updated NodalTriExp + LocalRegions/Project2D + some other modifications
*
* Revision 1.14  2007/05/15 05:18:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.13  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.12  2007/04/08 03:36:58  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.11  2007/04/06 08:44:43  sherwin
* Update to make 2D regions work at StdRegions level
*
* Revision 1.10  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.9  2007/04/05 11:40:21  pvincent
* *** empty log message ***
*
* Revision 1.8  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.7  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.6  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.5  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.4  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.3  2006/06/01 14:13:36  kirby
* *** empty log message ***
*
* Revision 1.2  2006/05/23 15:08:19  jfrazier
* Minor cleanup to correct compile warnings.
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.34  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.33  2006/03/05 23:17:53  sherwin
*
* Corrected to allow MMatrix1D and MMatrix2D to execute properly
*
* Revision 1.32  2006/03/04 20:26:55  bnelson
* Added comments after #endif.
*
* Revision 1.31  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.30  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/




