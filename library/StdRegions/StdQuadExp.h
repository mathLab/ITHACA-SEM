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
            ShapeType DetShapeType()
            {
                return eQuadrilateral;
            };


            /** \brief Fill outarray with mode \a mode of expansion
            *
            *	Note for quadrilateral expansions _base[0] (i.e. p)  modes run 
            *  fastest
            */
            void FillMode(const int mode, NekDoubleSharedArray &outarray);

            //////////////////////////////
            // Integration Methods
            //////////////////////////////

            NekDouble Integral(const NekDoubleSharedArray &inarray);

            void IProductWRTBase(const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &outarray);

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
            void IProductWRTBase(const NekDouble *base0, 
                const NekDouble *base1,
                const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &outarray,
                int coll_check);

            //----------------------------------
            // Local Matrix Routines
            //----------------------------------

            DNekMatSharedPtr GenMassMatrix();
            DNekMatSharedPtr GenLapMatrix();

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the derivative of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
            void PhysDeriv(const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &out_d0, 
                NekDoubleSharedArray &out_d1,
                NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray);


            //----------------------------
            // Evaluations Methods
            //-----------------------------

            void BwdTrans(const NekDoubleSharedArray &inarray,
                NekDoubleSharedArray &outarray);
            void FwdTrans(const NekDoubleSharedArray &inarray,
                NekDoubleSharedArray &outarray);

            double PhysEvaluate(const NekDoubleSharedArray &coords);

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

            const int GetEdgeNcoeffs(const int i)
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

            const LibUtilities::BasisType GetEdgeBasisType(const int i)
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

            void GetCoords(NekDoubleSharedArray &coords_0, 
                NekDoubleSharedArray &coords_1);

        protected:


        private:

            virtual int v_GetNverts()
            {
                return 4;
            }

            virtual int v_GetNedges()
            {
                return 4;
            }

            virtual int v_GetEdgeNcoeffs(const int i)
            {
                return GetEdgeNcoeffs(i);
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i)
            {
                return GetEdgeBasisType(i);
            }

            virtual ShapeType v_DetShapeType()
            {
                return DetShapeType();
            }

            virtual void v_GetCoords(NekDoubleSharedArray &coords_0,
                NekDoubleSharedArray &coords_1,
                NekDoubleSharedArray &coords_2)
            {
                GetCoords(coords_0,coords_1);
            }


            virtual void v_FillMode(const int mode, NekDoubleSharedArray &array)
            {
                FillMode(mode,array);
            }

            virtual double v_Integral(const NekDoubleSharedArray &inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const NekDoubleSharedArray &inarray,
                NekDoubleSharedArray &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual DNekMatSharedPtr v_GenMassMatrix()
            {
                return GenMassMatrix();
            }

            virtual DNekMatSharedPtr v_GenLapMatrix()
            {
                GenLapMatrix();
            }

            virtual void v_PhysDeriv(const NekDoubleSharedArray &inarray,
                NekDoubleSharedArray &out_d0,
                NekDoubleSharedArray &out_d1 = NullNekDoubleSharedArray,
                NekDoubleSharedArray &out_d2 = NullNekDoubleSharedArray)
            {
                PhysDeriv(inarray,out_d0, out_d1);
            }

            virtual void v_StdPhysDeriv(const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &out_d0,
                NekDoubleSharedArray &out_d1)
            {
                PhysDeriv(inarray, out_d0,  out_d1);
            }

            virtual void v_BwdTrans(const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &outarray)
            {
                BwdTrans(inarray, outarray);
            }

            virtual void v_FwdTrans(const NekDoubleSharedArray &inarray, 
                NekDoubleSharedArray &outarray)
            {
                FwdTrans(inarray, outarray);
            }

            virtual NekDouble v_PhysEvaluate(const NekDoubleSharedArray &Lcoords)
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




