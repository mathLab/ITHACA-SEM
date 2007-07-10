///////////////////////////////////////////////////////////////////////////////
//
// File StdTriExp.h
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
// Description: Header field for triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDTRIEXP_H
#define NEKTAR_LIB_STDREGIONS_STDTRIEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdMatrixKey.h>
#include <StdRegions/StdSegExp.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdTriExp: public StdExpansion2D
        {

        public:

            StdTriExp();

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition 
            */
            StdTriExp(const LibUtilities::BasisKey &Ba,
                const LibUtilities::BasisKey &Bb);

            /** \brief Constructor using BasisKey class for quadrature points
            *  and order definition where m_coeffs and m_phys are all set.
            */

            /** \brief Copy Constructor */
            StdTriExp(const StdTriExp &T);

            /** \brief Destructor */
            ~StdTriExp();

            /** \brief Return Shape of region, using  ShapeType enum list.
            *  i.e. Triangle
            */
            ShapeType DetShapeType()
            {
                return eTriangle;
            }

            //////////////////////////////
            // Integration Methods
            //////////////////////////////
            /**
            *  This is just a wrapper around StdExpansoin2D which multiplies by
            *  0.5 which is due to the factor \f$ (1-\xi_2)/2 \f$ in the 
            *  integral weight
            */
            NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray);

            void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray);

            /** \brief Fill outarray with mode \a mode of expansion
            *
            * Note for quadrilateral expansions _base[0] (i.e. p)  modes run 
            *  fastest
            */
            void FillMode(const int mode, 
                Array<OneD, NekDouble> &outarray);

            ///////////////////////////////////
            // Differentiation Methods
            //////////////////////////////////

            /** \brief Calculate the deritive of the physical points
            *
            *  \f$ \frac{\partial u}{\partial  x_1} =  \left . 
            *  \frac{2.0}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
            *  \right |_{\eta_2}\f$
            * 
            *  \f$ \frac{\partial u}{\partial  x_2} =  \left . 
            *  \frac{1+\eta_1}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
            *  \right |_{\eta_2}  + \left . \frac{\partial u}{\partial d\eta_2}
            *  \right |_{\eta_1}  \f$
            */
            void PhysDeriv(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            //-----------------------------
            // Evaluations Methods
            //-----------------------------

            /** \brief Backward tranform for triangular elements
            *
            *  \b Note: That 'q' (base[1]) runs fastest in this element 
            */
            void BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);
            void FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);

            /** \brief Single Point Evaluation */
            NekDouble PhysEvaluate(const ConstArray<OneD, NekDouble>& coords);

            void MapTo(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype, const int eid,
                const EdgeOrientation eorient, StdExpMap &Map);

            void MapTo_ModalFormat(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype,
                const int eid,
                const EdgeOrientation eorient,
                StdExpMap &Map);

            void WriteToFile(std::ofstream &outfile);
            void WriteCoeffsToFile(std::ofstream &outfile);

            int GetEdgeNcoeffs(const int i)
            {
                ASSERTL2((i >= 0) && (i <= 2), "edge id is out of range");

                if (i == 0)
                {
                    return GetBasisNumModes(0);
                }
                else
                {
                    return GetBasisNumModes(1);
                }

            }

            LibUtilities::BasisType GetEdgeBasisType(const int i)
            {
                ASSERTL2((i >= 0) && (i <= 2), "edge id is out of range");

                if (i == 0)
                {
                    return GetBasisType(0);
                }
                else
                {
                    return GetBasisType(1);
                }

            }

            void GetCoords(Array<OneD, NekDouble> &coords_0, 
                Array<OneD, NekDouble> &coords_1);

        protected:


            /** \brief Calculate the inner product of inarray with respect to
            *  the basis B=base0[p]*base1[pq] and put into outarray.
            *
            *  \f$ 
            *  \begin{array}{rcl}
            *  I_{pq} = (\phi^A_q \phi^B_{pq}, u) &=& 
            *  \sum_{i=0}^{nq_0}\sum_{j=0}^{nq_1}
            *  \phi^A_p(\eta_{0,i})\phi^B_{pq}(\eta_{1,j}) w^0_i w^1_j 
            *  u(\xi_{0,i} \xi_{1,j})\\
            *  & = & \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i})
            *  \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) \tilde{u}_{i,j} 
            *  \end{array}
            *  \f$ 
            *
            *  where
            *
            *  \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
            *
            *  which can be implemented as
            *
            *  \f$  f_{pj} = \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i}) 
            *  \tilde{u}_{i,j} 
            *  \rightarrow {\bf B_1 U}  \f$
            *  \f$  I_{pq} = \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) f_{pj} 
            *  \rightarrow {\bf B_2[p*skip] f[skip]}  \f$
            *
            *  \b Recall: \f$ \eta_{1} = \frac{2(1+\xi_1)}{(1-\xi_2)}-1, \, 
            *  \eta_2 = \xi_2\f$
            *
            *  \b Note: For the orthgonality of this expansion to be realised
            *  the 'q' ordering must run fastest in contrast to the Quad and
            *  Hex ordering where 'p' index runs fastest to be consistent with
            *  the quadrature ordering. 
            *
            *  In the triangular space the i (i.e. \f$\eta_1\f$ direction)
            *  ordering still runs fastest by convention.
            */
            void IProductWRTBase(const ConstArray<OneD, NekDouble>& base0, 
                const ConstArray<OneD, NekDouble>& base1,
                const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray);
        private:

            virtual int v_GetNverts()
            {
                return 3;
            }

            virtual int v_GetNedges()
            {
                return 3;
            }

            virtual int v_GetEdgeNcoeffs(const int i)
            {
                return GetEdgeNcoeffs(i);
            }

            /** \brief Virtual call to GenMassMatrix */
            virtual DNekMatSharedPtr v_GenMassMatrix() 
            {
                return StdExpansion::GenerateMassMatrix();
            }

            virtual DNekMatSharedPtr v_GenLaplacianMatrix() 
            {
                return GenLaplacianMatrix();
            }	    

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i)
            {
                return GetEdgeBasisType(i);
            }

            virtual ShapeType v_DetShapeType()
            {
                return DetShapeType();
            };

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
            {
                GetCoords(coords_0,coords_1);
            }

            virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray, outarray);
            }

            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                FillMode(mode, outarray);
            }

            virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0, out_d1);
            }

            virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0, 
                Array<OneD, NekDouble> &out_d1)
            {
                PhysDeriv(inarray, out_d0, out_d1);
            }

            /** \brief Virtual call to StdTriExp::BwdTrans */
            virtual void v_BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
            }

            /** \brief Virtual call to StdTriExp::FwdTrans */
            virtual void v_FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }

            virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
            {
                return PhysEvaluate(coords);
            }

            virtual void v_MapTo(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient,
                StdExpMap &Map)
            {
                MapTo(edge_ncoeffs, Btype, eid, eorient, Map);
            }

            virtual void v_MapTo_ModalFormat(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype,
                const int eid,
                const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                MapTo_ModalFormat(edge_ncoeffs, Btype, eid, eorient, Map);
            }

            virtual void v_WriteToFile(std::ofstream &outfile)
            {
                WriteToFile(outfile);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }
        };

    } //end of namespace
} //end of namespace
#endif //STDTRIEXP_H


/**
* $Log: StdTriExp.h,v $
* Revision 1.14  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.13  2007/05/31 19:13:12  pvos
* Updated NodalTriExp + LocalRegions/Project2D + some other modifications
*
* Revision 1.12  2007/05/15 05:18:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.11  2007/04/10 14:00:46  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.10  2007/04/06 08:44:43  sherwin
* Update to make 2D regions work at StdRegions level
*
* Revision 1.9  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.8  2007/03/31 08:18:07  bcarmo
* Changes in order to make the code compile
*
* Revision 1.7  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.6  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/07/02 17:16:19  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.2  2006/06/01 14:13:37  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:33  kirby
* *** empty log message ***
*
* Revision 1.39  2006/03/12 14:20:45  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.38  2006/03/05 23:17:53  sherwin
*
* Corrected to allow MMatrix1D and MMatrix2D to execute properly
*
* Revision 1.37  2006/03/05 22:11:03  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.36  2006/03/04 20:26:55  bnelson
* Added comments after #endif.
*
* Revision 1.35  2006/03/01 08:25:05  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.34  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/
