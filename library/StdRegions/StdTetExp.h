///////////////////////////////////////////////////////////////////////////////
//
// File StdTetExp.h
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDTETEXP_H
#define NEKTAR_LIB_STDREGIONS_STDTETEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdMatrixKey.h>


namespace Nektar
{
    namespace StdRegions
    {

        class StdTetExp: public StdExpansion3D
        {

        public:

            StdTetExp();
            /** \brief Constructor using BasisKey class for quadrature points and order definition 
            */
            StdTetExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature points and order definition where m_coeffs and m_phys are all set.
            */
            StdTetExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc, 
                      double *coeffs, double *phys);

            /** \brief Copy Constructor */
            StdTetExp(const StdTetExp &T);

            /** \brief Destructor */
            ~StdTetExp();

            /** \brief Return Shape of region, using  ShapeType enum list. 
             *  i.e. Tetrahedron
             */
            ShapeType DetShapeType() const
            {
                return eTetrahedron;
            }

          NekDouble Integral3D(const ConstArray<OneD, NekDouble>& inarray, const ConstArray<OneD, NekDouble>& w0,
                               const ConstArray<OneD, NekDouble>& w1, const ConstArray<OneD, NekDouble>&w2);
          NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray);
          void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray);
          void FillMode(const int mode, Array<OneD, NekDouble> &outarray);
          void PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
               Array<OneD, NekDouble> &out_d0, 
               Array<OneD, NekDouble> &out_d1,
               Array<OneD, NekDouble> &out_d2);
          void BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
               Array<OneD, NekDouble> &outarray);
          void FwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
               Array<OneD, NekDouble> &outarray);
          NekDouble PhysEvaluate(const ConstArray<OneD, NekDouble>& coords);
          NekDouble PhysEvaluate3D(const ConstArray<OneD, NekDouble>& coords);
          void WriteToFile(std::ofstream &outfile);
          void WriteCoeffsToFile(std::ofstream &outfile);
          void GetCoords(Array<OneD, NekDouble> &coords_0, 
               Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2);


            //----------------------------------
            // Generate Matrix Routine
            //----------------------------------

            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }


        //    void SetInvInfo(StdMatContainer *mat, MatrixType Mform);
        
        protected:

         void IProductWRTBase(const ConstArray<OneD, NekDouble>& base0, 
              const ConstArray<OneD, NekDouble>& base1, 
              const ConstArray<OneD, NekDouble>& base2, 
              const ConstArray<OneD, NekDouble>& inarray, 
              Array<OneD, NekDouble> & outarray);

          //  static StdMatrix s_elmtmats;

        private:

            virtual int v_GetNverts() const
            {
                return 4;
            }

            virtual int v_GetNedges() const
            {
                return 6;
            }

            virtual int v_GetNfaces() const
            {
                return 4;
            }

            virtual ShapeType v_DetShapeType() const
            {
                return DetShapeType();
            }

//             virtual void v_SetInvInfo(StdMatContainer *mat, MatrixType Mform)
//             {
//                 SetInvInfo(mat,Mform);
//             }


            // BEGIN ///

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                std::cout << "Implement me" << std::endl;
                return -1;
                //return GetEdgeNcoeffs(i);
            }


            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey) 
            {
                return GenMatrix(mkey);
            }


            virtual void v_GenMassMatrix(Array<OneD, NekDouble> & outarray) 
            {
                std::cout << "Implement me" << std::endl;
                return;
            } // TODO: Implement
            
            virtual void v_GenLapMatrix (Array<OneD, NekDouble> & outarray) 
            {
                std::cout << "Implement me" << std::endl;
                return;
            } // TODO: Implement
            
            virtual DNekMatSharedPtr v_GetMassMatrix() 
            {
                std::cout << "Implement me" << std::endl;
                int foo = 0;
                return DNekMatSharedPtr();
            } // TODO: Implement
            
            virtual DNekMatSharedPtr v_GetLapMatrix() 
            {
                std::cout << "Implement me" << std::endl;
                int foo = 0;
                return DNekMatSharedPtr();
            }  // TODO: Implement


            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                std::cout << "Implement me" << std::endl;
                return LibUtilities::eNoBasisType;
                //return GetEdgeBasisType(i);
            }


            virtual void v_GetCoords(
                Array<OneD, NekDouble> &coords_x,
                Array<OneD, NekDouble> &coords_y,
                Array<OneD, NekDouble> &coords_z)
            {
                GetCoords(coords_x, coords_y, coords_z);
            }

            virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                return Integral(inarray);
            }
            
            virtual NekDouble v_Evaluate(const ConstArray<OneD, NekDouble> &coords) {
                std::cout << "Implement me" << std::endl;
                return -1; // TODO: Implement
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

            virtual void v_PhysDeriv(
                const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &out_dx,
                Array<OneD, NekDouble> &out_dy,
                Array<OneD, NekDouble> &out_dz )
            {
                std::cout << "Implement me" << std::endl;
                PhysDeriv( inarray, out_dx, out_dy, out_dz );
            }

            virtual void v_StdPhysDeriv(
                const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &out_dx, 
                Array<OneD, NekDouble> &out_dy,
                Array<OneD, NekDouble> &out_dz )
            {
                std::cout << "Implement me" << std::endl;
                PhysDeriv( inarray, out_dx, out_dy, out_dz );
            }

            /** \brief Virtual call to StdTetExp::BwdTrans */
            virtual void v_BwdTrans(
                const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray )
            {
                BwdTrans(inarray,outarray);
            }

            /** \brief Virtual call to StdTetExp::FwdTrans */
            virtual void v_FwdTrans(
                const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray )
            {
                FwdTrans(inarray,outarray);
            }

            virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
            {
                std::cout << "Implement me" << std::endl;
                return PhysEvaluate(coords);
            }

            virtual void v_MapTo(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient,
                StdExpMap &Map)
            {
                std::cout << "Implement me" << std::endl;
                return;
                //MapTo(edge_ncoeffs, Btype, eid, eorient, Map);
            }

            virtual void v_MapTo_ModalFormat(const int edge_ncoeffs,
                const LibUtilities::BasisType Btype,
                const int eid,
                const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                std::cout << "Implement me" << std::endl;
                return;
                //MapTo_ModalFormat(edge_ncoeffs, Btype, eid, eorient, Map);
            }

            virtual void v_WriteToFile(std::ofstream &outfile)
            {
                std::cout << "Implement me" << std::endl;
                WriteToFile(outfile);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                std::cout << "Implement me" << std::endl;
                WriteCoeffsToFile(outfile);
            }

        };    
        typedef boost::shared_ptr<StdTetExp> StdTetExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDTETEXP_H

/**
 * $Log: StdTetExp.h,v $
 * Revision 1.11  2008/01/08 22:48:41  ehan
 * Fixed the call signature of a shadowed virtual function: Added a const qualifier to the passed parameter StdMatrixKey in the virtual function v_GenMatrix().  This enables Nektar to generate the correct standard mass matrix at initialization time.
 *
 * Revision 1.10  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.9  2007/10/28 18:32:34  bnelson
 * Fixed visual c++ compile error.
 *
 * Revision 1.8  2007/10/15 20:40:07  ehan
 * Completed Basis, Backward, and Forward transformation
 *
 * Revision 1.7  2007/07/20 02:16:55  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.6  2007/07/10 21:05:17  kirby
 * even more fixes
 *
 * Revision 1.5  2007/07/10 20:41:52  kirby
 * more fixes
 *
 * Revision 1.4  2007/01/17 16:05:41  pvos
 * updated doxygen documentation
 *
 * Revision 1.3  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.2  2006/07/02 17:16:19  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.24  2006/03/12 21:59:48  sherwin
 *
 * compiling version of LocalRegions
 *
 * Revision 1.23  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.22  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.21  2006/03/01 08:25:05  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

