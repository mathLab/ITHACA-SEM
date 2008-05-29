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
            
            virtual bool v_IsBoundaryInteriorExpansion()
            {
                bool returnval = false;
                
                if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
                {
                    if(m_base[1]->GetBasisType() == LibUtilities::eModified_B)
                    {
                       if(m_base[2]->GetBasisType() == LibUtilities::eModified_C)
                        {
                          returnval = true;
                        }
                    }
                }
                
                return returnval;
            }


          NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray, const Array<OneD, const NekDouble>& w0,
                               const Array<OneD, const NekDouble>& w1, const Array<OneD, const NekDouble>&w2);
          NekDouble Integral(const Array<OneD, const NekDouble>& inarray);
          void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, Array<OneD, NekDouble> &outarray);
          void FillMode(const int mode, Array<OneD, NekDouble> &outarray);
          void PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, NekDouble> &out_d0,
                         Array<OneD, NekDouble> &out_d1,
                         Array<OneD, NekDouble> &out_d2);
                         
          /** \brief Backward tranform for triangular elements
          *
          *  \b Note: That 'r' (base[2]) runs fastest in this element
          */
          void BwdTrans(const Array<OneD, const NekDouble>& inarray, 
                        Array<OneD, NekDouble> &outarray);
          void FwdTrans(const Array<OneD, const NekDouble>& inarray, 
                        Array<OneD, NekDouble> &outarray);
                        
          /** \brief Single Point Evaluation */
          NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coords);
          NekDouble PhysEvaluate3D(const Array<OneD, const NekDouble>& coords);
          
          void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);
          void WriteCoeffsToFile(std::ofstream &outfile);
          void GetCoords(Array<OneD, NekDouble> &coords_0, 
                         Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2);
                         
          void MapTo(const int edge_ncoeffs,
                     const LibUtilities::BasisType Btype, const int eid,
                     const EdgeOrientation eorient, StdExpMap &Map);
          void MapTo_ModalFormat(const int edge_ncoeffs,
                                 const LibUtilities::BasisType Btype,
                                 const int eid,
                                 const EdgeOrientation eorient,
                                 StdExpMap &Map);


            //----------------------------------
            // Generate Matrix Routine
            //----------------------------------

            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }

            int GetEdgeNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0) && (i <= 5), "edge id is out of range");

                if (i == 0)
                {
                    return GetBasisNumModes(0);
                }
                else if((i==1)||(i==2))
                {
                    return GetBasisNumModes(1);
                } else
                {
                    return GetBasisNumModes(2);
                }

            }

            LibUtilities::BasisType GetEdgeBasisType(const int i) const
            {
                ASSERTL2((i >= 0) && (i <= 5), "edge id is out of range");

                if (i == 0)
                {
                    return GetBasisType(0);
                }
                else if((i==1)||(i==2))
                {
                    return GetBasisType(1);
                }
                else
                {
                    return GetBasisType(2);
                }
            }
        
        protected:

         void IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                              const Array<OneD, const NekDouble>& base1,
                              const Array<OneD, const NekDouble>& base2,
                              const Array<OneD, const NekDouble>& inarray,
                              Array<OneD, NekDouble> & outarray);

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


            // BEGIN ///

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey) 
            {
                return GenMatrix(mkey);
            }

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }            

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            virtual void v_GetCoords( Array<OneD, NekDouble> &coords_x,
                                      Array<OneD, NekDouble> &coords_y,
                                      Array<OneD, NekDouble> &coords_z)
            {
                GetCoords(coords_x, coords_y, coords_z);
            }

            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray )
            {
                return Integral(inarray);
            }
            
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray, outarray);
            }

            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                FillMode(mode, outarray);
            }

            virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_dx,
                                     Array<OneD, NekDouble> &out_dy,
                                     Array<OneD, NekDouble> &out_dz )
            {

                PhysDeriv( inarray, out_dx, out_dy, out_dz );
            }

            virtual void v_StdPhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &out_dx,
                                        Array<OneD, NekDouble> &out_dy,
                                        Array<OneD, NekDouble> &out_dz )
            {

                PhysDeriv( inarray, out_dx, out_dy, out_dz );
            }

            /** \brief Virtual call to StdTetExp::BwdTrans */
            virtual void v_BwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray )
            {
                BwdTrans(inarray,outarray);
            }

            /** \brief Virtual call to StdTetExp::FwdTrans */
            virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray )
            {
                FwdTrans(inarray,outarray);
            }

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
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

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true)
            {                
                WriteToFile(outfile,format,dumpVar);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {               
                WriteCoeffsToFile(outfile);
            }

        };
        
        typedef boost::shared_ptr<StdTetExp> StdTetExpSharedPtr;
        
    } //end of namespace
} //end of namespace

#endif //STDTETEXP_H

/**
 * $Log: StdTetExp.h,v $
 * Revision 1.16  2008/05/15 22:42:29  ehan
 * Added WriteToFile() function and its virtual function
 *
 * Revision 1.15  2008/05/15 04:15:52  ehan
 * Added virtual function v_CreatStdMatrix()
 *
 * Revision 1.14  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.13  2008/03/25 08:40:20  ehan
 * Added GetEdgeNcoeffs() and GetEdgeBasisType().
 *
 * Revision 1.12  2008/03/17 10:37:58  pvos
 * Clean up of the code
 *
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

