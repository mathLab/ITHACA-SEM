///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.h
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
// Description: Hex routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGSION_STDHEXEXP_H
#define NEKTAR_LIBS_STDREGSION_STDHEXEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>


namespace Nektar
{
    namespace StdRegions
    {

        class StdHexExp: public StdExpansion3D
        {

        public:
        
            StdHexExp();

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition 
            */
            StdHexExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition where m_coeffs and m_phys 
            *  are all set. 
            */
            StdHexExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc,
                double *coeffs, double *phys);

            /** \brief Copy Constructor */
            StdHexExp(const StdHexExp &T);

            /** \brief Destructor */
            ~StdHexExp();
            

            /** \brief Return Shape of region, using  ShapeType enum list. 
            *  i.e. Hexahedron
            */
      
            ShapeType DetShapeType() const
            {
                return eHexahedron;
            }

            /** \brief Fill outarray with mode \a mode of expansion
            *
            *    Note for hexahedral expansions _base[0] (i.e. p)  modes run 
            *  fastest
            */
            //void FillMode(int mode, double *array);
            void FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

            NekDouble Integral3D(const ConstArray<OneD, NekDouble>& inarray, 
                                const ConstArray<OneD, NekDouble>& wx,
                                const ConstArray<OneD, NekDouble>& wy, 
                                const ConstArray<OneD, NekDouble>& wz);
            NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray);
            
            void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> & outarray);

            void IProductWRTBase(const ConstArray<OneD, NekDouble>& bx, 
                                 const ConstArray<OneD, NekDouble>& by, 
                                 const ConstArray<OneD, NekDouble>& bz, 
                                 const ConstArray<OneD, NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray );

            //----------------------------------
            // Local Matrix Routines
            //----------------------------------
            //void GenMassMatrix(double * outarray);
              DNekMatSharedPtr GenMatrixHex(MatrixType mtype);
           // DNekMatSharedPtr GenMatrix(MatrixType mtype);

             void GenLapMatrix(double * outarray);

//             StdMatContainer * GetMassMatrix();
//             StdMatContainer * GetLapMatrix();

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the deritive of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
           // void Deriv(double * outarray_d1, double *outarray_d2, double *outarray_d3);

            /** \brief Calculate the deritive of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
//             void Deriv(const double *inarray, double * outarray_d1,
//                 double *outarray_d2, double * outarray_d3);
                
         void PhysDeriv( Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2);

         void PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                   Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2);
                                   
          void GetCoords(Array<OneD, NekDouble> &coords_0, 
               Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2);

            //----------------------------
            // Evaluations Methods
            //---------------------------

            //void BwdTrans(double * outarray);
            void BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            void FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);
            NekDouble PhysEvaluate(ConstArray<OneD, NekDouble>& coords);


            //----------------------------------
            // Local Matrix Routines
            //----------------------------------

            DNekMatSharedPtr GenMassMatrix();

            DNekMatSharedPtr GenLaplacianMatrix();

            DNekMatSharedPtr GenLaplacianMatrix(const int i, const int j);

            DNekMatSharedPtr GenWeakDerivMatrix(const int i);

            DNekMatSharedPtr GenNBasisTransMatrix();

            DNekMatSharedPtr GenBwdTransMatrix();


            DNekMatSharedPtr GenMatrix(MatrixType mtype)
            {
                return StdExpansion::CreateGeneralMatrix(mtype);
            }

        protected:
                
            void IProductWRTBase(const ConstArray<OneD, NekDouble>&base0, 
                                 const ConstArray<OneD, NekDouble>& base1, 
                                 const ConstArray<OneD, NekDouble>& base2,
                                 const ConstArray<OneD, NekDouble>& inarray,
                                 Array<OneD, NekDouble>& outarray, int coll_check);

        private:

            virtual int v_GetNverts() const
            {
                return 8;
            }

            virtual int v_GetNedges() const
            {
                return 12;
            }

            virtual int v_GetNfaces() const
            {
                return 6;
            }

            virtual ShapeType v_DetShapeType() const
            {
                return DetShapeType();
            };
            
            virtual DNekMatSharedPtr v_GenMatrix(MatrixType mtype) 
            {
                return GenMatrix(mtype);
            }
           
            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                return FillMode(mode, outarray);
            }

            virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                return Integral(inarray);
            }
            
            virtual void v_GetCoords(
                Array<OneD, NekDouble> &coords_x,
                Array<OneD, NekDouble> &coords_y,
                Array<OneD, NekDouble> &coords_z)
            {
                GetCoords(coords_x, coords_y, coords_z);
            }
            
            virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray, outarray);
            }

            /** \brief Virtual call to GenMassMatrix */

            virtual void v_PhysDeriv( Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
            {
                    PhysDeriv(out_d0, out_d1, out_d2);
            }
            virtual void v_StdPhysDeriv( Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
            {
                                   
                    PhysDeriv(out_d0, out_d1, out_d2);                
            }

            virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                   Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
            {
                    PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }                                  
            virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                   Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
            {
                    PhysDeriv(inarray, out_d0, out_d1, out_d2);
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
            
             virtual void v_GenMassMatrix(Array<OneD, NekDouble> & outarray) 
            {
                 std::cout << "Implement me" << std::endl;
                 return;
            } 
            
            virtual void v_GenLapMatrix (Array<OneD, NekDouble> & outarray)
            {
                std::cout << "Implement me" << std::endl;
                return;
            }
            
            virtual DNekMatSharedPtr v_GetMassMatrix() 
            {
                std::cout << "Implement me" << std::endl;
                int foo = 0;
            } 
            
            virtual DNekMatSharedPtr v_GetLapMatrix()
            {
                std::cout << "Implement me" << std::endl;
                int foo = 0;
            }  

        };

    } //end of namespace
} //end of namespace

#endif //STDHEXEXP_H

/**
* $Log: StdHexExp.h,v $
* Revision 1.10  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/07/10 21:05:16  kirby
* even more fixes
*
* Revision 1.7  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.6  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
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
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.30  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.29  2006/03/06 17:12:45  sherwin
*
* Updated to properly execute all current StdRegions Demos.
*
* Revision 1.28  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.27  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.26  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
*
**/



