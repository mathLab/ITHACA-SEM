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

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition 
            */
            StdHexExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature
            *  points and order definition where m_coeffs and m_phys 
            *  are all set. 
            */
            StdHexExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc,
                double *coeffs, double *phys);

            /** \brief Copy Constructor */
            StdHexExp(const StdHexExp &T);

            /** \brief Destructor */
            ~StdHexExp();

            /** \brief Return Shape of region, using  ShapeType enum list. 
            *  i.e. Hexahedron
            */
            ShapeType DetShapeType()
            {
                return eHexahedron;
            };

            /** \brief Fill outarray with mode \a mode of expansion
            *
            *    Note for hexahedral expansions _base[0] (i.e. p)  modes run 
            *  fastest
            */
            void FillMode(int mode, double *array);

            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

            double Integral(const double *inarray);
            void IProductWRTBase(const double * inarray, double * outarray);

            //----------------------------------
            // Local Matrix Routines
            //----------------------------------

            void GenMassMatrix(double * outarray);

            void GenLapMatrix(double * outarray);

            StdMatContainer * GetMassMatrix();
            StdMatContainer * GetLapMatrix();

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the deritive of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
            void Deriv(double * outarray_d1, double *outarray_d2, double *outarray_d3);

            /** \brief Calculate the deritive of the physical points 
            *
            *  For quadrilateral region can use the Tensor_Deriv function
            *  defined under StdExpansion.
            */
            void Deriv(const double *inarray, double * outarray_d1,
                double *outarray_d2, double * outarray_d3);

            //----------------------------
            // Evaluations Methods
            //---------------------------

            void BwdTrans(double * outarray);
            void FwdTrans(const double * inarray);
            double Evaluate(const double * coords);

            // Matrix related methods 
            void SetInvInfo(StdMatContainer *mat, MatrixType Mform);

            //----------------------------------
            // Local Matrix Routines
            //----------------------------------

            DNekMatSharedPtr GenMassMatrix();

            DNekMatSharedPtr GenLaplacianMatrix();

            DNekMatSharedPtr GenLaplacianMatrix(const int i, const int j);

            DNekMatSharedPtr GenWeakDerivMatrix(const int i);

            DNekMatSharedPtr GenNBasisTransMatrix();

            DNekMatSharedPtr GenBwdTransMatrix();

        protected:

            static StdMatrix s_elmtmats;

            void IProductWRTBase(const double *base0, const double *base1,
                const double *base2, const double *inarray,
                double *outarray, int coll_check);

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

            virtual void v_FillMode(const int mode, double *array)
            {
                FillMode(mode,array);
            }

            virtual double v_Integral(const double *inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const double * inarray, double * outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            /** \brief Virtual call to GenMassMatrix */
            virtual DNekMatSharedPtr v_GenMassMatrix() 
            {
                return GenMassMatrix();
            }

            virtual DNekMatSharedPtr v_GenLaplacianMatrix() 
            {
                return GenLaplacianMatrix();
            }

            virtual DNekMatSharedPtr v_GenLaplacianMatrix(const int i, const int j) 
            {
                return GenLaplacianMatrix(i,j);
            }

            virtual DNekMatSharedPtr v_GenWeakDerivMatrix(const int i) 
            {
                return GenWeakDerivMatrix(i);
            }

            virtual DNekMatSharedPtr v_GenNBasisTransMatrix() 
            {
                return GenNBasisTransMatrix();
            }

            virtual DNekMatSharedPtr v_GenBwdTransMatrix() 
            {
                return GenBwdTransMatrix();
            }

            virtual void v_Deriv(double * outarray_d1, double *outarray_d2,
                double *outarray_d3)
            {
                Deriv(outarray_d1, outarray_d2, outarray_d3);
            }

            virtual void v_StdDeriv(double * outarray_d1, double *outarray_d2,
                double *outarray_d3)
            {
                Deriv(outarray_d1, outarray_d2, outarray_d3);
            }

            virtual void v_Deriv(const double *inarray, double * outarray_d1,
                double *outarray_d2, double *outarray_d3)
            {
                Deriv(inarray, outarray_d1, outarray_d2, outarray_d3);
            }

            virtual void v_StdDeriv(const double *inarray, double * outarray_d1,
                double *outarray_d2, double *outarray_d3)
            {
                Deriv(inarray, outarray_d1, outarray_d2, outarray_d3);
            }


            virtual void v_BwdTrans(double * outarray)
            {
                BwdTrans(outarray);
            }

            virtual void v_FwdTrans(const double * inarray)
            {
                FwdTrans(inarray);
            }

            virtual double v_Evaluate(const double * coords)
            {
                return Evaluate(coords);
            }

            virtual void v_SetInvInfo(StdMatContainer *mat, MatrixType Mform)
            {
                SetInvInfo(mat,Mform);
            }

        };

    } //end of namespace
} //end of namespace

#endif //STDHEXEXP_H

/**
* $Log: StdHexExp.h,v $
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



