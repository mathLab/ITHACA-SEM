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

#ifndef STDHEXEXP_H
#define STDHEXEXP_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdMatrix.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdHexExp: public StdExpansion3D
        {

        public:

            /** \brief Constructor using BasisKey class for quadrature
            points and order definition */
            StdHexExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature
            points and order definition where _coeffs and _phys are all
            set. */
            StdHexExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc,
                double *coeffs, double *phys);

            /// Copy Constructor
            StdHexExp(const StdHexExp &T);

            /// Destructor
            ~StdHexExp();

            /// Return Shape of region, using  ShapeType enum list. i.e. Hexahedron
            ShapeType DetShapeType()
            {
                return eHexahedron;
            };

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

            void Deriv(double * outarray_d1, double *outarray_d2, double *outarray_d3);
            void Deriv(const double *inarray, double * outarray_d1,
                double *outarray_d2, double * outarray_d3);

            //----------------------------
            // Evaluations Methods
            //---------------------------

            void BwdTrans(double * outarray);
            void FwdTrans(const double * inarray);
            double Evaluate(const double * coords);

        protected:

            static StdMatrix s_elmtmats;

            void IProductWRTBase(const double *base0, const double *base1,
                const double *base2, const double *inarray,
                double *outarray, int coll_check);

        private:

            virtual ShapeType v_DetShapeType()
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


            virtual void v_GenMassMatrix(double * outarray)
            {
                GenMassMatrix(outarray);
            }

            virtual void v_GenLapMatrix(double * outarray)
            {
                GenLapMatrix(outarray);
            }

            virtual StdMatContainer *v_GetMassMatrix()
            {
                return GetMassMatrix();
            }

            virtual StdMatContainer * v_GetLapMatrix()
            {
                return GetLapMatrix();
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

        };

    } //end of namespace
} //end of namespace

#endif //STDHEXEXP_H

/**
* $Log: StdHexExp.h,v $
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



