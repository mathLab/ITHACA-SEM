///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion3D.h
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
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP3D_H
#define STDEXP3D_H

#include <StdRegions/StdExpansion.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion3D: public StdExpansion
        {

        public:
            StdExpansion3D();
            StdExpansion3D(const BasisKey &Ba, const BasisKey &Bb,
                const BasisKey &Bc, int numcoeffs,
                double *coeffs, double *phys, bool spaceowner);
            StdExpansion3D(const StdExpansion3D &T);
            ~StdExpansion3D();

            // Differentiation
            void TensorDeriv(const double *inarray, double *outarray_d1,
                double *outarray_d2, double *outarray_d3);

            void TensorDeriv(double *outarray_d1, double *outarray_d2,
                double *outarray_d3);

            /** \brief Evaluate a function at points coords which is assumed
            to be in local collapsed coordinate format. The function is
            assumed to be in physical space */
            double PhysEvaluate(const double *coords);

            void  Deriv(double *outarray_d0, double *outarray_d1, 
                double *outarray_d2)
            {
                v_Deriv(outarray_d0, outarray_d1, outarray_d2);
            }

            void  StdDeriv(double *outarray_d0, double *outarray_d1, 
                double *outarray_d2)
            {
                v_StdDeriv(outarray_d0, outarray_d1, outarray_d2);
            }

            void Deriv(const double *inarray, double *outarray_d0, 
                double *outarray_d1, double *outarray_d2)
            {
                v_Deriv(inarray, outarray_d0, outarray_d1, outarray_d2);
            }

            void StdDeriv(const double *inarray, double *outarray_d0,
                double *outarray_d1,   double *outarray_d2)
            {
                v_StdDeriv(inarray, outarray_d0, outarray_d1, outarray_d2);
            }


        protected:

        private:

            // Virtual Functions ----------------------------------------

            virtual void v_GenMassMatrix(double * outarray)   = 0;
            virtual void v_GenLapMatrix (double * outarray)   = 0;
            virtual ShapeType v_DetShapeType()                = 0;

            virtual StdMatContainer *v_GetMassMatrix()        = 0;
            virtual StdMatContainer *v_GetLapMatrix()         = 0;

            virtual int v_get_nodalpoints(const double *x, const double *y)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return 0;
            }

            virtual void v_GenNBasisTransMatrix(double * outarray)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
            }

            virtual StdMatContainer *v_GetNBasisTransMatrix()
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return NULL;
            }

            virtual void   v_BwdTrans (double *outarray)      = 0;
            virtual void   v_FwdTrans (const double *inarray) = 0;

            virtual double v_Integral(const double *inarray ) = 0;
            virtual double v_Evaluate(const double * coords) = 0;

            virtual void   v_Deriv(double *outarray_d0, double *outarray_d1,
                double *outarray_d2) = 0;
            virtual void   v_StdDeriv(double *outarray_d0, double *outarray_d1,
                double *outarray_d2) = 0;
            virtual void   v_Deriv(const double *inarray, double *outarray_d0,
                double *outarray_d1, double *outarray_d2) = 0;
            virtual void   v_StdDeriv(const double *inarray, double *outarray_d0,
                double *outarray_d1, double *outarray_d2) = 0;

        };

    } //end of namespace
} //end of namespace

#endif //STDEXP3D_H

/**
* $Log: StdExpansion3D.h,v $
* Revision 1.12  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.11  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.10  2006/03/06 17:12:45  sherwin
*
* Updated to properly execute all current StdRegions Demos.
*
* Revision 1.9  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.8  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.7  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
**/

