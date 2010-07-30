///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3DHomogeneous1D.h
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
// Description: A 2D field which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FIELD3DHOMO1D_H
#define FIELD3DHOMO1D_H

#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <SpatialDomains/MeshGraph2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        class ExpList3DHomogeneous1D;

        /// Shared pointer to an ExpList3DHomogeneous1D object.
        typedef boost::shared_ptr<ExpList3DHomogeneous1D>      ExpList3DHomogeneous1DSharedPtr;
        /// Vector of pointers to ExpList3DHomogeneous1D objects.
        typedef std::vector< ExpList3DHomogeneous1DSharedPtr > ExpList3DHomogeneous1DVector;
        /// Iterator for the vector of ExpList3DHomogeneous1D pointers.
        typedef std::vector< ExpList3DHomogeneous1DSharedPtr >::iterator ExpList3DHomogeneous1DVectorIter;

        /// Abstraction of a two-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpList3DHomogeneous1D: public ExpList
        {
        public:
            /// Default constructor.
            ExpList3DHomogeneous1D();

            /// Sets up a list of local expansions based on an input mesh.
            ExpList3DHomogeneous1D(const int nzplanes,  const NekDouble lz, 
                                   SpatialDomains::MeshGraph2D &graph2D);

            /// Copy constructor.
            ExpList3DHomogeneous1D(const ExpList3DHomogeneous1D &In);

            /// Destructor.
            ~ExpList3DHomogeneous1D();

        protected:

            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            void      SetCoeffPhys(void);

            int       m_nzplanes; ///< Number of planes of data in expansion in z-direction
            NekDouble m_lz;       ///< Width of homogeneous direction
            Array<OneD, ExpList2DSharedPtr> m_planes; 

            //  virtual functions
            virtual void v_FwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                                     Array<OneD, NekDouble> &coord_1,
                                     Array<OneD, NekDouble> &coord_2);

        private:

        };
        

    } //end of namespace
} //end of namespace

#endif//FIELD3DHOMO1D_H

/**
* $Log: v $
*
**/

