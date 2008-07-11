///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionContField2D.h
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
// Description: Advection Field definition in tow-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_ADVECTIONCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_ADVECTIONCONTFIELD2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContField2D.h>

#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>


namespace Nektar
{
    namespace MultiRegions
    {        
        class AdvectionContField2D
        {
        public:           
          
            /**
             * \brief 
             */ 
            AdvectionContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs);

            void SetInitialConditions(SpatialDomains::BoundaryConditions &bcs, int initialtime = 0);

            void LinearAdvectionOperation(const Array<OneD,const NekDouble>& a, const Array<OneD,const NekDouble>& b);

            void   FwdTrans(const AdvectionContField2D &In);

            void   BwdTrans(const AdvectionContField2D &In);


         
        protected:
            Array<OneD, ExpListSharedPtr> m_fields;


        private:  
        };

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_ADVECTIONCONTFIELD2D_H
