///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.h
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
// Description: Field definition in one-dimension for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/BoundaryConditions.h>
#include <MultiRegions/LocalToGlobalBndryMap1D.h>
#include <MultiRegions/GlobalLinSys.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField1D: public ExpList1D
        {
        public:
            DisContField1D();

            DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
                SpatialDomains::BoundaryConditions &bcs, const int bc_loc = 0);

            DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
                SpatialDomains::BoundaryConditions &bcs, 
                const std::string variable);

            DisContField1D(const DisContField1D &In);

            ~DisContField1D();

            void HelmSolve(DisContField1D &Fce, NekDouble lambda);

            GlobalLinSysSharedPtr GetGlobalBndLinSys(const GlobalLinSysKey &mkey);

        protected:

        private:
            Array<OneD,LocalRegions::PointExpSharedPtr>        m_bndConstraint;
            Array<OneD,SpatialDomains::BoundaryConditionType>  m_bndTypes;
	    LocalToGlobalBndryMapSharedPtr                     m_lambdaMap;
            GlobalLinSysMapShPtr                               m_globalBndMat;
            
            void GenerateField1D(SpatialDomains::BoundaryConditions &bcs,  
                                 const std::string variable);
        };

        typedef boost::shared_ptr<DisContField1D>   DisContField1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD1D_H
