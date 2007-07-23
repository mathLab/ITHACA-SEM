///////////////////////////////////////////////////////////////////////////////
//
// File ContSolnField1D.h
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
// Description: Field definition in one-dimension
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTSOLNFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTSOLNFIELD1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>

#include <LocalRegions/PointExp.h>
#include <SpatialDomains/BoundaryConditions.h>


namespace Nektar
{
    namespace MultiRegions
    {

	class ContField1D:
	    public ContExpList1D
	    {
	    public:
		ContField1D();
                ContField1D(const LibUtilities::BasisKey &Ba, 
                            const SpatialDomains::Composite &cmps,
                            SpatialDomains::BoundaryConditions &bcs,
                            const int bc_loc = 0);
                ContField1D(const ContField1D &In);
		~ContField1D();
		
                void SetBoundaryCondition(const int loc, const NekDouble value)
                {
                    m_bndConstraint[loc]->SetValue(value);
                }

                void FwdTrans (const ExpList &In);
                void HelmSolve(const ExpList &In, NekDouble lambda);

	    protected:
		
	    private:
		LocalRegions::PointExpVector                         m_bndConstraint;
                std::vector<SpatialDomains::BoundaryConditionType>   m_bndTypes;

                GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);
                void GlobalSolve(const GlobalLinSysKey &key, const ExpList &Rhs);

	    };
        typedef boost::shared_ptr<ContField1D>      ContField1DSharedPtr;
	
    } //end of namespace
} //end of namespace
  
#endif // MULTIERGIONS_CONTSOLNFIELD1D_H
