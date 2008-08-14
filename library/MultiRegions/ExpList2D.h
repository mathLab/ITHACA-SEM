///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2D.h
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
// Description: Expansion list 2D header definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST2D_H
#define EXPLIST2D_H

#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>

namespace Nektar
{
    namespace MultiRegions
    {      
        /**
         * \brief This class is the abstraction of a two-dimensional multi-elemental 
         * expansions which is merely a collection of local expansions.
         * 
         * This multi-elemental expansion, which does not exhibit any coupling between the 
         * expansion on the separate elements, can be formulated as,
         * \f[u^{\delta}(\boldsymbol{x}_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and \f$N^{e}_m\f$ is the 
         * local elemental number of expansion modes.
         * This class inherits all its variables and member functions from the base class 
         * #ExpList.
         */  
    class ExpList2D:
        public ExpList
        {
        public:
            /**
             * \brief The default constructor.  
             */  
            ExpList2D(); 
            
            /**
             * \brief The copy constructor.
             */  
            ExpList2D(const ExpList2D &In);   
            
            /**
             * \brief 
             */  
            ExpList2D(const LibUtilities::BasisKey &TriBa, 
                      const LibUtilities::BasisKey &TriBb, 
                      const LibUtilities::BasisKey &QuadBa, 
                      const LibUtilities::BasisKey &QuadBb, 
                      const SpatialDomains::MeshGraph2D &graph2D,
                      const LibUtilities::PointsType 
                      TriNb = LibUtilities::SIZE_PointsType);

            /**
             * \brief This constructor sets up a list of local expansions based on an input 
             * mesh.
             * 
             * Given a mesh \a graph2D, containing information about the domain and the 
             * spectral/hp element expansion, this constructor fills the list of local 
             * expansions \texttt{m_exp} with the proper expansions, calculates the total 
             * number of quadrature points \f$\boldsymbol{x}_i\f$ and local expansion 
             * coefficients \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs 
             * and #m_phys.
             *
             * \param graph2D A mesh, containing information about the domain and the 
             * spectral/hp element expansion.
             */  
            ExpList2D(SpatialDomains::MeshGraph2D &graph2D);
            
            /**
             * \brief The default destructor.
             */  
            ~ExpList2D();
            
            /**
             * \brief 
             */  
            void   PhysDeriv  (ExpList &S0,
                               ExpList &S1, 
                               ExpList &S2 = NullExpList)
            {
                ExpList::PhysDeriv(S0,S1,S2);
            }
            
        protected:
            void SetBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                               SpatialDomains::BoundaryConditions &bcs, 
                                               const std::string variable,
                                               Array<OneD, ExpList1DSharedPtr> &bndCondExpansions,
                                               Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions);

            void EvaluateBoundaryConditions(const NekDouble time,
                                            Array<OneD, ExpList1DSharedPtr> &bndCondExpansions,
                                            Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions);
            
            void GetPeriodicEdges(SpatialDomains::MeshGraph2D &graph2D,
                                  SpatialDomains::BoundaryConditions &bcs, 
                                  const std::string variable,
                                  map<int,int>& periodicVertices,
                                  map<int,int>& periodicEdges);
        private:
            
        };
        
        typedef boost::shared_ptr<ExpList2D>      ExpList2DSharedPtr;
        typedef std::vector< ExpList2DSharedPtr > ExpList2DVector;
        typedef std::vector< ExpList2DSharedPtr >::iterator ExpList2DVectorIter;
    } //end of namespace
} //end of namespace

#endif//EXPLIST2D_H

/**
* $Log: ExpList2D.h,v $
* Revision 1.13  2008/06/05 15:06:58  pvos
* Added documentation
*
* Revision 1.12  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.11  2007/07/22 23:04:21  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.10  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

