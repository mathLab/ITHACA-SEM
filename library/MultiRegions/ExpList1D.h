///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.h
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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H
#define NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

#include <vector>
#include <fstream>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/GenSegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LibUtilities/Kernel/kernel.h>

namespace Nektar
{
    namespace MultiRegions
    {     
        /**
         * \brief This class is the abstraction of a one-dimensional multi-elemental 
         * expansions which is merely a collection of local expansions.
         * 
         * This multi-elemental expansion, which does not exhibit any coupling between the 
         * expansion on the separate elements, can be formulated as,
         * \f[u^{\delta}(x_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(x_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and \f$N^{e}_m\f$ is the 
         * local elemental number of expansion modes.
         * This class inherits all its variables and member functions from the base class 
         * #ExpList.
         */ 


        class ExpList1D: 
        public ExpList
        {
        public:

            /**
             * \brief The default constructor.  
             */  
            ExpList1D();
            
            /**
             * \brief The copy constructor.
             */  
            ExpList1D(const ExpList1D &In);   
            
            /**
             * \brief 
             */  
            ExpList1D(const LibUtilities::BasisKey &Ba, 
                      const SpatialDomains::MeshGraph1D &graph1D);

            /**
             * \brief This constructor sets up a list of local expansions based on an input 
             * mesh.
             * 
             * Given a mesh \a graph1D, containing information about the domain and the 
             * spectral/hp element expansion, this constructor fills the list of local 
             * expansions \texttt{m_exp} with the proper expansions, calculates the total 
             * number of quadrature points \f$x_i\f$ and local expansion 
             * coefficients \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs 
             * and #m_phys.
             *
             * \param graph1D A mesh, containing information about the domain and the 
             * spectral/hp element expansion.
             */  
            ExpList1D(SpatialDomains::MeshGraph1D &graph1D);
            
            /**
             * \brief 
             */  
            ExpList1D(const SpatialDomains::CompositeVector &domain, 
                      SpatialDomains::MeshGraph2D &graph2D); 
            
            /**
             * \brief 
             */  
            ExpList1D(const Array<OneD,const boost::shared_ptr<ExpList1D> > &bndConstraint,  const Array<OneD,const SpatialDomains::BoundaryConditionType>  &bndTypes, const StdRegions::StdExpansionVector &locexp, SpatialDomains::MeshGraph2D &graph2D);
            
            /**
             * \brief The default destructor.
             */  
            ~ExpList1D();
            
            /**
             * \brief 
             */  
            void   PhysDeriv  (ExpList &S0,
                               ExpList &S1 = NullExpList, 
                               ExpList &S2 = NullExpList)
            {
                ExpList::PhysDeriv(S0,S1,S2);
            }

            void SetBoundaryConditionExpansion(const SpatialDomains::MeshGraph1D &graph1D,
                                               SpatialDomains::BoundaryConditions &bcs, 
                                               const std::string variable,
                                               Array<OneD, LocalRegions::PointExpSharedPtr> &bndCondExpansions,
                                               Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions);
            
            void GetPeriodicVertices(const SpatialDomains::MeshGraph1D &graph1D,
                                     SpatialDomains::BoundaryConditions &bcs, 
                                     const std::string variable,
                                     map<int,int>& periodicVertices);
            
            void EvaluateBoundaryConditions(const NekDouble time,
                                            Array<OneD, LocalRegions::PointExpSharedPtr> &bndCondExpansions,
                                            Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions);

			
			/**
			 * \brief This function does the post-processing on a specified element
			 * \param kernel is the post-processing kernel
			 * \param inarray is the set of evaluation points
			 * \param outarray contains the resulting post-processed solution for the element \param elmId
			 * \param elmId specifies which element to perform the post-processing on
			 * \param h is the mesh spacing 
			 */
			void PostProcess(LibUtilities::KernelSharedPtr kernel, 
							 Array<OneD,NekDouble> &inarray,
							 Array<OneD,NekDouble> &outarray, 
							 int elmId,
							 NekDouble h);

			/**
			 * \brief This function does the post-processing on the entire domain
			 * \param kernel is the post-processing kernel
			 * \param inarray contains the set of evaluation points
			 * \param outarray contains the resulting post-procesed solution for the entire domain
			 * \param h is the mesh spacing 
			 */
			void PostProcess( LibUtilities::KernelSharedPtr kernel,
							  Array<OneD,NekDouble> &inarray,
							  Array<OneD,NekDouble> &outarray,
							  NekDouble h);

			
			/**
			 * \brief This function evaluates the global spectral/hp expansion
			 * at some arbitray set of points
			 
			 * Given the elemental coefficients \f$\hat{u}_n^e\f$ of an expansion, this 
             * function periodically evaluates the spectral/hp expansion 
             * \f$u^{\delta}(\boldsymbol{x})\f$ at arbitrary points

			 * * \param inarray1 An array of size \f$N_{\mathrm{eof}}\f$ containing the local 
             * coefficients \f$\hat{u}_n^e\f$.
			 * \param inarray2 contains the set of evaluation points
			 * \param h is the mesh spacing
			 * \param nmodes is the number of polynomial modes for each element
			 * (we consider that each element has the same number of polynomial modes)
			 * \param outarray contains the resulting values at the evaluation points
			  
			 */
			void PeriodicEval(Array<OneD,NekDouble> &inarray1, Array<OneD,NekDouble> &inarray2,
							  NekDouble h, int nmodes,
							  Array<OneD,NekDouble> &outarray);


			            
        protected:
            
        private:

        };

        typedef boost::shared_ptr<ExpList1D>      ExpList1DSharedPtr;
        typedef std::vector<ExpList1DSharedPtr>   ExpList1DVector;
        typedef std::vector<ExpList1DSharedPtr>::iterator ExpList1DVectorIter;

        const static Array<OneD, ExpList1DSharedPtr> NullExpList1DSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif//NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

/**
 * $Log: ExpList1D.h,v $
 * Revision 1.23  2008/09/16 13:36:06  pvos
 * Restructured the LocalToGlobalMap classes
 *
 * Revision 1.22  2008/08/14 22:15:51  sherwin
 * Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
 *
 * Revision 1.21  2008/07/29 22:27:33  sherwin
 * Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
 *
 * Revision 1.20  2008/06/23 14:21:01  pvos
 * updates for 1D ExpLists
 *
 * Revision 1.19  2008/05/10 18:27:33  sherwin
 * Modifications necessary for QuadExp Unified DG Solver
 *
 * Revision 1.18  2007/12/06 22:52:30  pvos
 * 2D Helmholtz solver updates
 *
 * Revision 1.17  2007/09/25 14:25:29  pvos
 * Update for helmholtz1D with different expansion orders
 *
 * Revision 1.16  2007/09/03 19:58:31  jfrazier
 * Formatting.
 *
 * Revision 1.15  2007/07/22 23:04:20  bnelson
 * Backed out Nektar::ptr.
 *
 * Revision 1.14  2007/07/20 02:04:12  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.13  2007/07/10 08:54:30  pvos
 * Updated ContField1D constructor
 *
 * Revision 1.12  2007/07/06 18:39:34  pvos
 * ContField1D constructor updates
 *
 * Revision 1.11  2007/06/05 16:36:55  pvos
 * Updated Explist2D ContExpList2D and corresponding demo-codes
 *
 **/
