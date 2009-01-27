///////////////////////////////////////////////////////////////////////////////
//
// File ContField2D.h
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
// Description: Field definition in tow-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>

#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>


namespace Nektar
{
    namespace MultiRegions
    {           
        /**
         * \brief This class is the abstraction of a global continuous two-dimensional 
         * spectral/hp element expansion which approximates the solution of a set of  
         * partial differential equations.
         *
         * As opposed to the class #ContExpList2D, the class #ContField2D is able to 
         * incorporate the boundary conditions imposed to the problem to be solved. 
         * Therefore, the class is equipped with three additional data members:
         * - #m_bndCondExpansions
         * - #m_bndTypes
         * - #m_bndCondEquations
         *
         * The first data structure, #m_bndCondExpansions, 
         * contains the one-dimensional spectral/hp expansion on the boundary,  #m_bndTypes
         * stores information about the type of boundary condition on the different parts
         * of the boundary while #m_bndCondEquations holds the equation of the imposed
         * boundary conditions.<BR>
         * Furthermore, in case of Dirichlet boundary conditions,
         * this class is capable of lifting a known solution satisfying these boundary 
         * conditions. If we denote the unknown solution by 
         * \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet boundary conditions
         *  by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$, the expansion then can be decomposed 
         * as
         * \f[ u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
         * u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
         * \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
         * \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}-1}\hat{u}_n^{\mathcal{H}}
         * \Phi_n(\boldsymbol{x}_i).\f]
         * This lifting is accomplished by ordering the known global degrees of freedom, 
         * prescribed by the Dirichlet boundary conditions, first in the global array 
         * \f$\boldsymbol{\hat{u}}\f$, that is,
         * \f[\boldsymbol{\hat{u}}=\left[ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{u}}^{\mathcal{H}}
         * \end{array} \right].\f]
         * Such kind of expansions are also referred to as continuoous fields.
         * This class should be used when solving 2D problems using a standard Galerkin 
         * approach.
         */ 
        class ContField2D: public ContExpList2D
        {
        public:           
            /**
             * \brief The default constructor. 
             */ 
            ContField2D();
          
            /**
             * \brief 
             */ 
            ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                const int bc_loc = 0);
          
            /**
             * \brief This constructor sets up global continuous field based on an 
             * input mesh and boundary conditions.
             *
             * Given a mesh \a graph2D, containing information about the domain and 
             * the spectral/hp element expansion, this constructor fills the list of 
             * local expansions #m_exp with the proper expansions, calculates the 
             * total number of quadrature points \f$\boldsymbol{x}_i\f$ and local 
             * expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory for the 
             * arrays #m_coeffs and #m_phys. Furthermore, it constructs the mapping 
             * array (contained in #m_locToGloMap) for the transformation between 
             * local elemental level and global level, it calculates the total 
             * number global expansion coefficients \f$\hat{u}_n\f$ and allocates
             *  memory for the array #m_contCoeffs. The constructor also discretises 
             * the boundary conditions, specified by the argument \a bcs, by 
             * expressing them in terms of the coefficient of the expansion on the 
             * boundary. 
             *
             * \param graph2D A mesh, containing information about the domain and 
             * the spectral/hp element expansion.
             * \param bcs The boundary conditions.
             * \param variable An optional parameter to indicate for which variable 
             * the field should be constructed.
             */ 
            ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                const std::string variable);
          
            /**
             * \brief 
             */ 
            ContField2D(const LibUtilities::BasisKey &TriBa, 
                        const LibUtilities::BasisKey &TriBb, 
                        const LibUtilities::BasisKey &QuadBa, 
                        const LibUtilities::BasisKey &QuadBb, 
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const int bc_loc = 0,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);
                      
            /**
             * \brief 
             */ 
            ContField2D(const LibUtilities::BasisKey &TriBa, 
                        const LibUtilities::BasisKey &TriBb, 
                        const LibUtilities::BasisKey &QuadBa, 
                        const LibUtilities::BasisKey &QuadBb, 
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);
          
            /**
             * \brief The copy constructor.
             */ 
            ContField2D(const ContField2D &In);
          
            /**
             * \brief The default destructor.
             */ 
            ~ContField2D();
          
            /**
             * \brief This function performs the global forward transformation of a 
             * function \f$f(\boldsymbol{x})\f$, subject to the boundary conditions 
             * specified.
             *
             * Given a function \f$f(\boldsymbol{x})\f$ defined at the quadrature 
             * points, this function determines the unknown global coefficients 
             * \f$\boldsymbol{\hat{u}}^{\mathcal{H}}\f$ employing a discrete 
             * Galerkin projection from physical space to coefficient 
             * space. The operation is evaluated by the function #GlobalSolve using 
             * the global mass matrix.
             *
             * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the 
             * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a Sin. The resulting global 
             * coefficients \f$\hat{u}_g\f$ are stored in the array #m_contCoeffs.
             *
             * \param Sin An ExpList, containing the discrete evaluation of 
             * \f$f(\boldsymbol{x})\f$ at the quadrature points in its array #m_phys.
             */ 
            void FwdTrans (const ExpList &In);
            
            void MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool ContinuousArrays = true, bool ZeroBCs = false);
                      
            /**
             * \brief This function solves the two-dimensional Helmholtz equation, 
             * subject to the boundary conditions specified.
             * 
             * Consider the two dimensional Helmholtz equation, 
             * \f[\nabla^2u(\boldsymbol{x})-\lambda u(\boldsymbol{x}) = 
             * f(\boldsymbol{x}),\f]
             * supplemented with appropriate boundary conditions (which are contained
             * in the data member #m_bndCondExpansions). Applying a \f$C^0\f$ continuous 
             * Galerkin discretisation, this equation leads to the following linear 
             * system:
             * \f[\left( \boldsymbol{L}+\lambda\boldsymbol{M}\right)
             * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
             * where \f$\boldsymbol{L}\f$ and \f$\boldsymbol{M}\f$ are the Laplacian and 
             * mass matrix respectively. This function solves the system above 
             * for the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to 
             * the function #GlobalSolve.
             *
             * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the 
             * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a Sin. The resulting global 
             * coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the array 
             * #m_contCoeffs.
             * 
             * \param Sin An ExpList, containing the discrete evaluation of the 
             * forcing function \f$f(\boldsymbol{x})\f$ at the quadrature points 
             * in its array #m_phys.
             * \param lambda The parameter \f$\lambda\f$ of the Helmholtz equation
             */ 
            void HelmSolve(const ExpList &In, 
                           NekDouble lambda,
                           Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray);

            /**
             * \brief This function solves the two-dimensional Laplace equation, 
             * subject to the boundary conditions specified.
             * 
             * Consider the two dimensional Laplace equation, 
             * \f[\nabla\cdot\left(\boldsymbol{\sigma}\nabla u(\boldsymbol{x})\right) = 
             * f(\boldsymbol{x}),\f]
             * supplemented with appropriate boundary conditions (which are contained
             * in the data member #m_bndCondExpansions). In the equation above
             *  \f$\boldsymbol{\sigma}\f$ is the (symmetric positive definite) diffusion tensor:
             * \f[ \sigma = \left[ \begin{array}{cc} 
             * \sigma_{00}(\boldsymbol{x},t) & \sigma_{01}(\boldsymbol{x},t) \\
             * \sigma_{01}(\boldsymbol{x},t) & \sigma_{11}(\boldsymbol{x},t) 
             * \end{array} \right]. \f]
             * Applying a \f$C^0\f$ continuous 
             * Galerkin discretisation, this equation leads to the following linear 
             * system:
             * \f[\boldsymbol{L}
             * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
             * where \f$\boldsymbol{L}\f$ is the
             * Laplacian matrix. This function solves the system above 
             * for the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to 
             * the function #GlobalSolve.
             *
             * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the 
             * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a Sin. The resulting global 
             * coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the array 
             * #m_contCoeffs.
             * 
             * \param Sin An ExpList, containing the discrete evaluation of the 
             * forcing function \f$f(\boldsymbol{x})\f$ at the quadrature points 
             * in its array #m_phys.
             * \param variablecoeffs The (optional) parameter containing the coefficients
             * evaluated at the quadrature points. It is an Array of (three) arrays which 
             * stores the laplacian coefficients in the following way
             * \f[
             * \mathrm{variablecoeffs} = \left[ \begin{array}{c}
             * \left[\sigma_{00}(\boldsymbol{x_i},t)\right]_i \\
             * \left[\sigma_{01}(\boldsymbol{x_i},t)\right]_i \\
             * \left[\sigma_{11}(\boldsymbol{x_i},t)\right]_i 
             * \end{array}\right]
             * \f]
             * If this argument is not passed to the function, the following equation will 
             * be solved:
             * \f[\nabla^2u(\boldsymbol{x}) = f(\boldsymbol{x}),\f]
             * \param time The time-level at which the coefficients are evaluated
             */ 
            void LaplaceSolve(const ExpList &In, 
                              const Array<OneD, Array<OneD,NekDouble> >& variablecoeffs = NullNekDoubleArrayofArray,
                              NekDouble time = 0.0,
                              Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray);
          
            /**
             * \brief This function evaluates the boundary conditions at a certain 
             * time-level.
             *
             * Based on the boundary condition \f$g(\boldsymbol{x},t)\f$ evaluated
             * at a given time-level \a t, this function transforms the boundary 
             * conditions onto the coefficients of the (one-dimensional) boundary 
             * expansion. Depending on the type of boundary conditions, these
             * expansion coefficients are calculated in different ways:
             * - <b>Dirichlet boundary conditions</b><BR>
             *   In order to ensure global \f$C^0\f$ continuity of the spectral/hp 
             *   approximation, the Dirichlet boundary conditions are projected onto 
             *   the boundary expansion by means of a modified \f$C^0\f$ continuous  
             *   Galerkin projection. This projection can be viewed as a collocation
             *   projection at the vertices, followed by an \f$L^2\f$ projection on 
             *   the interior modes of the edges. The resulting coefficients 
             *   \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ will be stored for the 
             *   boundary expansion.
             * - <b>Neumann boundary conditions</b>
             *   In the discrete Galerkin formulation of the problem to be solved, 
             *   the Neumann boundary conditions appear as the set of surface 
             *   integrals: \f[\boldsymbol{\hat{g}}=\int_{\Gamma}
             *   \phi^e_n(\boldsymbol{x})g(\boldsymbol{x})d(\boldsymbol{x})\quad
             *   \forall n \f]
             *   As a result, it are the coefficients \f$\boldsymbol{\hat{g}}\f$ 
             *   that will be stored in the boundary expansion
             *
             * \param time The time at which the boundary conditions should be 
             * evaluated
             */ 
            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ExpList2D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
            }
          
            /**
             * \brief This function return the boundary conditions expansion.
             */ 
            inline const Array<OneD,const MultiRegions::ExpList1DSharedPtr>& GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }
            
            inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& GetBndConditions()
            {
                return m_bndConditions;
            }

            void GenerateDirBndCondForcing(const GlobalLinSysKey &key, 
                                           Array<OneD, NekDouble> &inout, 
                                           Array<OneD, NekDouble> &outarray);
                
        protected:

        private:  
            /**
             * \brief The number of boundary segments on which
             * Dirichlet boundary conditions are imposed
             */ 
            int m_numDirBndCondExpansions;
        
            /**
             * \brief An object which contains the discretised boundary conditions. 
             *
             * It is an array of size equal to the number of boundary regions and 
             * consists of entries of the type MultiRegions#ExpList1D. Every entry corresponds to 
             * the one-dimensional spectral/hp expansion on a single boundary region. 
             * The values of the boundary conditions are stored as the coefficients 
             * of the one-dimensional expansion.
             */ 
            Array<OneD,MultiRegions::ExpList1DSharedPtr>       m_bndCondExpansions;
          
            /**
             * \brief An array which contains the information about
             * the boundary condition on the different boundary
             * regions.
             */ 
            Array<OneD,SpatialDomains::BoundaryConditionShPtr>  m_bndConditions;
          
            /**
             * \brief This function returns the linear system specified by the key 
             * \a mkey.
             * 
             * The function searches the map #m_globalMat to see if the global matrix 
             * has been created before. If not, it calls the function  
             #GenGlobalLinSys to generate the requested global system.
             *
             * \param mkey This key uniquely defines the requested linear system.
             */ 
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);
          
            /**
             * \brief This function solves the linear system specified by the key 
             * \a key. NEEDS UPDATING (SJS)
             * 
             * Given a linear system specified by the key \a key,
             * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}},\f]
             * this function solves this linear system taking into account the 
             * boundary conditions specified in the data member #m_bndCondExpansions. 
             * Therefore, it adds an array \f$\boldsymbol{\hat{g}}\f$ which 
             * represents the non-zero surface integral resulting from the weak 
             * boundary conditions (e.g. Neumann boundary conditions) to the right 
             * hand side, that is,
             * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}+
             * \boldsymbol{\hat{g}}.\f]
             * Furthermore, it lifts the known degrees of freedom which are 
             * prescribed by the Dirichlet boundary conditions. As these known 
             * coefficients \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ are numbered 
             * first in the global coefficient array \f$\boldsymbol{\hat{u}}_g\f$, 
             * the linear system can be decomposed as,
             * \f[\left[\begin{array}{cc}
             * \boldsymbol{M}^{\mathcal{DD}}&\boldsymbol{M}^{\mathcal{DH}}\\
             * \boldsymbol{M}^{\mathcal{HD}}&\boldsymbol{M}^{\mathcal{HH}}
             * \end{array}\right]
             * \left[\begin{array}{c}
             * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
             * \boldsymbol{\hat{u}}^{\mathcal{H}}
             * \end{array}\right]=
             * \left[\begin{array}{c}
             * \boldsymbol{\hat{f}}^{\mathcal{D}}\\
             * \boldsymbol{\hat{f}}^{\mathcal{H}}
             * \end{array}\right]+
             * \left[\begin{array}{c}
             * \boldsymbol{\hat{g}}^{\mathcal{D}}\\
             * \boldsymbol{\hat{g}}^{\mathcal{H}}
             * \end{array}\right]
             * \f]
             * which will then be solved for the unknown coefficients 
             * \f$\boldsymbol{\hat{u}}^{\mathcal{H}}\f$ as,
             * \f[
             * \boldsymbol{M}^{\mathcal{HH}}\boldsymbol{\hat{u}}^{\mathcal{H}}=
             * \boldsymbol{\hat{f}}^{\mathcal{H}}+\boldsymbol{\hat{g}}^{\mathcal{H}}-
             * \boldsymbol{M}^{\mathcal{HD}}\boldsymbol{\hat{u}}^{\mathcal{D}}\f]
             *
             * \param mkey This key uniquely defines the linear system to be solved.
             * \param Sin An ExpList, containing the discrete evaluation of the 
             * forcing function \f$f(\boldsymbol{x})\f$ at the quadrature points in 
             * its array #m_phys.
             * \param ScaleForcing An optional parameter with which the forcing 
             * vector \f$\boldsymbol{\hat{f}}\f$ should be multiplied.
             */ 
            void GlobalSolve(const GlobalLinSysKey &key, 
                             const Array<OneD, const  NekDouble> &rhs, 
                             Array<OneD, NekDouble> &inout,
                             Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray);


          
            /**
             * \brief This function discretises the boundary conditions by setting up
             * a list of one-dimensional boundary expansions.    
             *
             * According to their boundary region, the separate segmental boundary 
             * expansions are bundled together in an object of the class 
             * MultiRegions#ExpList1D. 
             * The list of expansions of the Dirichlet boundary regions are listed 
             * first in the array #m_bndCondExpansions.
             *
             * \param graph2D A mesh, containing information about the domain and 
             * the spectral/hp element expansion.
             * \param bcs An entity containing information about the boundaries and 
             * boundary conditions.
             * \param variable An optional parameter to indicate for which variable 
             * the boundary conditions should be discretised.
             */ 
            void GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                    SpatialDomains::BoundaryConditions &bcs, 
                                                    const std::string variable);
            
            virtual void v_MultiplyByInvMassMatrix(const Array<OneD,const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool GlobalArrays, bool ZeroBCs)
            {
                MultiplyByInvMassMatrix(inarray,outarray,GlobalArrays, ZeroBCs);
            }

            virtual void v_FwdTrans(const ExpList &Sin)
            {
                FwdTrans(Sin);
            }

            virtual void v_HelmSolve(const ExpList &In, 
                                     NekDouble lambda,
                                     Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray)
            {
                HelmSolve(In,lambda,dirForcing);
            }

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return GetBndConditions();
            }

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                EvaluateBoundaryConditions(time);
            }


        };
        typedef boost::shared_ptr<ContField2D>      ContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD2D_H
