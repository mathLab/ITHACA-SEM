///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.h
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

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/ExpList1D.h>

#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/BoundaryConditions.h>


namespace Nektar
{
    namespace MultiRegions
    {           
        /**
         * \brief This class is the abstraction of a global continuous one-dimensional 
         * spectral/hp element expansion which approximates the solution of a set of  
         * partial differential equations.
         *
         * As opposed to the class #ContExpList1D, the class #ContField1D is able to 
         * incorporate the boundary conditions imposed to the problem to be solved. 
         * Therefore, the class is equipped with three additional data members:
         * - #m_bndCondExpansions
         * - #m_bndTypes
         * - #m_bndCondEquations
         *
         * The first data structure, #m_bndCondExpansions, 
         * contains the point Expansion on the boundary,  #m_bndTypes
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
        class ContField1D: public ContExpList1D
        {
        public:           
            /**
             * \brief The default constructor. 
             */ 
            ContField1D();
          
            /**
             * \brief 
             */ 
            ContField1D(SpatialDomains::MeshGraph1D &graph1D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const int bc_loc = 0);
  
            /**
             * \brief This constructor sets up global continuous field based on an 
             * input mesh and boundary conditions.
             *
             * Given a mesh \a graph1D, containing information about the domain and 
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
             * \param graph1D A mesh, containing information about the domain and 
             * the spectral/hp element expansion.
             * \param bcs The boundary conditions.
             * \param variable An optional parameter to indicate for which variable 
             * the field should be constructed.
             */ 
            ContField1D(SpatialDomains::MeshGraph1D &graph1D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const std::string variable);
      
            /**
             * \brief 
             */ 
            ContField1D(const LibUtilities::BasisKey &Ba,
                const SpatialDomains::MeshGraph1D &graph1D,
                SpatialDomains::BoundaryConditions &bcs,
                const int bc_loc = 0);
      
            /**
             * \brief 
             */ 
            ContField1D(const LibUtilities::BasisKey &Ba,
                const SpatialDomains::MeshGraph1D &graph1D,
                SpatialDomains::BoundaryConditions &bcs,
                const std::string variable);

            /**
             * \brief The copy constructor.
             */ 
            ContField1D(const ContField1D &In);

            /**
             * \brief The default destructor.
             */ 
            ~ContField1D();

            /**
             * \brief This function performs the global forward transformation of a 
             * function \f$f(x)\f$, subject to the boundary conditions 
             * specified.
             *
             * Given a function \f$f(x)\f$ defined at the quadrature 
             * points, this function determines the unknown global coefficients 
             * \f$\boldsymbol{\hat{u}}^{\mathcal{H}}\f$ employing a discrete 
             * Galerkin projection from physical space to coefficient 
             * space. The operation is evaluated by the function #GlobalSolve using 
             * the global mass matrix.
             *
             * The values of the function \f$f(x)\f$ evaluated at the 
             * quadrature points \f$x_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a Sin. The resulting global 
             * coefficients \f$\hat{u}_g\f$ are stored in the array #m_contCoeffs.
             *
             * \param Sin An ExpList, containing the discrete evaluation of 
             * \f$f(x)\f$ at the quadrature points in its array #m_phys.
             */ 
            void FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray,
                          bool  UseContCoeffs = false);
            
            void MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD,       NekDouble> &outarray,
                                         bool  UseContCoeffs = false);
          
            /**
             * \brief This function solves the one-dimensional Helmholtz equation, 
             * subject to the boundary conditions specified.
             * 
             * Consider the one dimensional Helmholtz equation, 
             * \f[\frac{d^2u}{dx^2}-\lambda u(x) = 
             * f(x),\f]
             * supplemented with appropriate boundary conditions (which are contained
             * in the data member #m_bndCondExpansions). Applying a \f$C^0\f$ continuous 
             * Galerkin discretisation, this equation leads to the following linear 
             * system:
             * \f[\left( \boldsymbol{M}+\lambda\boldsymbol{L}\right)
             * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
             * where \f$\boldsymbol{M}\f$ and \f$\boldsymbol{L}\f$ are the mass and 
             * Laplacian matrix respectively. This function solves the system above 
             * for the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to 
             * the function #GlobalSolve.
             *
             * The values of the function \f$f(x)\f$ evaluated at the 
             * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a Sin. The resulting global 
             * coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the array 
             * #m_contCoeffs.
             * 
             * \param Sin An ExpList, containing the discrete evaluation of the 
             * forcing function \f$f(x)\f$ at the quadrature points 
             * in its array #m_phys.
             * \param lambda The parameter \f$\lambda\f$ of the Helmholtz equation
             */ 
            void HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                           NekDouble lambda,
                           bool UseContCoeffs = false,
                           Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray);

            /**
             * \brief This function evaluates the boundary conditions at a certain 
             * time-level.
             *
             * Based on the expression \f$g(x,t)\f$ for the boundary conditions, this
             * function evaluates the boundary conditions for all boundaries at 
             * time-level \a t.
             *
             * \param time The time at which the boundary conditions should be 
             * evaluated
             */ 
            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ExpList1D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
            };
          
            /**
             * \brief This function return the boundary conditions expansion.
             */ 
            inline const Array<OneD,const LocalRegions::PointExpSharedPtr>& GetBndCondExpansions()
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
             * \brief An object which contains the discretised boundary conditions. 
             *
             * It is an array of size equal to the number of boundary regions and 
             * consists of entries of the type LocalRegions#PointExp. 
             */ 
            Array<OneD,LocalRegions::PointExpSharedPtr>         m_bndCondExpansions;
          
            /**
             * \brief An array which contains the information about the boundary condition  
             * on the different boundary regions.
             */ 
            Array<OneD,SpatialDomains::BoundaryConditionShPtr>  m_bndConditions;
          
            /**
             * \brief This function returns the linear system specified by the key 
             * \a mkey.
             * 
             * The function searches the map #m_globalLinSys to see if the global matrix 
             * has been created before. If not, it calls the function  
             #GenGlobalLinSys to generate the requested global system.
             *
             * \param mkey This key uniquely defines the requested linear system.
             */ 
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

          
            /**
             * \brief This function solves the linear system specified by the key 
             * \a key.
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
                             Array<OneD, NekDouble> &dirForcing = NullNekDouble1DArray);

            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray,
                                    bool  UseContCoeffs)
            {
                FwdTrans(inarray,outarray,UseContCoeffs);
            }

            virtual void v_MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD,       NekDouble> &outarray,
                                                   bool  UseContCoeffs)
            {
                MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);
            }

            virtual void v_HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,       NekDouble> &outarray,
                                     NekDouble lambda,
                                     bool UseContCoeffs,
                                     Array<OneD, NekDouble>& dirForcing)
            {
                HelmSolve(inarray,outarray,lambda,UseContCoeffs,dirForcing);
            }
          
            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return GetBndConditions();
            }


            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                EvaluateBoundaryConditions(time);
            }


            /**
             * \brief This function discretises the boundary conditions by setting up
             * a list of point expansions.    
             *
             * The point expansions of the Dirichlet boundary regions are listed 
             * first in the array #m_bndCondExpansions.
             *
             * \param graph1D A mesh, containing information about the domain and 
             * the spectral/hp element expansion.
             * \param bcs An entity containing information about the boundaries and 
             * boundary conditions.
             * \param variable An optional parameter to indicate for which variable 
             * the boundary conditions should be discretised.
             */ 
            void GenerateBoundaryConditionExpansion(const SpatialDomains::MeshGraph1D        &graph1D,
                                                          SpatialDomains::BoundaryConditions &bcs, 
                                                    const std::string variable);
            
        };
        typedef boost::shared_ptr<ContField1D>      ContField1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTSOLNFIELD1D_H
