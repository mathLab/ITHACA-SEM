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

#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/ExpList1D.h>

#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSys.h>


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
        class ContField2D: public DisContField2D
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
             * \brief 
             */ 
            ContField2D(const ContField2D &In, 
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const int bc_loc = 0);

            /**
             * \brief 
             */ 
            ContField2D::ContField2D(SpatialDomains::MeshGraph2D &graph2D);
          
            /**
             * \brief This constructor sets up global continuous field
             * based on an input mesh and boundary conditions.
             *
             * Given a mesh \a graph2D, containing information about
             * the domain and the spectral/hp element expansion, this
             * constructor fills the list of local expansions #m_exp
             * with the proper expansions, calculates the total number
             * of quadrature points \f$\boldsymbol{x}_i\f$ and local
             * expansion coefficients \f$\hat{u}^e_n\f$ and allocates
             * memory for the arrays #m_coeffs and
             * #m_phys. Furthermore, it constructs the mapping array
             * (contained in #m_locToGloMap) for the transformation
             * between local elemental level and global level, it
             * calculates the total number global expansion
             * coefficients \f$\hat{u}_n\f$ and allocates memory for
             * the array #m_contCoeffs. The constructor also
             * discretises the boundary conditions, specified by the
             * argument \a bcs, by expressing them in terms of the
             * coefficient of the expansion on the boundary.
             *
             * \param graph2D A mesh, containing information about the
             * domain and the spectral/hp element expansion.  \param
             * bcs The boundary conditions.  \param variable An
             * optional parameter to indicate for which variable the
             * field should be constructed.
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
          

            bool SameTypeOfBoundaryConditions(const ContField2D &In);

 
            /**
             * \brief This function returns (a reference to) the array 
             * \f$\boldsymbol{\hat{u}}_g\f$ (implemented as #m_contCoeffs) containing all 
             * global expansion coefficients.
             *
             * If one wants to get hold of the underlying data without modifying them, 
             * rather use the function #GetContCoeffs instead.
             *
             * \return (A reference to) the array #m_contCoeffs.
             */  
            inline Array<OneD, NekDouble> &UpdateContCoeffs()
            {
                m_transState = eContinuous;
                return m_contCoeffs;
            }

            /**
             * \brief This function returns (a reference to) the array 
             * \f$\boldsymbol{\hat{u}}_g\f$ (implemented as #m_contCoeffs) containing all 
             * global expansion coefficients.
             *
             * As the function returns a constant reference to a <em>const Array</em>, it is not 
             * possible to modify the underlying data of the array #m_contCoeffs. In order to 
             * do so, use the function #UpdateContCoeffs instead.
             *
             * \return (A reference to) the array #m_contCoeffs.
             */  
            inline const Array<OneD, const NekDouble> &GetContCoeffs() const
            {
                return m_contCoeffs;
            }
         
            /**
             * \brief This function returns the total number of global degrees of freedom 
             * \f$N_{\mathrm{dof}}\f$.
             *
             * \return  #m_contNcoeffs, the total number of global degrees of 
             * freedom.
             */  
            inline int GetContNcoeffs()
            {
                return m_contNcoeffs;
            }
             
            /**
             * \brief This function scatters from the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             * 
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}^{e}[i] = \mbox{sign}[e][i] \cdot 
             * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]$ \\
             *  \>     \>      continue \\
             *  \> continue
             *\f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the scatter operation 
             * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where 
             * \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * Note that the array #m_contCoeffs should be filled with the global 
             * coefficients \f$\boldsymbol{\hat{u}}_g\f$ and that the resulting local 
             * coefficients \f$\boldsymbol{\hat{u}}_l\f$ will be stored in #m_coeffs.
             */  
            inline void GlobalToLocal()
            {
                m_locToGloMap->GlobalToLocal(m_contCoeffs,m_coeffs);
            }
 
            /**
             * \brief This function scatters from the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             * 
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}^{e}[i] = \mbox{sign}[e][i] \cdot 
             * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]$ \\
             *  \>     \>      continue \\
             *  \> continue
             *\f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the scatter operation 
             * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where 
             * \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * \param outarray The resulting local degrees of freedom \f$\boldsymbol{x}_l\f$
             * will be stored in this array of size \f$N_\mathrm{eof}\f$.
             */  
            inline const void GlobalToLocal(Array<OneD,NekDouble> &outarray) const
            {
                m_locToGloMap->GlobalToLocal(m_contCoeffs,outarray);
            }


            /**
             * \brief This function scatters from the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             * 
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}^{e}[i] = \mbox{sign}[e][i] \cdot 
             * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]$ \\
             *  \>     \>      continue \\
             *  \> continue
             *\f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the scatter operation 
             * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where 
             * \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * \param inarray An array of size \f$N_\mathrm{dof}\f$ containing the global 
             * degrees of freedom \f$\boldsymbol{x}_g\f$.
             * \param outarray The resulting local degrees of freedom \f$\boldsymbol{x}_l\f$
             * will be stored in this array of size \f$N_\mathrm{eof}\f$.
             */  
            inline const void GlobalToLocal(const Array<OneD, const NekDouble> &inarray,
                                                  Array<OneD,       NekDouble> &outarray) const 
            {
                m_locToGloMap->GlobalToLocal(inarray,outarray);
            }

            /**
             * \brief This function gathers the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ from the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             *
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] = 
             *  \mbox{sign}[e][i] \cdot \boldsymbol{\hat{u}}^{e}[i]$\\
             *  \>     \>      continue\\
             *  \> continue
             * \f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the gather operation 
             * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{-1}\boldsymbol{\hat{u}}_l\f$, 
             * where \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * Note that the array #m_coeffs should be filled with the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$ and that the resulting global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ will be stored in #m_contCoeffs.
            */
            inline void LocalToGlobal()
            {
                m_locToGloMap->LocalToGlobal(m_coeffs,m_contCoeffs);
            }        
         
            /**
             * \brief This function assembles the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ from the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             *
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] = 
             * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot 
             * \boldsymbol{\hat{u}}^{e}[i]$\\
             *  \>     \>      continue\\
             *  \> continue
             * \f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the gather operation 
             * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$, where 
             * \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * Note that the array #m_coeffs should be filled with the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$ and that the resulting global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ will be stored in #m_contCoeffs.
             */  
            inline void Assemble()
            {
                m_locToGloMap->Assemble(m_coeffs,m_contCoeffs);
            }
 
            /**
             * \brief This function assembles the global coefficients 
             * \f$\boldsymbol{\hat{u}}_g\f$ from the local coefficients 
             * \f$\boldsymbol{\hat{u}}_l\f$.
             *
             * This operation is evaluated as:
             * \f{tabbing}
             * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
             *  \>     \> Do \= $i=$  $0,N_m^e-1$ \\
             *  \>     \>      \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] = 
             * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot 
             * \boldsymbol{\hat{u}}^{e}[i]$\\ 
            *  \>     \>      continue\\
             *  \> continue
             * \f}
             * where \a map\f$[e][i]\f$ is the mapping array and \a sign\f$[e][i]\f$ is an 
             * array of similar dimensions ensuring the correct modal connectivity between 
             * the different elements (both these arrays are contained in the data member 
             * #m_locToGloMap). This operation is equivalent to the gather operation 
             * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$, where 
             * \f$\mathcal{A}\f$ is the \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ 
             * permutation matrix.
             *
             * \param inarray An array of size \f$N_\mathrm{eof}\f$ containing the local 
             * degrees of freedom \f$\boldsymbol{x}_l\f$.
             * \param outarray The resulting global degrees of freedom 
             * \f$\boldsymbol{x}_g\f$ will be stored in this array of size 
             * \f$N_\mathrm{dof}\f$.
             */  
            inline const void Assemble(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray) const
            {
                m_locToGloMap->Assemble(inarray,outarray);
            }
 
            /**
             * \brief This function returns the map from local to global level.
             */  
            inline const LocalToGlobalC0ContMapSharedPtr& GetLocalToGlobalMap() const
            {
                return  m_locToGloMap;
            }


        
            /**
             * \brief This function calculates the inner product of a function 
             * \f$f(\boldsymbol{x})\f$ with respect to all <em>global</em> expansion modes 
             * \f$\phi_n^e(\boldsymbol{x})\f$.
             * 
             * The operation is evaluated locally (i.e. with respect to all local expansion 
             * modes) by the function ExpList#IProductWRTBase. The inner product with 
             * respect to the global expansion modes is than obtained by a global assembly 
             * operation.
             *
             * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the 
             * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the variable 
             * #m_phys of the ExpList object \a in. The result is stored in the array 
             * #m_contCoeffs.
             *
             * \param In An ExpList, containing the discrete evaluation of 
             * \f$f(\boldsymbol{x})\f$ at the quadrature points in its array #m_phys.
             */  
            void IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                       Array<OneD, NekDouble> &outarray,
                                 bool  UseContCoeffs = false)
            {
                if(UseContCoeffs)
                {
                    bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(StdRegions::eIProductWRTBase);

                    if(doGlobalOp)
                    {
                        GlobalMatrixKey gkey(StdRegions::eIProductWRTBase,m_locToGloMap);
                        GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                        mat->Multiply(inarray,outarray);
                    }
                    else
                    {
                        Array<OneD, NekDouble> wsp(m_ncoeffs);
                        IProductWRTBase_IterPerExp(inarray,wsp);
                        Assemble(wsp,outarray);
                    }
                }
                else
                {
                    IProductWRTBase_IterPerExp(inarray,outarray);
                }
            }

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
            void FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray,
                          bool  UseContCoeffs = false);
            

            /**
             * \brief This function performs the backward transformation of the spectral/hp 
             * element expansion.
             *
             * Given the coefficients of an expansion, this function evaluates the 
             * spectral/hp expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the quadrature 
             * points \f$\boldsymbol{x}_i\f$. This operation is evaluated locally by the 
             * function ExpList#BwdTrans.
             *
             * The coefficients of the expansion should be contained in the variable 
             * #m_coeffs of the ExpList object \a In. The resulting physical values at the 
             * quadrature points \f$u^{\delta}(\boldsymbol{x}_i)\f$ are stored in the array 
             * #m_phys.
             *
             * \param In An ExpList, containing the local coefficients \f$\hat{u}_n^e\f$ 
             * in its array #m_coeffs.
             */  
            void BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                Array<OneD,       NekDouble> &outarray,
                          bool  UseContCoeffs = false)
            {
                if(UseContCoeffs)
                {
                    bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(StdRegions::eBwdTrans);

                    if(doGlobalOp)
                    {
                        GlobalMatrixKey gkey(StdRegions::eBwdTrans,m_locToGloMap);
                        GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                        mat->Multiply(inarray,outarray);
                    }
                    else
                    {
                        Array<OneD, NekDouble> wsp(m_ncoeffs);
                        GlobalToLocal(inarray,wsp);
                        BwdTrans_IterPerExp(wsp,outarray);
                    }
                }
                else
                {
                    BwdTrans_IterPerExp(inarray,outarray);
                }
            }


            void MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, 
                                         Array<OneD,  NekDouble> &outarray,
                                         bool  UseContCoeffs = false);
            /**
             * \brief This function solves the two-dimensional
             * Helmholtz equation, subject to the boundary conditions
             * specified.
             * 
             * Consider the two dimensional Helmholtz equation,
             * \f[\nabla^2u(\boldsymbol{x})-\lambda u(\boldsymbol{x})
             * = f(\boldsymbol{x}),\f] supplemented with appropriate
             * boundary conditions (which are contained in the data
             * member #m_bndCondExpansions). Applying a \f$C^0\f$
             * continuous Galerkin discretisation, this equation leads
             * to the following linear system: \f[\left(
             * \boldsymbol{L}+\lambda\boldsymbol{M}\right)
             * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f] where
             * \f$\boldsymbol{L}\f$ and \f$\boldsymbol{M}\f$ are the
             * Laplacian and mass matrix respectively. This function
             * solves the system above for the global coefficients
             * \f$\boldsymbol{\hat{u}}\f$ by a call to the function
             * #GlobalSolve.
             *
             * The values of the function \f$f(\boldsymbol{x})\f$
             * evaluated at the quadrature points
             * \f$\boldsymbol{x}_i\f$ should be contained in the
             * variable #m_phys of the ExpList object \a Sin. The
             * resulting global coefficients
             * \f$\boldsymbol{\hat{u}}_g\f$ are stored in the array
             * #m_contCoeffs.
             * 
             * \param Sin An ExpList, containing the discrete
             * evaluation of the forcing function
             * \f$f(\boldsymbol{x})\f$ at the quadrature points in its
             * array #m_phys.  \param lambda The parameter
             * \f$\lambda\f$ of the Helmholtz equation
             */ 
            void HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                           NekDouble lambda,
                           bool      UseContCoeffs = false,
                           const Array<OneD, const NekDouble>& dirForcing = NullNekDouble1DArray);

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
            void LaplaceSolve(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,       NekDouble> &outarray,
                              const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray,
                              const Array<OneD,       Array<OneD,NekDouble> >& variablecoeffs = NullNekDoubleArrayofArray,
                              NekDouble time = 0.0,
                              bool UseContCoeffs = false);


            void LinearAdvectionSolve(const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                      NekDouble ax,     NekDouble ay,
                                      bool        UseContCoeffs = false,
                                      const Array<OneD, const NekDouble>& dirForcing = NullNekDouble1DArray);
            
            void LinearAdvectionEigs(const NekDouble ax, 
                                     const NekDouble ay,
                                     Array<OneD, NekDouble> &Real, 
                                     Array<OneD, NekDouble> &Imag, 
                                     Array<OneD, NekDouble> &Evecs = NullNekDouble1DArray);
          
            /**
             * \brief This function evaluates the boundary conditions
             * at a certain time-level.
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
                

            /**
             * \brief This function calculates the result of the multiplication of a global 
             * matrix of type specified by \a mkey with a vector given by \a inarray.
             *
             * This is equivalent to the operation:
             * \f[\boldsymbol{M\hat{u}}_g\f]
             * where \f$\boldsymbol{M}\f$ is the global matrix of type specified by \a mkey.
             * After scattering the global array \a inarray to local level, this operation 
             * is evaluated locally by the function ExpList#GeneralMatrixOp. The global 
             * result is then obtained by a global assembly procedure.
             *
             * \param mkey This key uniquely defines the type matrix required for the
             * operation.
             * \param inarray The vector \f$\boldsymbol{\hat{u}}_g\f$ of size 
             * \f$N_{\mathrm{dof}}\f$.
             * \param outarray The resulting vector of size \f$N_{\mathrm{dof}}\f$.
             */  
            void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                 const Array<OneD,const NekDouble> &inarray, 
                                       Array<OneD,      NekDouble> &outarray,
                                 bool  UseContCoeffs = false)
            {
                if(UseContCoeffs)
                {
                    bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(gkey.GetMatrixType());

                    if(doGlobalOp)
                    {
                        GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                        mat->Multiply(inarray,outarray);
                    }
                    else
                    {
                        Array<OneD,NekDouble> tmp1(2*m_ncoeffs);
                        Array<OneD,NekDouble> tmp2(tmp1+m_ncoeffs);
                        GlobalToLocal(inarray,tmp1);
                        GeneralMatrixOp_IterPerExp(gkey,tmp1,tmp2);
                        Assemble(tmp2,outarray);
                    }
                }
                else
                {
                    GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
                }
            }

            inline int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
            {
                ASSERTL1(gkey.LocToGloMapIsDefined(),
                         "To use method must have a LocalToGlobalBaseMap "
                         "attached to key");
                
                GlobalMatrixMap::iterator matrixIter = m_globalMat->find(gkey);
                
                if(matrixIter == m_globalMat->end())
                {
                    return 0;
                }
                else
                {
                    return matrixIter->second->GetMatrix()->GetNumNonZeroEntries();
                }
                
                return 0;
            }

        protected:

        private:  
            /**
             * \brief (A shared pointer to) the object which contains
             * all the required information for the transformation
             * from local to global degrees of freedom.
             */  
            LocalToGlobalC0ContMapSharedPtr m_locToGloMap;
  
            /**
             * \brief The total number of global degrees of freedom. 
             * #m_contNcoeffs\f$=N_{\mathrm{dof}}\f$
             */  
      	    int                       m_contNcoeffs;

  
            /**
             * \brief The array of length
             * #m_ncoeffs\f$=N_{\mathrm{dof}}\f$ containing the global
             * expansion coefficients.
             */  
	    Array<OneD, NekDouble>    m_contCoeffs;

            /**
             * \brief (A shared pointer to) a list which collects all
             * the global matrices being assembled, such that they
             * should be constructed only once.
             */  
            GlobalMatrixMapShPtr      m_globalMat;
 
            /**
             * \brief (A shared pointer to) a list which collects all
             * the global linear system being assembled, such that they
             * should be constructed only once.
             */  
            GlobalLinSysMapShPtr      m_globalLinSys;


            GlobalMatrixSharedPtr GetGlobalMatrix(const GlobalMatrixKey &mkey);
          
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
             * \brief This function solves the linear system specified
             * by the key \a key. 
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
                                   Array<OneD,        NekDouble> &inout,
                             const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);
            
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
                                     bool      UseContCoeffs,
                                     const Array<OneD, const NekDouble>& dirForcing)
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

        };
        typedef boost::shared_ptr<ContField2D>      ContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD2D_H
