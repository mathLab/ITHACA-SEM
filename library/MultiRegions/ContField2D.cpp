///////////////////////////////////////////////////////////////////////////////
//
// File ContField2D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Field definition for 2D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <tuple>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ContField2D
         * The class #ContField2D is
         * able to incorporate the boundary conditions imposed to the problem
         * to be solved. Therefore, the class is equipped with three additional
         * data members:
         * - #m_bndCondExpansions
         * - #m_bndTypes
         * - #m_bndCondEquations
         *
         * The first data structure, #m_bndCondExpansions, contains the
         * one-dimensional spectral/hp expansion on the boundary,  #m_bndTypes
         * stores information about the type of boundary condition on the
         * different parts of the boundary while #m_bndCondEquations holds the
         * equation of the imposed boundary conditions.
         *
         * Furthermore, in case of Dirichlet boundary conditions, this class is
         * capable of lifting a known solution satisfying these boundary
         * conditions. If we denote the unknown solution by
         * \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet
         * boundary conditions by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$, the
         * expansion then can be decomposed as
         * \f[ u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
         * u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
         * \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
         * \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}-1}
         *  \hat{u}_n^{\mathcal{H}} \Phi_n(\boldsymbol{x}_i).\f]
         * This lifting is accomplished by ordering the known global degrees of
         * freedom, prescribed by the Dirichlet boundary conditions, first in
         * the global array
         * \f$\boldsymbol{\hat{u}}\f$, that is,
         * \f[\boldsymbol{\hat{u}}=\left[ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{u}}^{\mathcal{H}}
         * \end{array} \right].\f]
         * Such kind of expansions are also referred to as continuous fields.
         * This class should be used when solving 2D problems using a standard
         * Galerkin approach.
         */

        /**
         *
         */
        ContField2D::ContField2D():
            DisContField2D(),
            m_locToGloMap(),
            m_globalMat(),
            m_globalLinSysManager(
                std::bind(&ContField2D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
        }


        /**
         * Given a mesh \a graph2D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory
         * for the arrays #m_coeffs and #m_phys. Furthermore, it constructs the
         * mapping array (contained in #m_locToGloMap) for the transformation
         * between local elemental level and global level, it calculates the
         * total number global expansion coefficients \f$\hat{u}_n\f$ and
         * allocates memory for the array #m_contCoeffs. The constructor also
         * discretises the boundary conditions, specified by the argument \a
         * bcs, by expressing them in terms of the coefficient of the expansion
         * on the boundary.
         *
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   variable    An optional parameter to indicate for which
         *                      variable the field should be constructed.
         */
        ContField2D::ContField2D(
                         const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &graph2D,
                         const std::string &variable,
                         const bool DeclareCoeffPhysArrays,
                         const bool CheckIfSingularSystem,
                         const Collections::ImplementationType ImpType):
            DisContField2D(pSession,graph2D,variable,false,DeclareCoeffPhysArrays,ImpType),
            m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSysManager(
                std::bind(&ContField2D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
            m_locToGloMap = MemoryManager<AssemblyMapCG>
                ::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                    m_bndCondExpansions,
                                    m_bndConditions,
                                    CheckIfSingularSystem,
                                    variable,
                                    m_periodicVerts,
                                    m_periodicEdges);

            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                m_locToGloMap->PrintStats(std::cout, variable);
            }
        }


        /**
         * Given a mesh \a graph2D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory
         * for the arrays #m_coeffs and #m_phys. Furthermore, it constructs the
         * mapping array (contained in #m_locToGloMap) for the transformation
         * between local elemental level and global level, it calculates the
         * total number global expansion coefficients \f$\hat{u}_n\f$ and
         * allocates memory for the array #m_coeffs. The constructor also
         * discretises the boundary conditions, specified by the argument \a
         * bcs, by expressing them in terms of the coefficient of the expansion
         * on the boundary.
         *
         * @param   In          Existing ContField2D object used to provide the
         *                      local to global mapping information and
         *                      global solution type.
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc
         */
        ContField2D::ContField2D(const ContField2D &In,
                                 const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                 const std::string &variable,
                                 bool DeclareCoeffPhysArrays,
                                 const bool CheckIfSingularSystem):
            DisContField2D(In,graph2D,variable,false,DeclareCoeffPhysArrays),
            m_globalMat   (MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSysManager(
                std::bind(&ContField2D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
            if(!SameTypeOfBoundaryConditions(In) || CheckIfSingularSystem)
            {
                m_locToGloMap = MemoryManager<AssemblyMapCG>
                    ::AllocateSharedPtr(m_session, m_ncoeffs,*this,
                                        m_bndCondExpansions,
                                        m_bndConditions,
                                        CheckIfSingularSystem,
                                        variable,
                                        m_periodicVerts,
                                        m_periodicEdges);

                if (m_session->DefinesCmdLineArgument("verbose"))
                {
                    m_locToGloMap->PrintStats(std::cout, variable);
                }
            }
            else
            {
                m_locToGloMap = In.m_locToGloMap;
            }
        }


        /**
         * Initialises the object as a copy of an existing ContField2D object.
         * @param   In                       Existing ContField2D object.
         * @param DeclareCoeffPhysArrays     bool to declare if \a m_phys
         * and \a m_coeffs should be declared. Default is true
         */
        ContField2D::ContField2D(const ContField2D &In, bool DeclareCoeffPhysArrays):
            DisContField2D(In,DeclareCoeffPhysArrays),
            m_locToGloMap(In.m_locToGloMap),
            m_globalMat(In.m_globalMat),
            m_globalLinSysManager(In.m_globalLinSysManager)
        {
        }


        /**
         *
         */
        ContField2D::~ContField2D()
        {
        }


        /**
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
         * coefficients \f$\hat{u}_g\f$ are stored in the array #m_coeffs.
         *
         * @param   Sin         An ExpList, containing the discrete evaluation
         *                      of \f$f(\boldsymbol{x})\f$ at the quadrature
         *                      points in its array #m_phys.
         */
        void ContField2D::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray)

        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            IProductWRTBase(inarray,wsp);

            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
            
            GlobalSolve(key,wsp,outarray);
        }

        /**
         *
         */
        void ContField2D::v_SmoothField(Array<OneD,NekDouble> &field)
        {
            int Ncoeffs = m_locToGloMap->GetNumLocalCoeffs();
            Array<OneD,NekDouble> tmp1(Ncoeffs);
            Array<OneD,NekDouble> tmp2(Ncoeffs);

            IProductWRTBase(field,tmp1);
            MultiplyByInvMassMatrix(tmp1,tmp2);
            BwdTrans(tmp2,field);
        }


        /**
         * Computes the matrix vector product
         * @f$ \mathbf{y} = \mathbf{M}^{-1}\mathbf{x} @f$.
         *
         * @param   inarray     Input vector @f$\mathbf{x}@f$.
         * @param   outarray    Output vector @f$\mathbf{y}@f$.
         */
        void ContField2D::MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)

        {

            GlobalLinSysKey key(StdRegions::eMass,m_locToGloMap);
            GlobalSolve(key,inarray,outarray);
        }


        /**
         * Consider the two dimensional Laplace equation,
         * \f[\nabla\cdot\left(\boldsymbol{\sigma}\nabla
         * u(\boldsymbol{x})\right) = f(\boldsymbol{x}),\f] supplemented with
         * appropriate boundary conditions (which are contained in the data
         * member #m_bndCondExpansions). In the equation above
         * \f$\boldsymbol{\sigma}\f$ is the (symmetric positive definite)
         * diffusion tensor:
         * \f[ \sigma = \left[ \begin{array}{cc}
         * \sigma_{00}(\boldsymbol{x},t) & \sigma_{01}(\boldsymbol{x},t) \\
         * \sigma_{01}(\boldsymbol{x},t) & \sigma_{11}(\boldsymbol{x},t)
         * \end{array} \right]. \f]
         * Applying a \f$C^0\f$ continuous Galerkin discretisation, this
         * equation leads to the following linear system:
         * \f[\boldsymbol{L}
         * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
         * where \f$\boldsymbol{L}\f$ is the Laplacian matrix. This function
         * solves the system above for the global coefficients
         * \f$\boldsymbol{\hat{u}}\f$ by a call to the function #GlobalSolve.
         *
         * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the
         * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a Sin. The resulting global
         * coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the array
         * #m_coeffs.
         *
         * @param   Sin         An ExpList, containing the discrete evaluation
         *                      of the forcing function \f$f(\boldsymbol{x})\f$
         *                      at the quadrature points in its array #m_phys.
         * @param   variablecoeffs The (optional) parameter containing the
         *                      coefficients evaluated at the quadrature
         *                      points. It is an Array of (three) arrays which
         *                      stores the laplacian coefficients in the
         *                      following way
         * \f[\mathrm{variablecoeffs} = \left[ \begin{array}{c}
         * \left[\sigma_{00}(\boldsymbol{x_i},t)\right]_i \\
         * \left[\sigma_{01}(\boldsymbol{x_i},t)\right]_i \\
         * \left[\sigma_{11}(\boldsymbol{x_i},t)\right]_i
         * \end{array}\right]
         * \f]
         * If this argument is not passed to the function, the following
         * equation will be solved:
         * \f[\nabla^2u(\boldsymbol{x}) = f(\boldsymbol{x}),\f]
         *
         * @param   time        The time-level at which the coefficients are
         *                      evaluated
         */
        void ContField2D::LaplaceSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const Array<OneD, const NekDouble> &dirForcing,
                const Array<OneD,       Array<OneD,NekDouble> >& variablecoeffs,
                NekDouble time)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            IProductWRTBase(inarray,wsp);

            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_ncoeffs, wsp, 1);

            // Forcing function with weak boundary conditions
            int i,j;
            int bndcnt = 0;
            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();

            // Add weak boundary conditions to forcing
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eNeumann ||
                   m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eRobin)
                {
                    const Array<OneD, NekDouble> bndcoeff =
                        (m_bndCondExpansions[i])->GetCoeffs(); 

                    if(m_locToGloMap->GetSignChange())
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt + j]] += sign[bndcnt + j] * bndcoeff[j]; 
                        }
                    }
                    else
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt+j]] += bndcoeff[bndcnt + j]; 
                        }
                    }                    
                }

                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }

            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffD00] = variablecoeffs[0];
            varcoeffs[StdRegions::eVarCoeffD11] = variablecoeffs[3];
            varcoeffs[StdRegions::eVarCoeffD22] = variablecoeffs[5];
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorTime] = time;

            // Solve the system
            GlobalLinSysKey key(StdRegions::eLaplacian,m_locToGloMap,factors,
                                varcoeffs);

            GlobalSolve(key,wsp,outarray,dirForcing);
        }


        /**
         * Constructs the GlobalLinearSysKey for the linear advection operator
         * with the supplied parameters, and computes the eigenvectors and
         * eigenvalues of the associated matrix.
         * @param   ax          Advection parameter, x.
         * @param   ay          Advection parameter, y.
         * @param   Real        Computed eigenvalues, real component.
         * @param   Imag        Computed eigenvalues, imag component.
         * @param   Evecs       Computed eigenvectors.
         */
        void ContField2D::LinearAdvectionEigs(const NekDouble ax,
                                              const NekDouble ay,
                                              Array<OneD, NekDouble> &Real,
                                              Array<OneD, NekDouble> &Imag,
                                              Array<OneD, NekDouble> &Evecs)
        {
            // Solve the system
            Array<OneD, Array<OneD, NekDouble> > vel(2);
            Array<OneD, NekDouble> vel_x(m_npoints,ax);
            Array<OneD, NekDouble> vel_y(m_npoints,ay);
            vel[0] = vel_x;
            vel[1] = vel_y;

            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffVelX] = Array<OneD, NekDouble>(m_npoints,ax);
            varcoeffs[StdRegions::eVarCoeffVelY] = Array<OneD, NekDouble>(m_npoints,ay);
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorTime] = 0.0;
            GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,m_locToGloMap,
                                factors,varcoeffs);

            DNekMatSharedPtr   Gmat = GenGlobalMatrixFull(key,m_locToGloMap);
            Gmat->EigenSolve(Real,Imag,Evecs);
        }




        /**
         * Given a linear system specified by the key \a key,
         * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}},\f]
         * this function solves this linear system taking into account the
         * boundary conditions specified in the data member
         * #m_bndCondExpansions. Therefore, it adds an array
         * \f$\boldsymbol{\hat{g}}\f$ which represents the non-zero surface
         * integral resulting from the weak boundary conditions (e.g. Neumann
         * boundary conditions) to the right hand side, that is,
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
         * \boldsymbol{\hat{f}}^{\mathcal{H}}+
         * \boldsymbol{\hat{g}}^{\mathcal{H}}-
         * \boldsymbol{M}^{\mathcal{HD}}\boldsymbol{\hat{u}}^{\mathcal{D}}\f]
         *
         * @param   mkey        This key uniquely defines the linear system to
         *                      be solved.
         * @param   locrhs      contains the forcing term in local coefficient space
         * @note    inout contains initial guess and final output in local coeffs.
         */
        void ContField2D::GlobalSolve(
                                const GlobalLinSysKey &key,
                                const Array<OneD, const NekDouble>& locrhs,
                                      Array<OneD,       NekDouble>& inout,
                                const Array<OneD, const NekDouble>& dirForcing)
        {
            int NumDirBcs   = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();

            // STEP 1: SET THE DIRICHLET DOFS TO THE RIGHT VALUE
            //         IN THE SOLUTION ARRAY
            v_ImposeDirichletConditions(inout);
            
            // STEP 2: CALCULATE THE HOMOGENEOUS COEFFICIENTS
            if(contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                LinSys->Solve(locrhs,inout,m_locToGloMap,dirForcing);

            }
        }


        /**
         * Returns the global matrix associated with the given GlobalMatrixKey.
         * If the global matrix has not yet been constructed on this field,
         * it is first constructed using GenGlobalMatrix().
         * @param   mkey        Global matrix key.
         * @returns Assocated global matrix.
         */
        GlobalMatrixSharedPtr ContField2D::GetGlobalMatrix(
                                const GlobalMatrixKey &mkey)
        {
            ASSERTL1(mkey.LocToGloMapIsDefined(),
                     "To use method must have a AssemblyMap "
                     "attached to key");

            GlobalMatrixSharedPtr glo_matrix;
            auto matrixIter = m_globalMat->find(mkey);

            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalMatrix(mkey,m_locToGloMap);
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }


        /**
         * The function searches the map #m_globalLinSys to see if the
         * global matrix has been created before. If not, it calls the function
         * #GenGlobalLinSys to generate the requested global system.
         *
         * @param   mkey        This key uniquely defines the requested
         *                      linear system.
         */
        GlobalLinSysSharedPtr ContField2D::GetGlobalLinSys(
                                const GlobalLinSysKey &mkey)
        {
            return m_globalLinSysManager[mkey];
        }

        GlobalLinSysSharedPtr ContField2D::GenGlobalLinSys(
                                const GlobalLinSysKey &mkey)
        {
            ASSERTL1(mkey.LocToGloMapIsDefined(),
                     "To use method must have a AssemblyMap "
                     "attached to key");
            return ExpList::GenGlobalLinSys(mkey, m_locToGloMap);
        }


        /**
         *
         */
        void ContField2D::v_BwdTrans(
                                     const Array<OneD, const NekDouble>
                                     &inarray,
                                     Array<OneD,       NekDouble> &outarray)
        {
            BwdTrans(inarray,outarray);
        }


        /**
         *
         */
        void ContField2D::v_FwdTrans(
                                     const Array<OneD, const NekDouble>
                                     &inarray,
                                     Array<OneD,       NekDouble> &outarray)
        {
            FwdTrans(inarray,outarray);
        }

        void ContField2D::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            int i,j;
            int bndcnt=0;

            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eDirichlet)
                {
                    const Array<OneD, NekDouble> bndcoeff =
                        (m_bndCondExpansions[i])->GetCoeffs(); 

                    if(m_locToGloMap->GetSignChange())
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            outarray[map[bndcnt + j]] = sign[bndcnt + j] * bndcoeff[j]; 
                        }
                    }
                    else
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            outarray[map[bndcnt+j]] = bndcoeff[j]; 
                        }
                    }                    
                }
                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }
            
            m_locToGloMap->UniversalAbsMaxBnd(outarray);
            
            set<ExtraDirDof> &copyLocalDirDofs = m_locToGloMap->GetCopyLocalDirDofs();
            for (auto &it : copyLocalDirDofs)
            {
                outarray[std::get<0>(it)] =
                    outarray[std::get<1>(it)]*std::get<2>(it);
            }

        }

        void ContField2D::v_FillBndCondFromField(void)
        {
            int bndcnt = 0;

            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> bndmap= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(int i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if (m_bndConditions[i]->GetBoundaryConditionType() ==
                       SpatialDomains::ePeriodic)
                {
                    continue;
                }

                Array<OneD, NekDouble>& coeffs = m_bndCondExpansions[i]->UpdateCoeffs();
                
                if(m_locToGloMap->GetSignChange())
                {
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        coeffs[j] = sign[bndcnt+j] * m_coeffs[bndmap[bndcnt+j]];
                    }
                }
                else
                {
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        coeffs[j] = m_coeffs[bndmap[bndcnt+j]];
                    }
                }

                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }
        }


        void ContField2D::v_FillBndCondFromField(const int nreg)
        {
            int bndcnt = 0;

            ASSERTL1(nreg < m_bndCondExpansions.size(),
                     "nreg is out or range since this many boundary "
                     "regions to not exist");
            
            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> bndmap= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();

            // Now fill in all other Dirichlet coefficients.
            Array<OneD, NekDouble>& coeffs = m_bndCondExpansions[nreg]->UpdateCoeffs();

            for(int j = 0; j < nreg; ++j)
            {
                if (m_bndConditions[j]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                bndcnt += m_bndCondExpansions[j]->GetNcoeffs();
            }
            
            if(m_locToGloMap->GetSignChange())
            {
                for(int j = 0; j < (m_bndCondExpansions[nreg])->GetNcoeffs(); ++j)
                {
                    coeffs[j] = sign[bndcnt + j] * m_coeffs[bndmap[bndcnt + j]];
                }
            }
            else
            {
                for(int j = 0; j < (m_bndCondExpansions[nreg])->GetNcoeffs(); ++j)
                {
                    coeffs[j] = m_coeffs[bndmap[bndcnt + j]];
                }
            }
        }


        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}^{e}[i] = \mbox{sign}[e][i] \cdot
         * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]$ \\
         * \> \> continue \\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the scatter operation
         * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$,
         * where \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @note The array #m_coeffs should be filled with the global
         * coefficients \f$\boldsymbol{\hat{u}}_g\f$ and that the resulting
         * local coefficients \f$\boldsymbol{\hat{u}}_l\f$ will be stored in
         * #m_coeffs.
         */
        void ContField2D::v_GlobalToLocal(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->GlobalToLocal(inarray, outarray);
        }


        void ContField2D::v_GlobalToLocal(void)
        {
            m_locToGloMap->GlobalToLocal(m_coeffs,m_coeffs);
        }


        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
         * \mbox{sign}[e][i] \cdot \boldsymbol{\hat{u}}^{e}[i]$\\
         * \> \> continue\\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the gather operation
         * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{-1}\boldsymbol{\hat{u}}_l\f$,
         * where \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @note The array #m_coeffs should be filled with the local
         *          coefficients \f$\boldsymbol{\hat{u}}_l\f$ and that
         *          the resulting global coefficients
         *          \f$\boldsymbol{\hat{u}}_g\f$ will be stored in
         *          #m_coeffs. Also if useComm is set to false then no
         *          communication call will be made to check if all
         *          values are consistent over processors
         */

        void ContField2D::v_LocalToGlobal(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,NekDouble> &outarray,
            bool useComm)
        {
            m_locToGloMap->LocalToGlobal(inarray, outarray, useComm);
        }


        void ContField2D::v_LocalToGlobal(bool useComm)

        {
            m_locToGloMap->LocalToGlobal(m_coeffs,m_coeffs, useComm);
        }

        /**
         *
         */
        void ContField2D::v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            MultiplyByInvMassMatrix(inarray,outarray);
        }


        /**
         * Consider the two dimensional Helmholtz equation,
         * \f[\nabla^2u(\boldsymbol{x})-\lambda u(\boldsymbol{x})
         * = f(\boldsymbol{x}),\f] supplemented with appropriate boundary
         * conditions (which are contained in the data member
         * #m_bndCondExpansions). Applying a \f$C^0\f$ continuous Galerkin
         * discretisation, this equation leads to the following linear system:
         * \f[\left(\boldsymbol{L}+\lambda\boldsymbol{M}\right)
         * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f] where
         * \f$\boldsymbol{L}\f$ and \f$\boldsymbol{M}\f$ are the Laplacian and
         * mass matrix respectively. This function solves the system above for
         * the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to the
         * function #GlobalSolve. It is assumed #m_coeff contains an
         * initial estimate for the solution.
         *
         * The values of the function \f$f(\boldsymbol{x})\f$
         * evaluated at the quadrature points \f$\boldsymbol{x}_i\f$
         * should be contained in the variable #m_phys of the ExpList
         * object \a inarray. 
         *
         * @param   inarray     An ExpList, containing the discrete evaluation
         *                      of the forcing function \f$f(\boldsymbol{x})\f$
         *                      at the quadrature points in its array #m_phys.
         * @param   factors    The parameter \f$\lambda\f$ of the Helmholtz
         *                      equation is specified through the factors map
         */
        void ContField2D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing)

        {
            int i,j;

            //----------------------------------
            //  Setup RHS Inner product
            //----------------------------------
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            if(PhysSpaceForcing)
            {
                IProductWRTBase(inarray,wsp);
                // Note -1.0 term necessary to invert forcing function to
                // be consistent with matrix definition
                Vmath::Neg(m_ncoeffs, wsp, 1);
            }
            else
            {
                Vmath::Smul(m_ncoeffs,-1.0,inarray,1,wsp,1);
            }
           
            int bndcnt = 0;
            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            // Add weak boundary conditions to forcing
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eNeumann ||
                   m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eRobin)
                {
                    const Array<OneD, NekDouble> bndcoeff =
                        (m_bndCondExpansions[i])->GetCoeffs(); 

                    if(m_locToGloMap->GetSignChange())
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt + j]] += sign[bndcnt + j] * bndcoeff[j]; 
                        }
                    }
                    else
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt+j]] += bndcoeff[j]; 
                        }
                    }                    
                }
                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }

            GlobalLinSysKey key(StdRegions::eHelmholtz,m_locToGloMap,factors,
                                varcoeff,varfactors);
            
            GlobalSolve(key,wsp,outarray,dirForcing);
        }


        /**
         * This is equivalent to the operation:
         * \f[\boldsymbol{M\hat{u}}_g\f]
         * where \f$\boldsymbol{M}\f$ is the global matrix of type specified by
         * \a mkey. After scattering the global array \a inarray to local
         * level, this operation is evaluated locally by the function
         * ExpList#GeneralMatrixOp. The global result is then obtained by a
         * global assembly procedure.
         *
         * @param   mkey        This key uniquely defines the type matrix
         *                      required for the operation.
         * @param   inarray     The vector \f$\boldsymbol{\hat{u}}_g\f$ of size
         *                      \f$N_{\mathrm{dof}}\f$.
         * @param   outarray    The resulting vector of size
         *                      \f$N_{\mathrm{dof}}\f$.
         */
        void ContField2D::v_GeneralMatrixOp(
                const GlobalMatrixKey              &gkey,
                const Array<OneD,const NekDouble>  &inarray,
                Array<OneD,      NekDouble>  &outarray)
        {
            GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
        }

        /**
         * First compute the inner product of forcing function with respect to
         * base, and then solve the system with the linear advection operator.
         * @param   velocity    Array of advection velocities in physical space
         * @param   inarray     Forcing function.
         * @param   outarray    Result.
         * @param   lambda      reaction coefficient
         * @param   dirForcing  Dirichlet Forcing.
         */

        // could combine this with HelmholtzCG.
        void ContField2D::v_LinearAdvectionDiffusionReactionSolve(const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                       const Array<OneD, const NekDouble> &inarray,
                                                       Array<OneD, NekDouble> &outarray,
                                                       const NekDouble lambda,
                                                       const Array<OneD, const NekDouble>& dirForcing)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            IProductWRTBase(inarray,wsp);
            
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_ncoeffs, wsp, 1);

            // Forcing function with weak boundary conditions
            int i,j;
            int bndcnt=0;
            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            // Add weak boundary conditions to forcing
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eNeumann ||
                   m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eRobin)
                {
                    const Array<OneD, NekDouble> bndcoeff =
                        (m_bndCondExpansions[i])->GetCoeffs(); 

                    if(m_locToGloMap->GetSignChange())
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt + j]] += sign[bndcnt + j] * bndcoeff[j]; 
                        }
                    }
                    else
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            wsp[map[bndcnt+j]] += bndcoeff[bndcnt + j]; 
                        }
                    }                    
                }
                
                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }

            // Solve the system
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] = lambda;
            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffVelX] = velocity[0];
            varcoeffs[StdRegions::eVarCoeffVelY] = velocity[1];

            GlobalLinSysKey key(StdRegions::eLinearAdvectionDiffusionReaction,m_locToGloMap,factors,varcoeffs);

            GlobalSolve(key,wsp,outarray,dirForcing);
        }

        /**
         * First compute the inner product of forcing function with respect to
         * base, and then solve the system with the linear advection operator.
         * @param   velocity    Array of advection velocities in physical space
         * @param   inarray     Forcing function.
         * @param   outarray    Result.
         * @param   lambda      reaction coefficient
         * @param   dirForcing  Dirichlet Forcing.
         */
        void ContField2D::v_LinearAdvectionReactionSolve(
                            const Array<OneD, Array<OneD, NekDouble> > &velocity,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray,
                            const NekDouble lambda,
                            const Array<OneD, const NekDouble>& dirForcing)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            IProductWRTBase(inarray,wsp);

            // Solve the system
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] = lambda;
            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffVelX] = velocity[0];
            varcoeffs[StdRegions::eVarCoeffVelY] = velocity[1];
            GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,m_locToGloMap,factors,varcoeffs);

            GlobalSolve(key,wsp,outarray,dirForcing);
        }


        /**
         *
         */
        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                ContField2D::v_GetBndConditions()
        {
            return GetBndConditions();
        }


        /**
         * Reset the GlobalLinSys Manager
         */
        void ContField2D::v_ClearGlobalLinSysManager(void)
        {
            m_globalLinSysManager.ClearManager("GlobalLinSys");
        }

    } // end of namespace
} //end of namespace
