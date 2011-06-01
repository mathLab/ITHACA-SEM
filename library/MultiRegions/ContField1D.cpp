///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.cpp
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
// Description: Field definition for 1D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ContField1D
         * As opposed to the class #ContExpList1D, the class #ContField1D is
         * able to incorporate the boundary conditions imposed to the problem
         * to be solved. Therefore, the class is equipped with three additional
         * data members:
         * - #m_bndCondExpansions
         * - #m_bndTypes
         * - #m_bndCondEquations
         *
         * The first data structure, #m_bndCondExpansions,
         * contains the point Expansion on the boundary,  #m_bndTypes
         * stores information about the type of boundary condition on the
         * different parts of the boundary while #m_bndCondEquations holds the
         * equation of the imposed boundary conditions.
         *
         * Furthermore, in case of Dirichlet boundary conditions,
         * this class is capable of lifting a known solution satisfying these
         * boundary conditions. If we denote the unknown solution by
         * \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet
         * boundary conditions by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$, the
         * expansion then can be decomposed as
         * \f[ u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
         * u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
         * \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
         * \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}
         *   -1}\hat{u}_n^{\mathcal{H}}
         * \Phi_n(\boldsymbol{x}_i).\f]
         * This lifting is accomplished by ordering the known global degrees of
         * freedom, prescribed by the Dirichlet boundary conditions, first in
         * the global array \f$\boldsymbol{\hat{u}}\f$, that is,
         * \f[\boldsymbol{\hat{u}}=\left[ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{u}}^{\mathcal{H}}
         * \end{array} \right].\f]
         * Such kind of expansions are also referred to as continuoous fields.
         * This class should be used when solving 2D problems using a standard
         * Galerkin approach.
         */

        /**
         * Constructs an empty 1D continuous field.
         */
        ContField1D::ContField1D():
            DisContField1D(),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs()
        {
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory
         * for the arrays #m_coeffs and #m_phys. Furthermore, it constructs the
         * mapping array (contained in #m_locToGloMap) for the transformation
         * between local elemental level and global level, it calculates the
         * total number global expansion coefficients \f$\hat{u}_n\f$ and
         * allocates memory for the array #m_contCoeffs.
         *
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   solnType    Type of global system to use.
         */
        ContField1D::ContField1D(LibUtilities::CommSharedPtr &pComm,
                                 SpatialDomains::MeshGraph1D &graph1D,
                                 const GlobalSysSolnType solnType):
            DisContField1D(pComm,graph1D,solnType,false),
            m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            ApplyGeomInfo(graph1D);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                ::AllocateSharedPtr(m_comm,m_ncoeffs,*this,solnType);


            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates
         * memory for the arrays #m_coeffs and #m_phys. Furthermore, it
         * constructs the mapping array (contained in #m_locToGloMap) for the
         * transformation between local elemental level and global level, it
         * calculates the total number global expansion coefficients
         * \f$\hat{u}_n\f$ and allocates memory for the array #m_contCoeffs.
         * The constructor also discretises the boundary conditions, specified
         * by the argument \a bcs, by expressing them in terms of the
         * coefficient of the expansion on the boundary.
         *
         * SpatialDomains#MeshGraph1D object, describing the domain, and a
         * SpatialDomains#BoundaryConditions object.
         * @param   graph1D     The domain description.
         * @param   bcs         Boundary conditions.
         * @param   bc_loc      ? (optional)
         */
        ContField1D::ContField1D(LibUtilities::CommSharedPtr &pComm,
                                 SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const int bc_loc,
                                 const GlobalSysSolnType solnType):
            DisContField1D(pComm,graph1D,bcs,bc_loc,solnType),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,
                                               bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),
                                periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                ::AllocateSharedPtr(m_comm,m_ncoeffs,*this,
                                    solnType,
                                    m_bndCondExpansions,
                                    m_bndConditions,
                                    periodicVertices);

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates
         * memory for the arrays #m_coeffs and #m_phys. Furthermore, it
         * constructs the mapping array (contained in #m_locToGloMap) for the
         * transformation between local elemental level and global level, it
         * calculates the total number global expansion coefficients
         * \f$\hat{u}_n\f$ and allocates memory for the array #m_contCoeffs.
         * The constructor also discretises the boundary conditions, specified
         * by the argument \a bcs, by expressing them in terms of the
         * coefficient of the expansion on  the boundary.
         *
         * @param   graph1D     A 1D mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   variable    An optional parameter to indicate for which
         *                      variable the field should be constructed.
         */
        ContField1D::ContField1D(LibUtilities::CommSharedPtr &pComm,
                                 SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const std::string variable,
                                 const GlobalSysSolnType solnType):
            DisContField1D(pComm,graph1D,bcs,variable,solnType),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                ::AllocateSharedPtr(m_comm,m_ncoeffs,*this,
                                    solnType,
                                    m_bndCondExpansions,
                                    m_bndConditions,
                                    periodicVertices);

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates
         * memory for the arrays #m_coeffs and #m_phys. Furthermore, it
         * constructs the mapping array (contained in #m_locToGloMap) for the
         * transformation between local elemental level and global level, it
         * calculates the total number global expansion coefficients
         * \f$\hat{u}_n\f$ and allocates memory for the array #m_contCoeffs.
         * The constructor also discretises the boundary conditions, specified
         * by the argument \a bcs, by expressing them in terms of the
         * coefficient of the expansion on  the boundary.
         *
         * @deprecated  Basis definition is now provided by the
         * SpatialDomains#MeshGraph1D object making this constructor deprecated.
         *
         * @param   Ba          BasisKey for defining expansions.
         * @param   graph1D     A 1D mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc      ? (optional)
         * @param   solnType    Type of solution to use. By default Direct
         *                      Static Condensation is used.
         */
        ContField1D::ContField1D(LibUtilities::CommSharedPtr &pComm,
                                 const LibUtilities::BasisKey &Ba,
                                 SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const int bc_loc,
                                 const GlobalSysSolnType solnType):
            DisContField1D(pComm,graph1D,solnType,false),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,
                                               bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),
                                periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                ::AllocateSharedPtr(m_comm,m_ncoeffs,*this,
                                    solnType,
                                    m_bndCondExpansions,
                                    m_bndConditions,
                                    periodicVertices);

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates
         * memory for the arrays #m_coeffs and #m_phys. Furthermore, it
         * constructs the mapping array (contained in #m_locToGloMap) for the
         * transformation between local elemental level and global level, it
         * calculates the total number global expansion coefficients
         * \f$\hat{u}_n\f$ and allocates memory for the array #m_contCoeffs.
         * The constructor also discretises the boundary conditions, specified
         * by the argument \a bcs, by expressing them in terms of the
         * coefficient of the expansion on  the boundary.
         *
         * @deprecated  Basis definition is now provided by the
         * SpatialDomains#MeshGraph1D object making this constructor deprecated.
         *
         * @param   Ba          BasisKey for defining expansions.
         * @param   graph1D     A 1D mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   variable    An optional parameter to indicate for which
         *                      variable the field should be constructed.
         * @param   solnType    Type of solution to use. By default Direct
         *                      Static Condensation is used.
         */

        ContField1D::ContField1D(LibUtilities::CommSharedPtr &pComm,
                                 const LibUtilities::BasisKey &Ba,
                                 SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const std::string variable,
                                 const GlobalSysSolnType solnType):
            DisContField1D(pComm,graph1D,solnType,false),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                                ::AllocateSharedPtr(m_comm, m_ncoeffs,*this,
                                                    solnType,
                                                    m_bndCondExpansions,
                                                    m_bndConditions,
                                                    periodicVertices);

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        /**
         * Constructs a 1D continuous field as a copy of an existing field
         * including all boundary conditions.
         * @param   In          Existing continuous field to duplicate.
         */
        ContField1D::ContField1D(const ContField1D &In):
            DisContField1D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_contNcoeffs(In.m_contNcoeffs),
            m_contCoeffs(m_contNcoeffs,0.0),
            m_globalLinSys(In.m_globalLinSys)
        {
        }



        /**
         *
         */
        ContField1D::~ContField1D()
        {
        }


        /**
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
         * @param   inarray     Discrete \f$f(x)\f$.
         * @param   outarray    Storage for result.
         * @param   UseContCoeffs   .
         */
        void ContField1D::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray,
                                   bool  UseContCoeffs)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);
            IProductWRTBase(inarray,wsp,true);

            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass,
                                m_locToGloMap,
                                m_locToGloMap->GetGlobalSysSolnType());

            if(UseContCoeffs)
            {
                GlobalSolve(key,wsp,outarray);
            }
            else
            {
                Array<OneD,NekDouble> tmp(m_contNcoeffs,0.0);
                GlobalSolve(key,wsp,tmp);
                GlobalToLocal(tmp,outarray);
            }
        }


        /**
         * Given the coefficients of an expansion, this function evaluates the
         * spectral/hp expansion \f$u^{\delta}(x)\f$ at the quadrature
         * points \f$x_i\f$. This operation is evaluated locally by the
         * function ExpList#BwdTrans.
         *
         * The coefficients of the expansion should be contained in the
         * variable #m_coeffs of the ExpList object \a In. The resulting
         * physical values at the quadrature points \f$u^{\delta}(x_i)\f$ are
         * stored in the array #m_phys.
         *
         * @param   In          An ExpList, containing the local
         *                      coefficients \f$\hat{u}_n^e\f$ in its array
         *                      #m_coeffs.
         */
        void ContField1D::BwdTrans(
                                const Array<OneD, const NekDouble>  &inarray,
                                      Array<OneD,       NekDouble>  &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                GlobalToLocal(inarray,wsp);
                BwdTrans_IterPerExp(wsp,outarray);
            }
            else
            {
                BwdTrans_IterPerExp(inarray,outarray);
            }
        }


        /**
         *
         */
        void ContField1D::MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            GlobalLinSysKey key(StdRegions::eMass,
                                m_locToGloMap,
                                m_locToGloMap->GetGlobalSysSolnType());
            if(UseContCoeffs)
            {
                if(inarray.data() == outarray.data())
                {
                    Array<OneD, NekDouble> tmp(m_contNcoeffs,0.0);
                    Vmath::Vcopy(m_contNcoeffs,inarray,1,tmp,1);
                    GlobalSolve(key,tmp,outarray);
                }
                else
                {
                    GlobalSolve(key,inarray,outarray);
                }
            }
            else
            {
                Array<OneD, NekDouble> globaltmp(m_contNcoeffs,0.0);

                if(inarray.data() == outarray.data())
                {
                    Array<OneD,NekDouble> tmp(inarray.num_elements());
                    Vmath::Vcopy(inarray.num_elements(),inarray,1,tmp,1);
                    Assemble(tmp,outarray);
                }
                else
                {
                    Assemble(inarray,outarray);
                }

                GlobalSolve(key,outarray,globaltmp);
                GlobalToLocal(globaltmp,outarray);
            }
        }


        /**
         * Given a linear system specified by the key \a key,
         * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}},\f]
         * this function solves this linear system taking into account the
         * boundary conditions specified in the data member
         * #m_bndCondExpansions.
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
         * \boldsymbol{M}^{\mathcal{HH}}\boldsymbol{\hat{u}}^{\mathcal{H}}
         * = \boldsymbol{\hat{f}}^{\mathcal{H}}
         *   +\boldsymbol{\hat{g}}^{\mathcal{H}}
         *   -\boldsymbol{M}^{\mathcal{HD}}\boldsymbol{\hat{u}}^{\mathcal{D}}\f]
         *
         * @param   key         Specifes the linear system to solve.
         * @param   rhs         Forcing term \f$\boldsymbol{f}\f$.
         * @param   inout       Solution vector \f$\boldsymbol{\hat{u}}\f$.
         * @param   dirForcing  .
         */
        void ContField1D::GlobalSolve(const GlobalLinSysKey &key,
                                const Array<OneD, const NekDouble>& rhs,
                                      Array<OneD,       NekDouble>& inout,
                                const Array<OneD, const NekDouble>& dirForcing)
        {
            int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();

            // STEP 1: SET THE DIRICHLET DOFS TO THE RIGHT VALUE
            //         IN THE SOLUTION ARRAY
            for(int i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    inout[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(i)]
                        = m_bndCondExpansions[i]->GetCoeff(0);
                }
            }

            // STEP 2: CALCULATE THE HOMOGENEOUS COEFFICIENTS
            if(m_contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                LinSys->Solve(rhs,inout,m_locToGloMap,dirForcing);
            }
        }


        /**
         * The function searches the map #m_globalLinSys to see if the global
         * matrix has been created before. If not, it calls the function
         * #GenglobalLinSys to generate the requested global system.
         *
         * @param   mkey        Key specifying the linear system.
         * @returns             Pointer to the required linear system.
         */
        GlobalLinSysSharedPtr ContField1D::GetGlobalLinSys(
                                const GlobalLinSysKey &mkey)
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalLinSys->find(mkey);

            if(matrixIter == m_globalLinSys->end())
            {
                glo_matrix = GenGlobalLinSys(mkey,m_locToGloMap);
                (*m_globalLinSys)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

        /**
         * The operation is evaluated locally (i.e. with respect to all local
         * expansion modes) by the function ExpList#IProductWRTBase. The inner
         * product with respect to the global expansion modes is than obtained
         * by a global assembly operation.
         *
         * The values of the function \f$f(x)\f$ evaluated at the quadrature
         * points \f$x_i\f$ should be contained in the variable #m_phys of the
         * ExpList object \a in. The result is stored in the array
         * #m_contCoeffs.
         *
         * @param   In          An ExpList, containing the discrete evaluation
         *                      of \f$f(x)\f$ at the quadrature points in its
         *                      array #m_phys.
         */
        void ContField1D::IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                IProductWRTBase_IterPerExp(inarray,wsp);
                Assemble(wsp,outarray);
            }
            else
            {
                IProductWRTBase_IterPerExp(inarray,outarray);
            }
        }
		
		int ContField1D::v_GetContNcoeffs() const
        {
            return m_contNcoeffs;
        }
		
		/**
         *
         */
        void ContField1D::v_SetContCoeffsArray(Array<OneD, NekDouble> &inarray)
        {
            m_contCoeffs = inarray;
        }

        void ContField1D::v_FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            FwdTrans(inarray,outarray,UseContCoeffs);
        }

        void ContField1D::v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);
        }

        void ContField1D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveCG(inarray, outarray, lambda, varLambda, varCoeff,
                                false, NullNekDouble1DArray);
        }

        /**
         * Consider the one dimensional Helmholtz equation,
         * \f[\frac{d^2u}{dx^2}-\lambda u(x) = f(x),\f]
         * supplemented with appropriate boundary conditions (which are
         * contained in the data member #m_bndCondExpansions). Applying a
         * \f$C^0\f$ continuous Galerkin discretisation, this equation leads to
         * the following linear system:
         * \f[\left( \boldsymbol{M}+\lambda\boldsymbol{L}\right)
         * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
         * where \f$\boldsymbol{M}\f$ and \f$\boldsymbol{L}\f$ are the mass and
         * Laplacian matrix respectively. This function solves the system above
         * for the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to
         * the function #GlobalSolve.
         *
         * The values of the function \f$f(x)\f$ evaluated at the
         * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a inarray. The resulting
         * global coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the
         * array #m_contCoeffs.
         *
         * @param   inarray     Input containing forcing function
         *                      \f$\boldsymbol{f}\f$ at the quadrature points.
         * @param   outarray    Output containing the coefficients
         *                      \f$\boldsymbol{u}_g\f$
         * @param   lambda      Parameter value.
         * @param   Sigma       Coefficients of lambda.
         * @param   varcoeff    Variable diffusivity coefficients.
         * @param   UseContCoeffs   Default: false
         * @param   dirForcing  Dirichlet Forcing.
         */
        void ContField1D::v_HelmSolveCG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          bool UseContCoeffs,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);
            IProductWRTBase(inarray,wsp,true);
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_contNcoeffs, wsp, 1);

            // Forcing function with weak boundary conditions
            int i;
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    wsp[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(i)]
                        += m_bndCondExpansions[i]->GetCoeff(0);
                }
            }

            // Solve the system
            GlobalLinSysKey key(StdRegions::eHelmholtz,
                                m_locToGloMap,lambda,
                                m_locToGloMap->GetGlobalSysSolnType());

            if(UseContCoeffs)
            {
                GlobalSolve(key,wsp,outarray,dirForcing);
            }
            else
            {
                Array<OneD,NekDouble> tmp(m_contNcoeffs,0.0);
                GlobalSolve(key,wsp,tmp,dirForcing);
                GlobalToLocal(tmp,outarray);
            }
        }

        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                ContField1D::v_GetBndConditions()
        {
            return GetBndConditions();
        }

        const Array<OneD, const NekDouble> &ContField1D::v_GetContCoeffs() const
        {
            return m_contCoeffs;
        }

        void ContField1D::v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            BwdTrans(inarray,outarray,UseContCoeffs);
        }

        void ContField1D::v_IProductWRTBase(
                                const Array<OneD, const NekDouble>   &inarray,
                                      Array<OneD,       NekDouble>   &outarray,
                                bool  UseContCoeffs)
        {
            IProductWRTBase(inarray,outarray,UseContCoeffs);
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
        void ContField1D::v_GeneralMatrixOp(
                                const GlobalMatrixKey                &gkey,
                                const Array<OneD,const NekDouble>    &inarray,
                                      Array<OneD,      NekDouble>    &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD,NekDouble> tmp1(2*m_ncoeffs);
                Array<OneD,NekDouble> tmp2(tmp1+m_ncoeffs);
                GlobalToLocal(inarray,tmp1);
                GeneralMatrixOp_IterPerExp(gkey,tmp1,tmp2);
                Assemble(tmp2,outarray);
            }
            else
            {
                GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
            }
        }

    } // end of namespace
} //end of namespace
