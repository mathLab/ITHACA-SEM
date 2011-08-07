///////////////////////////////////////////////////////////////////////////////
//
// File ContField3D.cpp
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
// Description: Field definition for 3D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D():
            DisContField3D(),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalMat()
        {
        }


        /**
         * Given a mesh \a graph3D, containing information about the domain and
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
         * @param   graph3D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   solnType    Type of global system to use.
         */
        ContField3D::ContField3D(LibUtilities::SessionReaderSharedPtr &pSession,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D):
            DisContField3D(pSession,graph3D,false),
            m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>
                ::AllocateSharedPtr(m_session,m_ncoeffs,*this);


            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        ContField3D::ContField3D(LibUtilities::SessionReaderSharedPtr &pSession,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const int bc_loc):
                DisContField3D(pSession,graph3D,bcs,bc_loc,false),
                m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
                m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_session, m_ncoeffs,*this,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces);

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
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
         * @param   In          Existing ContField2D object used to provide the
         *                      local to global mapping information and
         *                      global solution type.
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc
         */
        ContField3D::ContField3D(const ContField3D &In,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const int bc_loc):
            DisContField3D(In),
            m_globalMat   (MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            if(!SameTypeOfBoundaryConditions(In))
            {
                map<int,int> periodicFaces;
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);

                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                                                                         m_bndCondExpansions,
                                                                                         m_bndConditions,
                                                                                         periodicVertices,
                                                                                         periodicEdges,
                                                                                         periodicFaces);

            }
            else
            {
                m_locToGloMap = In.m_locToGloMap;
            }

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        ContField3D::ContField3D(LibUtilities::SessionReaderSharedPtr &pSession,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs,
                                 const std::string variable):
                DisContField3D(pSession,graph3D,bcs,variable,false),
                m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
                m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces);
            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        ContField3D::ContField3D(LibUtilities::SessionReaderSharedPtr &pSession,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string variable):
                DisContField3D(pSession,graph3D,variable,false),
                m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
                m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces);
            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
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
         * @param   In          Existing ContField2D object used to provide the
         *                      local to global mapping information and
         *                      global solution type.
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc
         */
        ContField3D::ContField3D(const ContField3D &In,
                                 SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string variable):
            DisContField3D(In),
            m_globalMat   (MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            if(!SameTypeOfBoundaryConditions(In))
            {
                map<int,int> periodicFaces;
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
                GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                                                                         m_bndCondExpansions,
                                                                                         m_bndConditions,
                                                                                         periodicVertices,
                                                                                         periodicEdges,
                                                                                         periodicFaces);

            }
            else
            {
                m_locToGloMap = In.m_locToGloMap;
            }

            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }



        ContField3D::ContField3D(const ContField3D &In):
                DisContField3D(In),
                m_locToGloMap(In.m_locToGloMap),
                m_contNcoeffs(In.m_contNcoeffs),
                m_contCoeffs(m_contNcoeffs,0.0),
                m_globalMat(In.m_globalMat),
                m_globalLinSys(In.m_globalLinSys)
        {
        }


        ContField3D::~ContField3D()
        {
        }


        /**
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * @param   In          ContField2D to compare with.
         * @returns True if boundary conditions match.
         */
        bool ContField3D::SameTypeOfBoundaryConditions(const ContField3D &In)
        {
            int i;
            bool returnval = true;

            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {

                // check to see if boundary condition type is the same
                // and there are the same number of boundary
                // conditions in the boundary definition.
                if((m_bndConditions[i]->GetBoundaryConditionType()
                    != In.m_bndConditions[i]->GetBoundaryConditionType())||
                   (m_bndCondExpansions[i]->GetExpSize()
                                    != In.m_bndCondExpansions[i]->GetExpSize()))
                {
                    returnval = false;
                    break;
                }
            }

            return returnval;
        }




        /**
         * Given the coefficients of an expansion, this function evaluates the
         * spectral/hp expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the
         * quadrature points \f$\boldsymbol{x}_i\f$. This operation is
         * evaluated locally by the function ExpList#BwdTrans.
         *
         * The coefficients of the expansion should be contained in the variable
         * #m_coeffs of the ExpList object \a In. The resulting physical values
         * at the quadrature points \f$u^{\delta}(\boldsymbol{x}_i)\f$ are
         * stored in the array #m_phys.
         *
         * @param   In          An ExpList, containing the local coefficients
         *                      \f$\hat{u}_n^e\f$ in its array #m_coeffs.
         */
        void ContField3D::v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(
                                                        StdRegions::eBwdTrans);

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


        /**
         * The operation is evaluated locally (i.e. with respect to all local
         * expansion modes) by the function ExpList#IProductWRTBase. The inner
         * product with respect to the global expansion modes is than obtained
         * by a global assembly operation.
         *
         * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the
         * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a in. The result is stored
         * in the array #m_contCoeffs.
         *
         * @param   In          An ExpList, containing the discrete evaluation
         *                      of \f$f(\boldsymbol{x})\f$ at the quadrature
         *                      points in its array #m_phys.
         */
        void ContField3D::v_IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD, NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(
                                                StdRegions::eIProductWRTBase);

                if(doGlobalOp)
                {
                    GlobalMatrixKey gkey(StdRegions::eIProductWRTBase,
                                         m_locToGloMap);
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


        void ContField3D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
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


        void ContField3D::v_MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray,
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


        void ContField3D::GenerateDirBndCondForcing(const GlobalLinSysKey &key,
                                                    Array<OneD, NekDouble> &inout,
                                                    Array<OneD, NekDouble> &outarray)
        {
            int bndcnt=0;
            const Array<OneD,const int>& map  = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
            NekDouble sign;

            for(int i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    const Array<OneD,const NekDouble>& coeffs = m_bndCondExpansions[i]->GetCoeffs();
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                        inout[map[bndcnt++]] = sign * coeffs[j];
                    }
                }
                else
                {
                    bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
                }
            }
            GeneralMatrixOp(*(key.GetGlobalMatrixKey()),inout,outarray,true);
        }

        // Note inout contains initial guess and final output.
        void ContField3D::GlobalSolve(const GlobalLinSysKey &key,
                                      const Array<OneD, const NekDouble>& rhs,
                                            Array<OneD,       NekDouble>& inout,
                                      const Array<OneD, const NekDouble>& dirForcing)
        {
            int i,j;
            int bndcnt=0;
            int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();

            // STEP 1: SET THE DIRICHLET DOFS TO THE RIGHT VALUE
            //         IN THE SOLUTION ARRAY
            const Array<OneD,const int>& map  = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
            NekDouble sign;

            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    const Array<OneD,const NekDouble>& coeffs = m_bndCondExpansions[i]->GetCoeffs();
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                        inout[map[bndcnt++]] = sign * coeffs[j];
                    }
                }
                else
                {
                    bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
                }
            }

            // STEP 2: CALCULATE THE HOMOGENEOUS COEFFICIENTS
            if(m_contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                LinSys->Solve(rhs,inout,m_locToGloMap,dirForcing);
            }
        }

        GlobalLinSysSharedPtr ContField3D::GetGlobalLinSys(const GlobalLinSysKey &mkey)
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
         * Returns the global matrix associated with the given GlobalMatrixKey.
         * If the global matrix has not yet been constructed on this field,
         * it is first constructed using GenGlobalMatrix().
         * @param   mkey        Global matrix key.
         * @returns Assocated global matrix.
         */
        GlobalMatrixSharedPtr ContField3D::GetGlobalMatrix(
                                const GlobalMatrixKey &mkey)
        {
            ASSERTL1(mkey.LocToGloMapIsDefined(),
                     "To use method must have a LocalToGlobalBaseMap "
                     "attached to key");

            GlobalMatrixSharedPtr glo_matrix;
            GlobalMatrixMap::iterator matrixIter = m_globalMat->find(mkey);

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


        void ContField3D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveCG(inarray, outarray, lambda, varLambda, varCoeff,
                              false, NullNekDouble1DArray);
        }

        // Solve the helmholtz problem assuming that m_contCoeff vector
        // contains an intial estimate for solution
        void ContField3D::v_HelmSolveCG(
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
            int i,j;
            int bndcnt = 0;
            NekDouble sign;
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                    {
                        sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                        wsp[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(bndcnt++)] +=
                            sign * (m_bndCondExpansions[i]->GetCoeffs())[j];
                    }
                }
                else
                {
                    bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
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

        void ContField3D::v_GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(
                                                        gkey.GetMatrixType());

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

        int ContField3D::GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
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

  } //end of namespace
} //end of namespace
