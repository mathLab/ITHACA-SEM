///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSys::lookupIds[12] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectFull",
                MultiRegions::eDirectFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectStaticCond",
                MultiRegions::eDirectStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "DirectMultiLevelStaticCond",
                MultiRegions::eDirectMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeFull",
                MultiRegions::eIterativeFull),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeStaticCond",
                MultiRegions::eIterativeStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "IterativeMultiLevelStaticCond",
                MultiRegions::eIterativeMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtFull",
                MultiRegions::eXxtFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtStaticCond",
                MultiRegions::eXxtStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "XxtMultiLevelStaticCond",
                MultiRegions::eXxtMultiLevelStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScFull",
                MultiRegions::ePETScFullMatrix),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScStaticCond",
                MultiRegions::ePETScStaticCond),
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalSysSoln", "PETScMultiLevelStaticCond",
                MultiRegions::ePETScMultiLevelStaticCond)
        };

#ifdef NEKTAR_USE_SCOTCH
        std::string GlobalLinSys::def = LibUtilities::SessionReader::
            RegisterDefaultSolverInfo("GlobalSysSoln",
                                      "DirectMultiLevelStaticCond");
#else
        std::string GlobalLinSys::def = LibUtilities::SessionReader::
            RegisterDefaultSolverInfo("GlobalSysSoln",
                                      "DirectStaticCond");
#endif

        /**
         * @class GlobalLinSys
         *
         * Consider the linear system
         * \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{\hat{f}}\f$.
         * Distinguishing between the boundary and interior components of
         * \f$\boldsymbol{\hat{u}}_g\f$ and \f$\boldsymbol{\hat{f}}\f$ using
         * \f$\boldsymbol{\hat{u}}_b\f$,\f$\boldsymbol{\hat{u}}_i\f$ and
         * \f$\boldsymbol{\hat{f}}_b\f$,\f$\boldsymbol{\hat{f}}_i\f$
         * respectively, this system can be split into its constituent parts as
         * \f[\left[\begin{array}{cc}
         * \boldsymbol{M}_b&\boldsymbol{M}_{c1}\\
         * \boldsymbol{M}_{c2}&\boldsymbol{M}_i\\
         * \end{array}\right]
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{u}_b}\\
         * \boldsymbol{\hat{u}_i}\\
         * \end{array}\right]=
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{f}_b}\\
         * \boldsymbol{\hat{f}_i}\\
         * \end{array}\right]\f]
         * where \f$\boldsymbol{M}_b\f$ represents the components of
         * \f$\boldsymbol{M}\f$ resulting from boundary-boundary mode
         * interactions,
         * \f$\boldsymbol{M}_{c1}\f$ and \f$\boldsymbol{M}_{c2}\f$ represent the
         * components resulting from coupling between the boundary-interior
         * modes, and \f$\boldsymbol{M}_i\f$ represents the components of
         * \f$\boldsymbol{M}\f$ resulting from interior-interior mode
         * interactions.
         *
         * The solution of the linear system can now be determined in two steps:
         * \f{eqnarray*}
         * \mathrm{step 1:}&\quad&(\boldsymbol{M}_b-\boldsymbol{M}_{c1}
         * \boldsymbol{M}_i^{-1}\boldsymbol{M}_{c2}) \boldsymbol{\hat{u}_b} =
         * \boldsymbol{\hat{f}}_b - \boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
         * \boldsymbol{\hat{f}}_i,\nonumber \\
         * \mathrm{step 2:}&\quad&\boldsymbol{\hat{u}_i}=\boldsymbol{M}_i^{-1}
         * \left( \boldsymbol{\hat{f}}_i
         *      - \boldsymbol{M}_{c2}\boldsymbol{\hat{u}_b}
         * \right). \nonumber \\ \f}
         * As the inverse of \f$\boldsymbol{M}_i^{-1}\f$ is
         * \f[ \boldsymbol{M}_i^{-1} = \left [\underline{\boldsymbol{M}^e_i}
         * \right ]^{-1} = \underline{[\boldsymbol{M}^e_i]}^{-1} \f]
         * and the following operations can be evaluated as,
         * \f{eqnarray*}
         * \boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}\boldsymbol{\hat{f}}_i &
         * =& \boldsymbol{\mathcal{A}}_b^T \underline{\boldsymbol{M}^e_{c1}}
         * \underline{[\boldsymbol{M}^e_i]}^{-1} \boldsymbol{\hat{f}}_i \\
         * \boldsymbol{M}_{c2} \boldsymbol{\hat{u}_b} &=&
         * \underline{\boldsymbol{M}^e_{c2}} \boldsymbol{\mathcal{A}}_b
         * \boldsymbol{\hat{u}_b}.\f}
         * where \f$\boldsymbol{\mathcal{A}}_b \f$ is the permutation matrix
         * which scatters from global to local degrees of freedom, only the
         * following four matrices should be constructed:
         * - \f$\underline{[\boldsymbol{M}^e_i]}^{-1}\f$
         * - \f$\underline{\boldsymbol{M}^e_{c1}}
         *                          \underline{[\boldsymbol{M}^e_i]}^{-1}\f$
         * - \f$\underline{\boldsymbol{M}^e_{c2}}\f$
         * - The Schur complement: \f$\boldsymbol{M}_{\mathrm{Schur}}=
         *   \quad\boldsymbol{M}_b-\boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
         *   \boldsymbol{M}_{c2}\f$
         *
         * The first three matrices are just a concatenation of the
         * corresponding local matrices and they can be created as such. They
         * also allow for an elemental evaluation of the operations concerned.
         *
         * The global Schur complement however should be assembled from the
         * concatenation of the local elemental Schur complements, that is,
         * \f[ \boldsymbol{M}_{\mathrm{Schur}}=\boldsymbol{M}_b
         *          - \boldsymbol{M}_{c1}
         * \boldsymbol{M}_i^{-1} \boldsymbol{M}_{c2} =
         * \boldsymbol{\mathcal{A}}_b^T \left [\underline{\boldsymbol{M}^e_b -
         * \boldsymbol{M}^e_{c1} [\boldsymbol{M}^e_i]^{-1}
         * (\boldsymbol{M}^e_{c2})} \right ] \boldsymbol{\mathcal{A}}_b \f]
         * and it is the only matrix operation that need to be evaluated on a
         * global level when using static condensation.
         * However, due to the size and sparsity of the matrix
         * \f$\boldsymbol{\mathcal{A}}_b\f$, it is more efficient to assemble
         * the global Schur matrix using the mapping array bmap\f$[e][i]\f$
         * contained in the input argument \a locToGloMap. The global Schur
         * complement is then constructed as:
         * \f[\boldsymbol{M}_{\mathrm{Schur}}\left[\mathrm{\a bmap}[e][i]\right]
         * \left[\mathrm{\a bmap}[e][j]\right]=\mathrm{\a bsign}[e][i]\cdot
         * \mathrm{\a bsign}[e][j]
         * \cdot\boldsymbol{M}^e_{\mathrm{Schur}}[i][j]\f]
         * All four matrices are stored in the \a GlobalLinSys returned by this
         * function.
         */

        /**
         * Given a block matrix, construct a global matrix system according to
         * a local to global mapping. #m_linSys is constructed by
         * AssembleFullMatrix().
         * @param   pkey        Associated linear system key.
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSys::GlobalLinSys(const GlobalLinSysKey &pKey,
                const std::weak_ptr<ExpList> &pExpList,
                const std::shared_ptr<AssemblyMap>
                                   &pLocToGloMap):
            m_linSysKey(pKey),
            m_expList(pExpList),
            m_robinBCInfo(m_expList.lock()->GetRobinBCInfo()),
            m_verbose(m_expList.lock()->GetSession()->
                      DefinesCmdLineArgument("verbose"))
        {
            boost::ignore_unused(pLocToGloMap);
        }

        /**
         *
         */
        GlobalLinSysFactory& GetGlobalLinSysFactory()
        {
            static GlobalLinSysFactory instance;
            return instance;
        }

        /**
         * @brief Create a preconditioner object from the parameters defined in
         * the supplied assembly map.
         *
         * @param asmMap  Assembly map used to construct the global system.
         */
        PreconditionerSharedPtr GlobalLinSys::CreatePrecon(AssemblyMapSharedPtr
                                                           asmMap)
        {
            PreconditionerType pType = asmMap->GetPreconType();
            std::string PreconType = MultiRegions::PreconditionerTypeMap[pType];
            return GetPreconFactory().CreateInstance(
                PreconType, GetSharedThisPtr(), asmMap);
        }


        /**
         * @brief Get the number of blocks in this system. 
         *
         * At the top level this corresponds to the number of elements in the
         * expansion list.
         */
        int GlobalLinSys::v_GetNumBlocks()
        {
            return m_expList.lock()->GetExpSize();
        }


        /** 
            Assemble the matrix key for each block n 
        **/
        
        LocalRegions::MatrixKey GlobalLinSys::GetBlockMatrixKey(unsigned int n)
        {            
            
            std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            int cnt = 0;
            
            LocalRegions::ExpansionSharedPtr vExp = expList->GetExp( n );
            
            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            StdRegions::ConstFactorMap vConstFactorMap =
                m_linSysKey.GetConstFactors();

            // setup variable factors
            if(m_linSysKey.GetNVarFactors() > 0)
            {
                if(m_linSysKey.GetVarFactors().
                   count(StdRegions::eFactorSVVDiffCoeff) != 0)
                {
                    vConstFactorMap[StdRegions::eFactorSVVDiffCoeff] =
                        m_linSysKey.GetVarFactors(
                                      StdRegions::eFactorSVVDiffCoeff)[n];
 
                    ASSERTL1(m_linSysKey.GetConstFactors().
                             count(StdRegions::eFactorSVVCutoffRatio),
                             "VarCoeffSVVCuroffRatio is set but "
                             " not FactorSVVCutoffRatio");

                    vConstFactorMap[StdRegions::eFactorSVVCutoffRatio] =
                        m_linSysKey.GetVarFactors(
                                   StdRegions::eFactorSVVCutoffRatio)[n];
                    
                }
                
                if(m_linSysKey.GetVarFactors().
                   count(StdRegions::eFactorSVVPowerKerDiffCoeff) != 0)
                {
                    vConstFactorMap[StdRegions::eFactorSVVPowerKerDiffCoeff] =
                        m_linSysKey.GetVarFactors(
                               StdRegions::eFactorSVVPowerKerDiffCoeff)[n];
                }

                if(m_linSysKey.GetVarFactors().
                   count(StdRegions::eFactorSVVDGKerDiffCoeff) != 0)
                {
                    vConstFactorMap[StdRegions::eFactorSVVDGKerDiffCoeff] =
                        m_linSysKey.GetVarFactors(
                                  StdRegions::eFactorSVVDGKerDiffCoeff)[n];
                }
            }

            // retrieve variable coefficients
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                cnt = expList->GetPhys_Offset(n);
                
                for (auto &x : m_linSysKey.GetVarCoeffs())
                {
                    vVarCoeffMap[x.first] = x.second + cnt;
                }
            }
            

            LocalRegions::MatrixKey matkey(m_linSysKey.GetMatrixType(),
                                           vExp->DetShapeType(),
                                           *vExp,
                                           vConstFactorMap,
                                           vVarCoeffMap);
            return matkey;
        }

        /**
         * @brief Retrieves the block matrix from n-th expansion using the
         * matrix key provided by the #m_linSysKey.
         *
         * @param   n           Number of the expansion.
         * @return              Block matrix for the specified expansion.
         */
        DNekScalMatSharedPtr GlobalLinSys::v_GetBlock(unsigned int n)
        {
            LocalRegions::ExpansionSharedPtr vExp = m_expList.lock()->GetExp( n );
            DNekScalMatSharedPtr loc_mat;
            loc_mat = vExp->GetLocMatrix(GetBlockMatrixKey(n));

            // apply robin boundary conditions to the matrix.
            if(m_robinBCInfo.count(n) != 0) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;

                // declare local matrix from scaled matrix.
                int rows = loc_mat->GetRows();
                int cols = loc_mat->GetColumns();
                const NekDouble *dat = loc_mat->GetRawPtr();
                DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(rows,cols,dat);
                Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                // add local matrix contribution
                for(rBC = m_robinBCInfo.find(n)->second;rBC; rBC = rBC->next)
                {
                    vExp->AddRobinMassMatrix(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs, new_mat);
                }

                // redeclare loc_mat to point to new_mat plus the scalar.
                loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    1.0, new_mat);
            }

            // finally return the matrix.
            return loc_mat;
        }

        /**
         * @brief Retrieves a the static condensation block matrices from n-th
         * expansion using the matrix key provided by the #m_linSysKey.
         * 
         * @param   n           Number of the expansion
         * @return              2x2 Block matrix holding the static condensation
         *                      matrices for the n-th expansion.
         */
        DNekScalBlkMatSharedPtr GlobalLinSys::v_GetStaticCondBlock(
            unsigned int n)
        {
            
            LocalRegions::ExpansionSharedPtr vExp = m_expList.lock()->GetExp( n );
            DNekScalBlkMatSharedPtr loc_mat;
            loc_mat = vExp->GetLocStaticCondMatrix(GetBlockMatrixKey(n));

            if(m_robinBCInfo.count(n) != 0) // add robin mass matrix
            {
                DNekScalMatSharedPtr    tmp_mat;
                RobinBCInfoSharedPtr rBC;

                tmp_mat = loc_mat->GetBlock(0,0);

                // declare local matrix from scaled matrix.
                int rows = tmp_mat->GetRows();
                int cols = tmp_mat->GetColumns();
                const NekDouble *dat = tmp_mat->GetRawPtr();
                DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(rows, cols, dat);
                Blas::Dscal(rows*cols,tmp_mat->Scale(),new_mat->GetRawPtr(),1);

                // add local matrix contribution
                for(rBC = m_robinBCInfo.find(n)->second;rBC; rBC = rBC->next)
                {
                    vExp->AddRobinMassMatrix(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs, new_mat);
                }

                // redeclare loc_mat to point to new_mat plus the scalar.
                tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    1.0, new_mat);
                DNekScalBlkMatSharedPtr new_loc_mat;
                unsigned int exp_size[] = {tmp_mat->GetRows(), loc_mat->GetBlock(1,1)->GetRows()};
                unsigned int nblks = 2;
                new_loc_mat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks, nblks, exp_size, exp_size);


                new_loc_mat->SetBlock(0,0,tmp_mat);
                new_loc_mat->SetBlock(0,1,loc_mat->GetBlock(0,1));
                new_loc_mat->SetBlock(1,0,loc_mat->GetBlock(1,0));
                new_loc_mat->SetBlock(1,1,loc_mat->GetBlock(1,1));
                loc_mat = new_loc_mat;
            }

            return loc_mat;
        }

        /**
         * @brief Releases the static condensation block matrices from NekManager
         * of n-th expansion using the matrix key provided by the #m_linSysKey.
         * 
         * @param   n           Number of the expansion
         */
        void GlobalLinSys::v_DropStaticCondBlock(unsigned int n)
        {
            LocalRegions::ExpansionSharedPtr vExp = m_expList.lock()->GetExp( n );
            vExp->DropLocStaticCondMatrix(GetBlockMatrixKey(n));
        }

        void GlobalLinSys::v_InitObject()
        {
            NEKERROR(ErrorUtil::efatal, "Method does not exist" );
        }

        void GlobalLinSys::v_Initialise(
            const std::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            boost::ignore_unused(pLocToGloMap);
            NEKERROR(ErrorUtil::efatal, "Method does not exist" );
        }
    } //end of namespace
} //end of namespace

