///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.h
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
// Description: Expansion list top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H
#define NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H

#include <MultiRegions/MultiRegions.hpp>
#include <StdRegions/StdExpansion.h>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalOptimizationParameters.h>

#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/SegGeom.h>

#include <SpatialDomains/MeshGraph.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class GlobalLinSys;
        class LocalToGlobalC0ContMap;
        class LocalToGlobalBaseMap;
        class LocalToGlobalDGMap;
        class ExpList1D;
        class ExpList1D;

        /// A map between global matrix keys and their associated block
        /// matrices.
        typedef map<GlobalMatrixKey,DNekScalBlkMatSharedPtr> BlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef boost::shared_ptr<BlockMatrixMap> BlockMatrixMapShPtr;

        /// Base class for all multi-elemental spectral/hp expansions.
        class ExpList
        {
        public:
            /// The default constructor.
            ExpList();

            /// The copy constructor.
            ExpList(const ExpList &in);

            /// The default destructor.
            virtual ~ExpList();

            /// Copy coefficients from concatenated list to expansion list.
            void PutCoeffsInToElmtExp(void);

            /// Copy coefficients from expansion list to concatenated list.
            void PutElmtExpInToCoeffs(void);

            /// Copy one elements coefficients from the concatenated list
            /// to the expansion list.
            void PutCoeffsInToElmtExp(int eid);

            /// Copy one elements coefficients from the expansion list to
            /// the concatenated list.
            void PutElmtExpInToCoeffs(int eid);

            /// Copy physical data from \a m_phys to expansion list.
            void PutPhysInToElmtExp(void);

            /// Copy physical data from given array to expansion list.
            void PutPhysInToElmtExp(Array<OneD, const NekDouble> &in);

            /// Copy expansion list physical data to given array.
            void PutElmtExpInToPhys(Array<OneD,NekDouble> &out);

            /// Copy expansion list physical data from one element to array.
            void PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out);

            /// Returns the total number of local degrees of freedom
            /// \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
            // inline
            int GetNcoeffs(void) const;

            /// Evaulates the maximum number of modes in the elemental basis
            /// order over all elements
            // inline
            int EvalBasisNumModesMax(void) const;

            /// Returns the vector of the number of modes in the elemental
            /// basis order over all elements.
            const Array<OneD,int> EvalBasisNumModesMaxPerExp(void) const;

            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            // inline
            int GetTotPoints(void) const;

            /// Returns the total number of quadrature points for eid's element
            /// \f$=Q_{\mathrm{tot}}\f$.
            // inline
            int GetTotPoints(const int eid) const;

            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            // inline
            int GetNpoints(void) const;

            /// Sets the transformed state #m_transState of the coefficient
            /// arrays.
            // inline
            void SetTransState(const TransState transState);

            /// This function returns the transformed state #m_transState of
            /// the coefficient arrays.
            // inline
            TransState GetTransState(void) const;

            /// Fills the array #m_phys
            // inline
            void SetPhys(const Array<OneD, const NekDouble> &inarray);

            /// This function manually sets whether the array of physical
            /// values \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is
            /// filled or not.
            // inline
            void SetPhysState(const bool physState);

            /// This function indicates whether the array of physical values
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is filled or
            /// not.
            // inline
            bool GetPhysState(void) const;

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            NekDouble PhysIntegral (void);

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            NekDouble PhysIntegral(const Array<OneD, const NekDouble> &inarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all \emph{local}
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            void   IProductWRTBase_IterPerExp(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            // inline
            void IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to the derivative (in
            /// direction \param dir) of all \emph{local} expansion modes
            /// \f$\phi_n^e(\boldsymbol{x})\f$.
            void   IProductWRTDerivBase(const int dir,
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            /// This function elementally evaluates the forward transformation
            /// of a function \f$u(\boldsymbol{x})\f$ onto the global
            /// spectral/hp expansion.
            void   FwdTrans_IterPerExp (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            // inline
            void FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// This function elementally mulplies the coefficient space of
            /// Sin my the elemental inverse of the mass matrix.
            void  MultiplyByElmtInvMass (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            // inline
            void MultiplyByInvMassMatrix(
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            // inline
            void HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                           NekDouble lambda,
                           NekDouble tau = 1);

            void HelmSolve(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD,       NekDouble> &outarray,
                           const Array<OneD, const Array<OneD, NekDouble> > &varcoeffs,
                           NekDouble lambda,
                           NekDouble tau = 1);

            void FwdTrans_BndConstrained(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);


            /// This function elementally evaluates the backward transformation
            /// of the global spectral/hp element expansion.
            void BwdTrans_IterPerExp (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            // inline
            void BwdTrans (const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray,
                             bool  UseContCoeffs = false);

            /// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            void GetCoords(Array<OneD, NekDouble> &coord_0,
                        Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                        Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

            /// This function calculates Surface Normal vector of a smooth
            /// manifold.
            void GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,
                                const int k);

            /// This function writes the spectral/hp element solution to the
            /// file \a out.
            void WriteToFile(std::ofstream &out,
                                OutputFormat format = eTecplot,
                                std::string var = "v");

            void ReadFromFile(std::ifstream &in,
                                OutputFormat format = eTecplot);

            /// This function returns the dimension of the coordinates of the
            /// element \a eid.
            // inline
            int GetCoordim(int eid);

            /// Set the \a i th coefficiient in \a m_coeffs to value \a val
            // inline
            void SetCoeff(int i, NekDouble val);

            /// Set the \a i th coefficiient in  #m_coeffs to value \a val
            // inline
            void SetCoeffs(int i, NekDouble val);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            // inline
            const Array<OneD, const NekDouble> &GetCoeffs() const;

            // inline
            const Array<OneD, const NekDouble> &GetContCoeffs() const;

            /// Get the \a i th value  (coefficient) of #m_coeffs
            // inline
            NekDouble GetCoeff(int i);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            // inline
            NekDouble GetCoeffs(int i);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            // inline
            const Array<OneD, const NekDouble> &GetPhys()  const;

            /// This function calculates the \f$L_\infty\f$ error of the global
            /// spectral/hp element approximation.
            NekDouble Linf (const Array<OneD, const NekDouble> &soln);

            /// This function calculates the \f$L_2\f$ error with
            /// respect to soln of the global
            /// spectral/hp element approximation.
            NekDouble L2 (const Array<OneD, const NekDouble> &soln);

            /// This function calculates the \f$L_2\f$ measure of the global
            /// spectral/hp element approximation.
            NekDouble L2 (void);

            /// Calculates the \f$H^1\f$ error of the global spectral/hp
            /// element approximation.
            NekDouble H1 (const Array<OneD, const NekDouble> &soln);

            /// This function returns the number of elements in the expansion.
            // inline
            int GetExpSize(void);

            /// This function returns (a shared pointer to) the local elemental
            /// expansion of the \f$n^{\mathrm{th}}\f$ element.
            // inline
            StdRegions::StdExpansionSharedPtr& GetExp(int n);

            /// Get the start offset position for a global list of #m_coeffs
            /// correspoinding to element n.
            // inline
            int GetCoeff_Offset(int n);

            /// Get the start offset position for a global list of m_phys
            /// correspoinding to element n.
            // inline
            int GetPhys_Offset(int n);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            // inline
            Array<OneD, NekDouble> &UpdateCoeffs();

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            // inline
            Array<OneD, NekDouble> &UpdatePhys();

            /// This function discretely evaluates the derivative of a function
            /// \f$f(\boldsymbol{x})\f$ on the domain consisting of all
            /// elements of the expansion.
            void PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                        Array<OneD, NekDouble> &out_d0,
                        Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                        Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            void PhysDeriv(const int dir,
                           const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &out_d);


            // functions associated with DisContField
            // inline
            const Array<OneD, const  boost::shared_ptr<ExpList1D> >
                                                &GetBndCondExpansions();

            // inline
            boost::shared_ptr<ExpList1D> &GetTrace();

            // inline
            boost::shared_ptr<LocalToGlobalDGMap> &GetTraceMap(void);

            // inline
            void AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fx,
                                const Array<OneD, const NekDouble> &Fy,
                                      Array<OneD, NekDouble> &outarray);

            // inline
            void AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fn,
                                      Array<OneD, NekDouble> &outarray);

            // inline
            void AddTraceBiIntegral(
                                const Array<OneD, const NekDouble> &Fwd,
                                const Array<OneD, const NekDouble> &Bwd,
                                      Array<OneD, NekDouble> &outarray);

            // inline
            void GetFwdBwdTracePhys( Array<OneD,NekDouble> &Fwd,
                                            Array<OneD,NekDouble> &Bwd);

            // inline
            void GetFwdBwdTracePhys(
                                const Array<OneD,const NekDouble> &field,
                                      Array<OneD,NekDouble> &Fwd,
                                      Array<OneD,NekDouble> &Bwd);

            // inline
            virtual void ExtractTracePhys(
                                Array<OneD,NekDouble> &outarray);

            // inline
            void ExtractTracePhys(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray);

            // inline
            const Array<OneD, const SpatialDomains
                                ::BoundaryConditionShPtr>& GetBndConditions();

            // inline
            void EvaluateBoundaryConditions(const NekDouble time = 0.0);


            // Routines for continous matrix solution
            /// This function calculates the result of the multiplication of a
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            // inline
            void GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            void GeneralMatrixOp_IterPerExp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray);


            std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                                                     GetFieldDefinitions(void);

            /// Append the element data listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(
                        SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                        std::vector<NekDouble> &fielddata);

            /// Extract the data in fielddata into the m_coeff list
            void ExtractDataToCoeffs(
                        SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                        std::vector<NekDouble> &fielddata,
                        std::string &field);

            /// load global optimisation parameters
            // inline
            void ReadGlobalOptimizationParameters(
                                const std::string &infilename);

            // inline
            void SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp);

            // inline
            void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                Array<OneD,int> &EdgeID);

        protected:

            boost::shared_ptr<DNekMat> GenGlobalMatrixFull(
                const GlobalLinSysKey &mkey,
                const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            void SetCoeffPhys(void);

            /// The total number of local degrees of freedom. #m_ncoeffs
            /// \f$=N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_l\f$
            int m_ncoeffs;

            /// The total number of quadrature points. #m_npoints
            /// \f$=Q_{\mathrm{tot}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_Q\f$
            int m_npoints;

            /**
             * \brief Concatenation of all local expansion coefficients.
             *
             * The array of length #m_ncoeffs\f$=N_{\mathrm{eof}}\f$ which is
             * the concatenation of the local expansion coefficients
             * \f$\hat{u}_n^e\f$ over all \f$N_{\mathrm{el}}\f$ elements
             * \f[\mathrm{\texttt{m\_coeffs}}=\boldsymbol{\hat{u}}_{l} =
             * \underline{\boldsymbol{\hat{u}}}^e = \left [ \begin{array}{c}
             * \boldsymbol{\hat{u}}^{1} \       \
             * \boldsymbol{\hat{u}}^{2} \       \
             * \vdots \                                                 \
             * \boldsymbol{\hat{u}}^{{{N_{\mathrm{el}}}}} \end{array} \right ],
             * \quad
             * \mathrm{where}\quad \boldsymbol{\hat{u}}^{e}[n]=\hat{u}_n^{e}\f]
             */
            Array<OneD, NekDouble> m_coeffs;

            /**
             * \brief The global expansion evaluated at the quadrature points
             *
             * The array of length #m_npoints\f$=Q_{\mathrm{tot}}\f$ containing
             * the evaluation of \f$u^{\delta}(\boldsymbol{x})\f$ at the
             * quadrature points \f$\boldsymbol{x}_i\f$.
             * \f[\mathrm{\texttt{m\_phys}}=\boldsymbol{u}_{l} =
             * \underline{\boldsymbol{u}}^e = \left [ \begin{array}{c}
             * \boldsymbol{u}^{1} \             \
             * \boldsymbol{u}^{2} \             \
             * \vdots \                                                 \
             * \boldsymbol{u}^{{{N_{\mathrm{el}}}}} \end{array} \right ],\quad
             * \mathrm{where}\quad
             * \boldsymbol{u}^{e}[i]=u^{\delta}(\boldsymbol{x}_i)\f]
             */
            Array<OneD, NekDouble> m_phys;

            /**
             * \brief The transformed state of the array of coefficients of the
             * expansion.
             *
             * #TransState is the enumeration which can attain the following
             * values:
             * - <em>eNotSet</em>: The coefficients are not set.
             * - <em>eLocal</em>: The array #m_coeffs is filled with the proper
             *   local coefficients.
             * - <em>eContinuous</em>: The array #m_contCoeffs is filled with
             *   the proper global coefficients.
             * - <em>eLocalCont</em>: Both the arrays #m_coeffs and
             *   #m_contCoeffs are filled with the proper coefficients.
             */
            TransState m_transState;

            /**
             * \brief The state of the array #m_phys.
             *
             * Indicates whether the array #m_phys, created to contain the
             * evaluation of \f$u^{\delta}(\boldsymbol{x})\f$ at the quadrature
             * points \f$\boldsymbol{x}_i\f$, is filled with these values.
             */
            bool       m_physState;

            /**
             * \brief The list of local expansions.
             *
             * The (shared pointer to the) vector containing (shared pointers
             * to) all local expansions. The fact that the local expansions are
             * all stored as a (pointer to a) #StdExpansion, the abstract base
             * class for all local expansions, allows a general implementation
             * where most of the routines for the derived classes are defined
             * in the #ExpList base class.
             */
            boost::shared_ptr<StdRegions::StdExpansionVector> m_exp;

            /// Offset of elemental data into the array #m_coeffs
            Array<OneD, int>  m_coeff_offset;

            /// Offset of elemental data into the array #m_phys
            Array<OneD, int>  m_phys_offset;

            NekOptimize::GlobalOptParamSharedPtr m_globalOptParam;

            BlockMatrixMapShPtr  m_blockMat;

            /// This function assembles the block diagonal matrix of local
            /// matrices of the type \a mtype.
            const DNekScalBlkMatSharedPtr GenBlockMatrix(
                                const GlobalMatrixKey &gkey);

            const DNekScalBlkMatSharedPtr& GetBlockMatrix(
                                const GlobalMatrixKey &gkey);

            void MultiplyByBlockMatrix(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray);

            boost::shared_ptr<GlobalMatrix>  GenGlobalMatrix(
                const GlobalMatrixKey &mkey,
                const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);


            /// This operation constructs the global linear system of type \a
            /// mkey.
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSysFullDirect(
                const GlobalLinSysKey &mkey,
                const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);


            void GlobalEigenSystem(const boost::shared_ptr<DNekMat> &Gmat,
                                   Array<OneD, NekDouble> &EigValsReal,
                                   Array<OneD, NekDouble> &EigValsImag,
                                   Array<OneD, NekDouble> &EigVecs
                                                    = NullNekDouble1DArray);

             /// This function constructs the necessary global matrices
             /// required for solving the linear system of type \a mkey by
             /// static condensation.
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSysStaticCond(
                const GlobalLinSysKey &mkey,
                const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /// This operation constructs the global linear system of type \a
            /// mkey.
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSys(
                const GlobalLinSysKey &mkey,
                const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /// Generate a GlobalLinSys from information provided by the key
            /// "mkey" and the mapping provided in LocToGloBaseMap.
            boost::shared_ptr<GlobalLinSys> GenGlobalBndLinSys(
                const GlobalLinSysKey     &mkey,
                const LocalToGlobalBaseMapSharedPtr &locToGloMap);

            // functions associated with DisContField
            virtual const
                Array<OneD,const boost::shared_ptr<ExpList1D> >
                &v_GetBndCondExpansions(void);

            virtual boost::shared_ptr<ExpList1D> &v_GetTrace();

            virtual boost::shared_ptr<LocalToGlobalDGMap>
                &v_GetTraceMap();

            virtual void v_AddTraceIntegral(
                                    const Array<OneD, const NekDouble> &Fx,
                                    const Array<OneD, const NekDouble> &Fy,
                                          Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceIntegral(
                                    const Array<OneD, const NekDouble> &Fn,
                                          Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceBiIntegral(
                                    const Array<OneD, const NekDouble> &Fwd,
                                    const Array<OneD, const NekDouble> &Bwd,
                                          Array<OneD, NekDouble> &outarray);

            virtual void v_GetFwdBwdTracePhys(
                                          Array<OneD,NekDouble> &Fwd,
                                          Array<OneD,NekDouble> &Bwd);

            virtual void v_GetFwdBwdTracePhys(
                                    const Array<OneD,const NekDouble>  &field,
                                          Array<OneD,NekDouble> &Fwd,
                                          Array<OneD,NekDouble> &Bwd);

            virtual void v_ExtractTracePhys(
                                          Array<OneD,NekDouble> &outarray);

            virtual void v_ExtractTracePhys(
                                    const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray);

            virtual void v_MultiplyByInvMassMatrix(
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);
            
            virtual void v_HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,       NekDouble> &outarray,
                                     NekDouble lambda,
                                     NekDouble tau);
            
            virtual void v_HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray,
                                   const Array<OneD, const Array<OneD, NekDouble> > &varcoeffs,
                                   NekDouble lambda,
                                   NekDouble tau);

            // wrapper functions about virtual functions
            virtual const
                Array<OneD, const NekDouble> &v_GetContCoeffs() const;

            virtual void v_BwdTrans(
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_FwdTrans(
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_IProductWRTBase(
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_GeneralMatrixOp(
                                    const GlobalMatrixKey             &gkey,
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

            virtual void v_SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp);

            virtual void v_GetBoundaryToElmtMap(
                                    Array<OneD, int> &ElmtID,
                                    Array<OneD,int> &EdgeID);

        private:
            virtual const
                Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions();

            virtual void v_EvaluateBoundaryConditions(
                                    const NekDouble time = 0.0);

        };

        /// Shared pointer to an ExpList object.
        typedef boost::shared_ptr<ExpList>      ExpListSharedPtr;
        /// An empty ExpList object.
        static ExpList NullExpList;


        // Inline routines follow.

        /**
         * Returns the total number of local degrees of freedom
         * \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
         */
        inline int ExpList::GetNcoeffs() const
        {
            return m_ncoeffs;
        }

        /**
         * Evaulates the maximum number of modes in the elemental basis
         * order over all elements
         */
        inline int ExpList::EvalBasisNumModesMax() const
        {
            int i;
            int returnval = 0;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval = max(returnval,
                                    (*m_exp)[i]->EvalBasisNumModesMax());
            }

            return returnval;
        }

        /**
         *
         */
        inline const Array<OneD,int> ExpList::EvalBasisNumModesMaxPerExp(void)
                                                                           const
        {
            int i;
            Array<OneD,int> returnval((*m_exp).size(),0);

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval[i] 
                    = max(returnval[i],(*m_exp)[i]->EvalBasisNumModesMax());
            }

            return returnval;
        }


        /**
         *
         */
        inline int ExpList::GetTotPoints() const
        {
            return m_npoints;
        }

        inline int ExpList::GetTotPoints(const int eid) const
        {
            return (*m_exp)[eid]->GetTotPoints();
        }

        /**
         *
         */
        inline int ExpList::GetNpoints() const
        {
            return m_npoints;
        }

        /**
         *
         */
        inline void ExpList::SetTransState(const TransState transState)
        {
            m_transState = transState;
        }

        /**
         *
         */
        inline TransState ExpList::GetTransState() const
        {
            return m_transState;
        }

        /**
         * This function fills the array \f$\boldsymbol{u}_l\f$, the evaluation
         * of the expansion at the quadrature points (implemented as #m_phys),
         * with the values of the array \a inarray.
         *
         * @param   inarray         The array containing the values where
         *                          #m_phys should be filled with.
         */
        inline void ExpList::SetPhys(
                                const Array<OneD, const NekDouble> &inarray)
        {
            Vmath::Vcopy(m_npoints,&inarray[0],1,&m_phys[0],1);
            m_physState = true;
        }


        /**
         * @param   physState       \a true (=filled) or \a false (=not filled).
         */
        inline void ExpList::SetPhysState(const bool physState)
        {
            m_physState = physState;
        }


        /**
         * @return  physState       \a true (=filled) or \a false (=not filled).
         */
        inline bool ExpList::GetPhysState() const
        {
            return m_physState;
        }

        /**
         *
         */
        inline void ExpList::IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            v_IProductWRTBase(inarray,outarray,UseContCoeffs);
        }


        /**
         *
         */
        inline void ExpList::FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            v_FwdTrans(inarray,outarray,UseContCoeffs);
        }


        /**
         *
         */
        inline void ExpList::MultiplyByInvMassMatrix(
                        const Array<OneD,const NekDouble> &inarray,
                              Array<OneD,      NekDouble> &outarray,
                        bool  UseContCoeffs)
        {
            v_MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);
        }

        /**
         *
         */
        inline void ExpList::HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       NekDouble tau)
        {
            v_HelmSolve(inarray,outarray,lambda,tau);
        }
        
        inline void ExpList::HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeffs,
                                       NekDouble lambda,
                                       NekDouble tau)
        {
            v_HelmSolve(inarray,outarray,varcoeffs,lambda,tau);
        }


        /**
         *
         */
        inline void ExpList::BwdTrans (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            v_BwdTrans(inarray,outarray,UseContCoeffs);
        }


        /**
         * @param   eid         The index of the element to be checked.
         * @return  The dimension of the coordinates of the specific element.
         */
        inline int ExpList::GetCoordim(int eid)
        {
            ASSERTL2(eid <= (*m_exp).size(),
                     "eid is larger than number of elements");
            return (*m_exp)[eid]->GetCoordim();
        }

        /**
         * @param   i           The index of m_coeffs to be set
         * @param   val         The value which m_coeffs[i] is to be set to.
         */
        inline void ExpList::SetCoeff(int i, NekDouble val)
        {
            m_coeffs[i] = val;
        }


        /**
         * @param   i           The index of #m_coeffs to be set.
         * @param   val         The value which #m_coeffs[i] is to be set to.
         */
        inline void ExpList::SetCoeffs(int i, NekDouble val)
        {
            m_coeffs[i] = val;
        }

        /**
         * As the function returns a constant reference to a
         * <em>const Array</em>, it is not possible to modify the
         * underlying data of the array #m_coeffs. In order to do
         * so, use the function #UpdateCoeffs instead.
         *
         * @return  (A constant reference to) the array #m_coeffs.
         */
        inline const Array<OneD, const NekDouble> &ExpList::GetCoeffs() const
        {
            return m_coeffs;
        }

        inline const Array<OneD, const NekDouble> &ExpList::GetContCoeffs()
                                                                        const
        {
            return v_GetContCoeffs();
        }

        /**
         * @param   i               The index of #m_coeffs to be returned
         * @return  The NekDouble held in #m_coeffs[i].
         */
        inline NekDouble ExpList::GetCoeff(int i)
        {
            return m_coeffs[i];
        }

        /**
         * @param   i               The index of #m_coeffs to be returned
         * @return  The NekDouble held in #m_coeffs[i].
         */
        inline NekDouble ExpList::GetCoeffs(int i)
        {
            return m_coeffs[i];
        }

        /**
         * As the function returns a constant reference to a
         * <em>const Array</em> it is not possible to modify the
         * underlying data of the array #m_phys. In order to do
         * so, use the function #UpdatePhys instead.
         *
         * @return  (A constant reference to) the array #m_phys.
         */
        inline const Array<OneD, const NekDouble> &ExpList::GetPhys()  const
        {
            return m_phys;
        }

        /**
         * @return  \f$N_{\mathrm{el}}\f$, the number of elements in the
         *          expansion.
         */
        inline int ExpList::GetExpSize(void)
        {
            return (*m_exp).size();
        }

        /**
         * @param   n               The index of the element concerned.
         *
         * @return  (A shared pointer to) the local expansion of the
         *          \f$n^{\mathrm{th}}\f$ element.
         */
        inline StdRegions::StdExpansionSharedPtr& ExpList::GetExp(int n)
        {
            return (*m_exp)[n];
        }

        /**
         *
         */
        inline int ExpList::GetCoeff_Offset(int n)
        {
            return m_coeff_offset[n];
        }

        /**
         *
         */
        inline int ExpList::GetPhys_Offset(int n)
        {
            return m_phys_offset[n];
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them, rather use the function #GetCoeffs insead.
         *
         * @return  (A reference to) the array #m_coeffs.
         */
        inline Array<OneD, NekDouble> &ExpList::UpdateCoeffs()
        {
            m_transState = eLocal;
            return m_coeffs;
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them,  rather use the function #GetPhys instead.
         *
         * @return  (A reference to) the array #m_phys.
         */
        inline Array<OneD, NekDouble> &ExpList::UpdatePhys()
        {
            m_physState = true;
            return m_phys;
        }


        // functions associated with DisContField
        inline const Array<OneD, const  boost::shared_ptr<ExpList1D> >
                                            &ExpList::GetBndCondExpansions()
        {
            return v_GetBndCondExpansions();
        }


        inline boost::shared_ptr<ExpList1D> &ExpList::GetTrace()
        {
            return v_GetTrace();
        }

        inline boost::shared_ptr<LocalToGlobalDGMap> &ExpList::GetTraceMap()
        {
            return v_GetTraceMap();
        }

        inline void ExpList::AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fx,
                                const Array<OneD, const NekDouble> &Fy,
                                      Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceIntegral(Fx,Fy,outarray);
        }

        inline void ExpList::AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fn,
                                      Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceIntegral(Fn,outarray);
        }

        inline void ExpList::AddTraceBiIntegral(
                                    const Array<OneD, const NekDouble> &Fwd,
                                    const Array<OneD, const NekDouble> &Bwd,
                                          Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceBiIntegral(Fwd,Bwd,outarray);
        }

        inline void ExpList::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
                                         Array<OneD,NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(Fwd,Bwd);
        }

        inline void ExpList::GetFwdBwdTracePhys(
                                const Array<OneD,const NekDouble>  &field,
                                      Array<OneD,NekDouble> &Fwd,
                                      Array<OneD,NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(field,Fwd,Bwd);
        }

        inline void ExpList::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {
            v_ExtractTracePhys(outarray);
        }


        inline void ExpList::ExtractTracePhys(
                            const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray)
        {
            v_ExtractTracePhys(inarray,outarray);
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                            &ExpList::GetBndConditions()
        {
            return v_GetBndConditions();
        }

        inline void ExpList::EvaluateBoundaryConditions(const NekDouble time)
        {
            v_EvaluateBoundaryConditions(time);
        }


        // Routines for continous matrix solution
        /**
         * This operation is equivalent to the evaluation of
         * \f$\underline{\boldsymbol{M}}^e\boldsymbol{\hat{u}}_l\f$, that is,
         * \f[ \left[
         * \begin{array}{cccc}
         * \boldsymbol{M}^1 & 0 & \hspace{3mm}0 \hspace{3mm}& 0 \\
         * 0 & \boldsymbol{M}^2 & 0 & 0 \\
         * 0 &  0 & \ddots &  0 \\
         * 0 &  0 & 0 & \boldsymbol{M}^{N_{\mathrm{el}}} \end{array} \right]
         *\left [ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{1} \\
         * \boldsymbol{\hat{u}}^{2} \\
         * \vdots \\
         * \boldsymbol{\hat{u}}^{{{N_{\mathrm{el}}}}} \end{array} \right ]\f]
         * where \f$\boldsymbol{M}^e\f$ are the local matrices of type
         * specified by the key \a mkey. The decoupling of the local matrices
         * allows for a local evaluation of the operation. However, rather than
         * a local matrix-vector multiplication, the local operations are
         * evaluated as implemented in the function
         * StdRegions#StdExpansion#GeneralMatrixOp.
         *
         * @param   mkey            This key uniquely defines the type matrix
         *                          required for the operation.
         * @param   inarray         The vector \f$\boldsymbol{\hat{u}}_l\f$ of
         *                          size \f$N_{\mathrm{eof}}\f$.
         * @param   outarray        The resulting vector of size
         *                          \f$N_{\mathrm{eof}}\f$.
         */
        inline void ExpList::GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            v_GeneralMatrixOp(gkey,inarray,outarray,UseContCoeffs);
        }


        // load global optimisation parameters
        inline void ExpList::ReadGlobalOptimizationParameters(
                                const std::string &infilename)
        {
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                                            ::AllocateSharedPtr(infilename);
        }

        inline void ExpList::SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {
            v_SetUpPhysNormals(locexp);
        }

        inline void ExpList::GetBoundaryToElmtMap( Array<OneD, int> &ElmtID,
                                            Array<OneD,int> &EdgeID)
        {
            v_GetBoundaryToElmtMap(ElmtID,EdgeID);
        }

  } //end of namespace
} //end of namespace

#endif // EXPLIST_H

/**
* $Log: ExpList.h,v $
* Revision 1.76  2009/11/07 17:15:17  sehunchun
* Add GetTotPoints(idx)
*
* Revision 1.75  2009/11/06 21:51:18  sherwin
* Added L2_DGDeriv method
*
* Revision 1.74  2009/11/04 20:30:15  cantwell
* Added documentation to ExpList and ExpList1D and tidied up code.
*
* Revision 1.73  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.72  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
* Revision 1.71  2009/10/25 18:54:38  sherwin
* Added H1 norm for error evaluation
*
* Revision 1.70  2009/10/22 16:46:15  cbiotto
* Adding function EvalBasisNumModesMaxPerExp
*
* Revision 1.69  2009/10/22 16:40:35  cbiotto
* *** empty log message ***
*
* Revision 1.68  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.67  2009/07/08 17:22:47  sehunchun
* Deleting GetTanBasis
*
* Revision 1.66  2009/07/08 11:13:54  sehunchun
* Adding GetSurfaceNormal function
*
* Revision 1.65  2009/07/07 16:36:45  sehunchun
* Adding AddTraceBiIntegral...
*
* Revision 1.64  2009/07/03 15:38:25  sehunchun
* Adding GetTanBasis function
*
* Revision 1.63  2009/05/14 14:26:41  pvos
* Updates to apply the dirichlet boundary condition forcing inside the static condensation algorithm
*
* Revision 1.62  2009/04/27 21:34:58  sherwin
* Modified WriteToField Method
*
* Revision 1.61  2009/04/27 15:02:04  pvos
* From h-to-p efficiently updates
*
* Revision 1.60  2009/04/22 22:32:10  sherwin
* Added in method to read dat file
*
* Revision 1.59  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.58  2009/04/03 20:33:57  sherwin
* Update for Eigenfunction evaluation
*
* Revision 1.57  2009/03/23 11:52:15  pvos
* NekMatrix updates
*
* Revision 1.56  2009/03/23 10:51:52  pvos
* Added BlockMatrix support
*
* Revision 1.55  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.54  2009/03/04 05:58:49  bnelson
* Fixed visual studio compile errors.
*
* Revision 1.53  2009/02/08 09:11:49  sherwin
* General updates to introduce multiple matrix definitions based on different boundary types
*
* Revision 1.52  2009/02/03 14:33:08  pvos
* Modifications for solvers with time-dependent dirichlet BC's
*
* Revision 1.51  2009/02/02 16:43:26  claes
* Added virtual functions for solver access
*
* Revision 1.50  2009/01/06 21:05:57  sherwin
* Added virtual function calls for BwdTrans, FwdTrans and IProductWRTBase from the class ExpList. Introduced _IterPerExp versions of these methods in ExpList.cpp§
*
* Revision 1.49  2008/11/19 10:52:55  pvos
* Changed MultiplyByInvMassMatrix + added some virtual functions
*
* Revision 1.48  2008/11/01 22:06:45  bnelson
* Fixed Visual Studio compile error.
*
* Revision 1.47  2008/10/29 22:46:35  sherwin
* Updates for const correctness and a few other bits
*
* Revision 1.46  2008/10/19 15:57:52  sherwin
* Added method EvalBasisNumModesMax
*
* Revision 1.45  2008/10/16 10:21:42  sherwin
* Updates to make methods consisten with AdvectionDiffusionReactionsSolver. Modified MultiplyByInvMassMatrix to take local or global arrays
*
* Revision 1.44  2008/10/04 20:04:26  sherwin
* Modifications for solver access
*
* Revision 1.43  2008/09/16 13:36:06  pvos
* Restructured the LocalToGlobalMap classes
*
* Revision 1.42  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.41  2008/07/29 22:27:33  sherwin
* Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
*
* Revision 1.40  2008/07/15 13:00:04  sherwin
* Updates for DG advection solver - not yet debugged
*
* Revision 1.39  2008/07/12 19:08:29  sherwin
* Modifications for DG advection routines
*
* Revision 1.38  2008/07/12 17:31:39  sherwin
* Added m_phys_offset and rename m_exp_offset to m_coeff_offset
*
* Revision 1.37  2008/07/11 15:48:32  pvos
* Added Advection classes
*
* Revision 1.36  2008/06/06 23:27:20  ehan
* Added doxygen documentation
*
* Revision 1.35  2008/06/05 15:06:58  pvos
* Added documentation
*
* Revision 1.34  2008/05/29 21:35:03  pvos
* Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
*
* Revision 1.33  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.32  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.31  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.30  2008/01/23 21:50:52  sherwin
* Update from EdgeComponents to SegGeoms
*
* Revision 1.29  2007/12/17 13:05:04  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.28  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.27  2007/11/20 16:27:16  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
* Revision 1.26  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.25  2007/09/03 19:58:31  jfrazier
* Formatting.
*
* Revision 1.24  2007/07/27 03:10:49  bnelson
* Fixed g++ compile error.
*
* Revision 1.23  2007/07/26 08:40:49  sherwin
* Update to use generalised i/o hooks in Helmholtz1D
*
* Revision 1.22  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.21  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.20  2007/07/17 07:11:05  sherwin
* Chaned definition of NullExpList
*
* Revision 1.19  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.18  2007/07/13 09:02:24  sherwin
* Mods for Helmholtz solver
*
* Revision 1.17  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.16  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
