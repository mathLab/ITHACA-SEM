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

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Communication/Transposition.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/Expansion.h>
#include <Collections/Collection.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalOptimizationParameters.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <tinyxml.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class GlobalLinSys;
        class AssemblyMapDG;

        class AssemblyMapCG;
        class GlobalLinSysKey;
        class GlobalMatrix;

        enum Direction
        {
            eX,
            eY,
            eZ,
            eS,
            eN
        };

        enum ExpansionType
        {
            e0D,
            e1D,
            e2D,
            e2DH1D,
            e3DH1D,
            e3DH2D,
            e3D,
            eNoType
        };

        MultiRegions::Direction const DirCartesianMap[] =
        {
            eX,
            eY,
            eZ
        };

        /// A map between global matrix keys and their associated block
        /// matrices.
        typedef std::map<GlobalMatrixKey,DNekScalBlkMatSharedPtr> BlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef std::shared_ptr<BlockMatrixMap> BlockMatrixMapShPtr;


        /// Base class for all multi-elemental spectral/hp expansions.
        class ExpList: public std::enable_shared_from_this<ExpList>
        {
        public:
            /// The default constructor.
            MULTI_REGIONS_EXPORT ExpList();

            /// The default constructor.
            MULTI_REGIONS_EXPORT ExpList(
                    const LibUtilities::SessionReaderSharedPtr &pSession);

            /// The default constructor.
            MULTI_REGIONS_EXPORT ExpList(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const SpatialDomains::MeshGraphSharedPtr &pGraph);

            /// Constructor copying only elements defined in eIds.
            MULTI_REGIONS_EXPORT ExpList(
                const ExpList &in,
                const std::vector<unsigned int> &eIDs,
                const bool DeclareCoeffPhysArrays = true);

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList(
                const ExpList &in,
                const bool DeclareCoeffPhysArrays = true);

            /// The default destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList();

            /// Returns the total number of local degrees of freedom
            /// \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
            inline int GetNcoeffs(void) const;

            /// Returns the total number of local degrees of freedom
            /// for element eid
            MULTI_REGIONS_EXPORT int GetNcoeffs(const int eid) const;

            /// Returns the type of the expansion
            MULTI_REGIONS_EXPORT ExpansionType GetExpType(void);

            /// Returns the type of the expansion
            MULTI_REGIONS_EXPORT void SetExpType(ExpansionType Type);

            /// Evaulates the maximum number of modes in the elemental basis
            /// order over all elements
            inline int EvalBasisNumModesMax(void) const;

            /// Returns the vector of the number of modes in the elemental
            /// basis order over all elements.
            MULTI_REGIONS_EXPORT const Array<OneD,int>
                EvalBasisNumModesMaxPerExp(void) const;
            
            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetTotPoints(void) const;

            /// Returns the total number of quadrature points for eid's element
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetTotPoints(const int eid) const;

            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetNpoints(void) const;


            /// Returns the total number of qudature points scaled by
            /// the factor scale on each 1D direction
            inline int Get1DScaledTotPoints(const NekDouble scale) const;

            /// Sets the wave space to the one of the possible configuration
            /// true or false
            inline void SetWaveSpace(const bool wavespace);


            ///Set Modified Basis for the stability analysis
            inline void SetModifiedBasis(const bool modbasis);

            /// Set the \a i th value of \a m_phys to value \a val
            inline void SetPhys(int i, NekDouble val);

            /// This function returns the third direction expansion condition,
            /// which can be in wave space (coefficient) or not
            /// It is stored in the variable m_WaveSpace.
            inline bool GetWaveSpace(void) const;

            /// Fills the array #m_phys
            inline void SetPhys(const Array<OneD, const NekDouble> &inarray);

            /// Sets the array #m_phys
            inline void SetPhysArray(Array<OneD, NekDouble> &inarray);

            /// This function manually sets whether the array of physical
            /// values \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is
            /// filled or not.
            inline void SetPhysState(const bool physState);

            /// This function indicates whether the array of physical values
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is filled or
            /// not.
            inline bool GetPhysState(void) const;

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            MULTI_REGIONS_EXPORT NekDouble PhysIntegral (void);

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            MULTI_REGIONS_EXPORT NekDouble PhysIntegral(
                const Array<OneD,
                const NekDouble> &inarray);

            /// multiply the metric jacobi and quadrature weights
            MULTI_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                const Array<OneD, const NekDouble>  &inarray,
                Array<OneD, NekDouble>              &outarray);

            /// Divided by the metric jacobi and quadrature weights
            MULTI_REGIONS_EXPORT void DivideByQuadratureMetric(
                const Array<OneD, const NekDouble>  &inarray,
                Array<OneD, NekDouble>              &outarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all \em local
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            inline void   IProductWRTBase_IterPerExp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,       NekDouble> &outarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to the derivative (in
            /// direction \param dir) of all \em local expansion modes
            /// \f$\phi_n^e(\boldsymbol{x})\f$.
            MULTI_REGIONS_EXPORT void   IProductWRTDerivBase(
                const int dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            MULTI_REGIONS_EXPORT void   IProductWRTDirectionalDerivBase(
                const Array<OneD, const NekDouble> &direction,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to the derivative (in
            /// direction \param dir) of all \em local expansion modes
            /// \f$\phi_n^e(\boldsymbol{x})\f$.
            MULTI_REGIONS_EXPORT void   IProductWRTDerivBase
                (const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                 Array<OneD,       NekDouble> &outarray);

            /// This function elementally evaluates the forward transformation
            /// of a function \f$u(\boldsymbol{x})\f$ onto the global
            /// spectral/hp expansion.
            inline void  FwdTrans_IterPerExp (
                const Array<OneD,
                const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            ///
            inline void FwdTrans(
                const Array<OneD,
                const NekDouble> &inarray,
                Array<OneD,       NekDouble> &outarray);

            MULTI_REGIONS_EXPORT void   ExponentialFilter(
                Array<OneD, NekDouble> &array,
                const NekDouble        alpha,
                const NekDouble        exponent,
                const NekDouble        cutoff);

            /// This function elementally mulplies the coefficient space of
            /// Sin my the elemental inverse of the mass matrix.
            MULTI_REGIONS_EXPORT void  MultiplyByElmtInvMass (
                 const Array<OneD,
                 const NekDouble> &inarray,
                 Array<OneD,       NekDouble> &outarray);

            ///
            inline void MultiplyByInvMassMatrix(
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);

            /// Smooth a field across elements
            inline void SmoothField(Array<OneD,NekDouble> &field);

            /// Solve helmholtz problem
            inline void HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff =
                                StdRegions::NullVarCoeffMap,
                const MultiRegions::VarFactorsMap &varfactors =
                                 MultiRegions::NullVarFactorsMap,
                const Array<OneD, const NekDouble> &dirForcing =
                NullNekDouble1DArray,
                const bool PhysSpaceForcing = true);

            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionDiffusionReactionSolve(
                const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       const Array<OneD, const NekDouble>&
                       dirForcing = NullNekDouble1DArray);


            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionReactionSolve(
                const Array<OneD, Array<OneD, NekDouble> > &velocity,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                const NekDouble lambda,
                const Array<OneD, const NekDouble>&
                      dirForcing = NullNekDouble1DArray);

            ///
            MULTI_REGIONS_EXPORT void FwdTrans_BndConstrained(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);


            /// This function elementally evaluates the backward transformation
            /// of the global spectral/hp element expansion.
            inline void BwdTrans_IterPerExp (
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            ///
            inline void BwdTrans (
                const Array<OneD,
                const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);

            /// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(
                Array<OneD, NekDouble> &coord_0,
                Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

            // Homogeneous transforms
            inline void HomogeneousFwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                bool Shuff = true,
                bool UnShuff = true);

            inline void HomogeneousBwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                bool Shuff = true,
                bool UnShuff = true);

            inline void DealiasedProd(
                const Array<OneD, NekDouble> &inarray1,
                const Array<OneD, NekDouble> &inarray2,
                Array<OneD, NekDouble> &outarray);

            inline void DealiasedDotProd(
                const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                Array<OneD, Array<OneD, NekDouble> > &outarray);

            inline void GetBCValues(
                      Array<OneD, NekDouble> &BndVals,
                const Array<OneD, NekDouble> &TotField,
                int BndID);

            inline void NormVectorIProductWRTBase(
                Array<OneD, const NekDouble> &V1,
                Array<OneD, const NekDouble> &V2,
                Array<OneD, NekDouble> &outarray,
                int BndID);

            inline void NormVectorIProductWRTBase(
                Array<OneD, Array<OneD, NekDouble> > &V,
                Array<OneD, NekDouble> &outarray);

            /// Apply geometry information to each expansion.
            MULTI_REGIONS_EXPORT void ApplyGeomInfo();

            /// Reset geometry information and reset matrices
            MULTI_REGIONS_EXPORT void Reset()
            {
                v_Reset();
            }

            void WriteTecplotHeader(std::ostream &outfile,
                                    std::string var = "")
            {
                v_WriteTecplotHeader(outfile, var);
            }

            void WriteTecplotZone(
                std::ostream &outfile,
                int expansion = -1)
            {
                v_WriteTecplotZone(outfile, expansion);
            }

            void WriteTecplotField(std::ostream &outfile,
                                   int expansion = -1)
            {
                v_WriteTecplotField(outfile, expansion);
            }

            void WriteTecplotConnectivity(std::ostream &outfile,
                                          int expansion = -1)
            {
                v_WriteTecplotConnectivity(outfile, expansion);
            }

            MULTI_REGIONS_EXPORT void WriteVtkHeader(std::ostream &outfile);
            MULTI_REGIONS_EXPORT void WriteVtkFooter(std::ostream &outfile);

            void WriteVtkPieceHeader(std::ostream &outfile, int expansion,
                                     int istrip = 0)
            {
                v_WriteVtkPieceHeader(outfile, expansion, istrip);
            }

            MULTI_REGIONS_EXPORT void WriteVtkPieceFooter(
                std::ostream &outfile,
                int expansion);

            void WriteVtkPieceData  (
                std::ostream &outfile,
                int expansion,
                std::string var = "v")
            {
                v_WriteVtkPieceData(outfile, expansion, var);
            }

            /// This function returns the dimension of the coordinates of the
            /// element \a eid.
            // inline
            MULTI_REGIONS_EXPORT int GetCoordim(int eid);

            /// Set the \a i th coefficiient in \a m_coeffs to value \a val
            inline void SetCoeff(int i, NekDouble val);

            /// Set the \a i th coefficiient in  #m_coeffs to value \a val
            inline void SetCoeffs(int i, NekDouble val);

            /// Set the  #m_coeffs array to inarray
            inline void SetCoeffsArray(Array<OneD, NekDouble> &inarray);

            /// This function returns the dimension of the shape of the
            /// element \a eid.
            // inline
            MULTI_REGIONS_EXPORT int GetShapeDimension();

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline const Array<OneD, const NekDouble> &GetCoeffs() const;

            /// Impose Dirichlet Boundary Conditions onto Array
            inline void ImposeDirichletConditions(
                Array<OneD,NekDouble>& outarray);


            /// Fill Bnd Condition expansion from the values stored in expansion
            inline void FillBndCondFromField(void);

            /// Fill Bnd Condition expansion in nreg from the values stored in expansion
            inline void FillBndCondFromField(const int nreg);

            /// Gathers the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT inline void LocalToGlobal(bool useComm = true);

            MULTI_REGIONS_EXPORT inline void LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool useComm = true);

            /// Scatters from the global coefficients
            /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
            /// \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT inline void GlobalToLocal(void);

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
             * @param   inarray     An array of size \f$N_\mathrm{dof}\f$
             *                      containing the global degrees of freedom
             *                      \f$\boldsymbol{x}_g\f$.
             * @param   outarray    The resulting local degrees of freedom
             *                      \f$\boldsymbol{x}_l\f$ will be stored in this
             *                      array of size \f$N_\mathrm{eof}\f$.
             */
            MULTI_REGIONS_EXPORT inline void GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeff(int i);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeffs(int i);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD, const NekDouble>
                &GetPhys()  const;

            /// This function calculates the \f$L_\infty\f$ error of the global
            /// spectral/hp element approximation.
            MULTI_REGIONS_EXPORT NekDouble Linf (
                const Array<OneD, const NekDouble> &inarray,
                const Array<OneD, const NekDouble> &soln = NullNekDouble1DArray);

            /// This function calculates the \f$L_2\f$ error with
            /// respect to soln of the global
            /// spectral/hp element approximation.
            NekDouble L2(
                const Array<OneD, const NekDouble> &inarray,
                const Array<OneD, const NekDouble> &soln = NullNekDouble1DArray)
            {
                return v_L2(inarray, soln);
            }

            /// Calculates the \f$H^1\f$ error of the global spectral/hp
            /// element approximation.
            MULTI_REGIONS_EXPORT NekDouble H1 (
                const Array<OneD, const NekDouble> &inarray,
                const Array<OneD, const NekDouble> &soln = NullNekDouble1DArray);

            NekDouble Integral (const Array<OneD, const NekDouble> &inarray)
            {
                return v_Integral(inarray);
            }

            NekDouble VectorFlux(const Array<OneD, Array<OneD, NekDouble> > &inarray)
            {
                return v_VectorFlux(inarray);
            }

            /// This function calculates the energy associated with
            /// each one of the modesof a 3D homogeneous nD expansion
            Array<OneD, const NekDouble> HomogeneousEnergy (void)
            {
                return v_HomogeneousEnergy();
            }

            /// This function sets the Spectral Vanishing Viscosity
            /// in homogeneous1D expansion.
            void SetHomo1DSpecVanVisc(Array<OneD, NekDouble> visc)
            {
                v_SetHomo1DSpecVanVisc(visc);
            }

            /// This function returns a vector containing the wave
            /// numbers in z-direction associated
            /// with the 3D homogenous expansion. Required if a
            /// parellelisation is applied in the Fourier direction
            Array<OneD, const unsigned int> GetZIDs(void)
            {
                return v_GetZIDs();
            }

            /// This function returns the transposition class
            /// associaed with the homogeneous expansion.
            LibUtilities::TranspositionSharedPtr GetTransposition(void)
            {
                return v_GetTransposition();
            }

            /// This function returns the Width of homogeneous direction
            /// associaed with the homogeneous expansion.
            NekDouble GetHomoLen(void)
            {
                return v_GetHomoLen();
            }

            /// This function sets the Width of homogeneous direction
            /// associaed with the homogeneous expansion.
            void SetHomoLen(const NekDouble lhom)
            {
                return v_SetHomoLen(lhom);
            }

            /// This function returns a vector containing the wave
            /// numbers in y-direction associated
            /// with the 3D homogenous expansion. Required if a
            /// parellelisation is applied in the Fourier direction
            Array<OneD, const unsigned int> GetYIDs(void)
            {
                return v_GetYIDs();
            }

            /// This function interpolates the physical space points in
            /// \a inarray to \a outarray using the same points defined in the
            /// expansion but where the number of points are rescaled
            /// by \a 1DScale
            void PhysInterp1DScaled(
                const NekDouble scale,
                const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble>  &outarray)
            {
                v_PhysInterp1DScaled(scale, inarray,outarray);
            }

            /// This function Galerkin projects the physical space points in
            /// \a inarray to \a outarray where inarray is assumed to
            /// be defined in the expansion but where the number of
            /// points are rescaled by \a 1DScale
            void PhysGalerkinProjection1DScaled(
                const NekDouble scale,
                const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray)
            {
                v_PhysGalerkinProjection1DScaled(scale, inarray, outarray);
            }

            /// This function returns the number of elements in the expansion.
            inline int GetExpSize(void);


            /// This function returns the number of elements in the
            /// expansion which may be different for a homogeoenous extended
            /// expansionp.
            inline int GetNumElmts(void)
            {
                return v_GetNumElmts();
            }

            /// This function returns the vector of elements in the expansion.
            inline const std::shared_ptr<LocalRegions::ExpansionVector>
                    GetExp() const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion of the \f$n^{\mathrm{th}}\f$ element.
            inline LocalRegions::ExpansionSharedPtr& GetExp(int n) const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion containing the arbitrary point given by \a gloCoord.
            MULTI_REGIONS_EXPORT LocalRegions::ExpansionSharedPtr& GetExp(
                const Array<OneD, const NekDouble> &gloCoord);

            /** This function returns the index of the local elemental
             * expansion containing the arbitrary point given by \a gloCoord.
             **/
            MULTI_REGIONS_EXPORT int GetExpIndex(
                const Array<OneD, const NekDouble> &gloCoord,
                NekDouble tol = 0.0,
                bool returnNearestElmt = false);

            /** This function returns the index and the Local
             * Cartesian Coordinates \a locCoords of the local
             * elemental expansion containing the arbitrary point
             * given by \a gloCoords.
             **/
            MULTI_REGIONS_EXPORT int GetExpIndex(
                const Array<OneD, const NekDouble> &gloCoords,
                Array<OneD, NekDouble>       &locCoords,
                NekDouble tol = 0.0,
                bool returnNearestElmt = false);

            /** This function return the expansion field value
             * at the coordinates given as input.
             **/
            MULTI_REGIONS_EXPORT NekDouble PhysEvaluate(
                const Array<OneD, const NekDouble> &coords,
                const Array<OneD, const NekDouble> &phys);

            /// Get the start offset position for a global list of #m_coeffs
            /// correspoinding to element n.
            inline int GetCoeff_Offset(int n) const;

            /// Get the start offset position for a global list of m_phys
            /// correspoinding to element n.
            inline int GetPhys_Offset(int n) const;

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline Array<OneD, NekDouble> &UpdateCoeffs();

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            inline Array<OneD, NekDouble> &UpdatePhys();

            inline void PhysDeriv(
                Direction edir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &out_d);

            /// This function discretely evaluates the derivative of a function
            /// \f$f(\boldsymbol{x})\f$ on the domain consisting of all
            /// elements of the expansion.
            inline void PhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &out_d0,
                      Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                      Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            inline void PhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &out_d);

            inline void CurlCurl(
                Array<OneD, Array<OneD, NekDouble> > &Vel,
                Array<OneD, Array<OneD, NekDouble> > &Q);

            inline void PhysDirectionalDeriv(
                const Array<OneD, const NekDouble> &direction,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);

            inline void GetMovingFrames(
                const SpatialDomains::GeomMMF MMFdir,
                const Array<OneD, const NekDouble> &CircCentre,
                      Array<OneD, Array<OneD, NekDouble> > &outarray);

            // functions associated with DisContField
            inline const Array<OneD, const  std::shared_ptr<ExpList> >
                &GetBndCondExpansions();

            /// Get the weight value for boundary conditions
            inline const Array<OneD, const NekDouble>
                &GetBndCondBwdWeight();

            /// Set the weight value for boundary conditions
            inline void SetBndCondBwdWeight(
                const int index, 
                const NekDouble value);
      
            inline std::shared_ptr<ExpList> &UpdateBndCondExpansion(int i);

            inline void Upwind(
                const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
                const Array<OneD,                   const NekDouble>   &Fwd,
                const Array<OneD,                   const NekDouble>   &Bwd,
                      Array<OneD,                         NekDouble>   &Upwind);

            inline void Upwind(
                const Array<OneD, const NekDouble> &Vn,
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &Upwind);

            /**
             * Return a reference to the trace space associated with this
             * expansion list.
             */
            inline std::shared_ptr<ExpList> &GetTrace();

            inline std::shared_ptr<AssemblyMapDG> &GetTraceMap(void);

            inline const Array<OneD, const int> &GetTraceBndMap(void);

            inline void GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals);

            /// Get the length of elements in boundary normal direction
            inline void GetElmtNormalLength(
                Array<OneD, NekDouble>  &lengthsFwd,
                Array<OneD, NekDouble>  &lengthsBwd);

            /// Get the weight value for boundary conditions
            /// for boundary average and jump calculations
            MULTI_REGIONS_EXPORT void GetBwdWeight(
                Array<OneD, NekDouble>  &weightAver,
                Array<OneD, NekDouble>  &weightJump);

            inline void AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fx,
                const Array<OneD, const NekDouble> &Fy,
                Array<OneD, NekDouble> &outarray);

            inline void AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD, NekDouble> &outarray);

            inline void AddFwdBwdTraceIntegral(
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD, NekDouble> &outarray);

            inline void GetFwdBwdTracePhys(
                Array<OneD,NekDouble> &Fwd,
                Array<OneD,NekDouble> &Bwd);

            inline void GetFwdBwdTracePhys(
                const Array<OneD,const NekDouble> &field,
                      Array<OneD,NekDouble> &Fwd,
                      Array<OneD,NekDouble> &Bwd);

            /// GetFwdBwdTracePhys without parallel communication
            inline void GetFwdBwdTracePhysSerial(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);

            /// GetFwdBwdTracePhys of derivatives
            inline void GetFwdBwdTracePhysDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);
            
            /// GetFwdBwdTracePhysDeriv without parallel communication
            inline void GetFwdBwdTracePhysDerivSerial(
                const int                          Dir,
                const Array<OneD, const NekDouble> &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);
            
            /// GetFwdBwdTracePhys without filling boundary conditions
            inline void GetFwdBwdTracePhysNoBndFill(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);

            /// Add Fwd and Bwd value to field, 
            /// a reverse procedure of GetFwdBwdTracePhys
            inline void AddTraceQuadPhysToField(
                const Array<OneD, const NekDouble>  &Fwd,
                const Array<OneD, const NekDouble>  &Bwd,
                Array<OneD,       NekDouble>        &field);

            /// Fill Bwd with boundary conditions
            inline void FillBwdWithBound(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            
            /// Fill Bwd with boundary conditions for derivatives 
            inline void FillBwdWithBoundDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);

            /// Fill Bwd with boundary conditions
            inline void FillBwdWithBwdWeight(
                    Array<OneD,       NekDouble> &weightave,
                    Array<OneD,       NekDouble> &weightjmp);

            /// Copy and fill the Periodic boundaries
            inline void PeriodicBwdCopy(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);

            inline const std::vector<bool> &GetLeftAdjacentFaces(void) const;

            inline void ExtractTracePhys(Array<OneD,NekDouble> &outarray);

            inline void ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            inline const Array<OneD, const SpatialDomains::
                BoundaryConditionShPtr>& GetBndConditions();

            inline Array<OneD, SpatialDomains::
                BoundaryConditionShPtr>& UpdateBndConditions();

            inline void EvaluateBoundaryConditions(
                const NekDouble   time      = 0.0,
                const std::string varName   = "",
                const             NekDouble = NekConstants::kNekUnsetDouble,
                const             NekDouble = NekConstants::kNekUnsetDouble);

            // Routines for continous matrix solution
            /// This function calculates the result of the multiplication of a
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            inline void GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);

            MULTI_REGIONS_EXPORT void GeneralMatrixOp_IterPerExp(
                const GlobalMatrixKey      &gkey,
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,      NekDouble> &outarray);

            inline void SetUpPhysNormals();

            inline void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                             Array<OneD,int> &EdgeID);

            inline void GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays = true);

            inline void ExtractElmtToBndPhys(int i,
                                             const Array<OneD, NekDouble> &elmt,
                                             Array<OneD, NekDouble> &boundary);

            inline void ExtractPhysToBndElmt(int i,
                            const Array<OneD, const NekDouble> &phys,
                            Array<OneD, NekDouble> &bndElmt);

            inline void ExtractPhysToBnd(int i,
                            const Array<OneD, const NekDouble> &phys,
                            Array<OneD, NekDouble> &bnd);

            inline void GetBoundaryNormals(int i,
                            Array<OneD, Array<OneD, NekDouble> > &normals);

            MULTI_REGIONS_EXPORT void  GeneralGetFieldDefinitions(
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef,
                int NumHomoDir = 0,
                Array<OneD, LibUtilities::BasisSharedPtr> &HomoBasis =
                    LibUtilities::NullBasisSharedPtr1DArray,
                std::vector<NekDouble> &HomoLen =
                    LibUtilities::NullNekDoubleVector,
                bool  homoStrips = false,
                std::vector<unsigned int> &HomoSIDs =
                    LibUtilities::NullUnsignedIntVector,
                std::vector<unsigned int> &HomoZIDs =
                    LibUtilities::NullUnsignedIntVector,
                std::vector<unsigned int> &HomoYIDs =
                    LibUtilities::NullUnsignedIntVector);


            std::map<int, RobinBCInfoSharedPtr> GetRobinBCInfo()
            {
                return v_GetRobinBCInfo();
            }

            void GetPeriodicEntities(
                PeriodicMap &periodicVerts,
                PeriodicMap &periodicEdges,
                PeriodicMap &periodicFaces = NullPeriodicMap)
            {
                v_GetPeriodicEntities(periodicVerts, periodicEdges, periodicFaces);
            }

            std::vector<LibUtilities::FieldDefinitionsSharedPtr>
                GetFieldDefinitions()
            {
                return v_GetFieldDefinitions();
            }


            void GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
            {
                v_GetFieldDefinitions(fielddef);
            }



            /// Append the element data listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata)
            {
                v_AppendFieldData(fielddef,fielddata);
            }


            /// Append the data in coeffs listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata,
                Array<OneD, NekDouble> &coeffs)
            {
                v_AppendFieldData(fielddef,fielddata,coeffs);
            }


            /** \brief Extract the data in fielddata into the coeffs
             * using the basic ExpList Elemental expansions rather
             * than planes in homogeneous case
             */
            MULTI_REGIONS_EXPORT void ExtractElmtDataToCoeffs(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata,
                std::string &field,
                Array<OneD, NekDouble> &coeffs);


            /** \brief Extract the data from fromField using
             * fromExpList the coeffs using the basic ExpList
             * Elemental expansions rather than planes in homogeneous
             * case
             */
            MULTI_REGIONS_EXPORT  void ExtractCoeffsToCoeffs(
                const std::shared_ptr<ExpList> &fromExpList,
                const Array<OneD, const NekDouble> &fromCoeffs,
                      Array<OneD, NekDouble> &toCoeffs);


            //Extract data in fielddata into the m_coeffs_list for the 3D stability analysis (base flow is 2D)
            MULTI_REGIONS_EXPORT void ExtractDataToCoeffs(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata,
                std::string &field,
                Array<OneD, NekDouble> &coeffs);

            MULTI_REGIONS_EXPORT void GenerateElementVector(
                const int ElementID,
                const NekDouble scalar1,
                const NekDouble scalar2,
                Array<OneD, NekDouble> &outarray);

            /// Returns a shared pointer to the current object.
            std::shared_ptr<ExpList> GetSharedThisPtr()
            {
                return shared_from_this();
            }

            /// Returns the session object
            std::shared_ptr<LibUtilities::SessionReader> GetSession() const
            {
                return m_session;
            }

            /// Returns the comm object
            std::shared_ptr<LibUtilities::Comm> GetComm()
            {
                return m_comm;
            }

            SpatialDomains::MeshGraphSharedPtr GetGraph()
            {
                return m_graph;
            }

            // Wrapper functions for Homogeneous Expansions
            inline LibUtilities::BasisSharedPtr  GetHomogeneousBasis(void)
            {
                return v_GetHomogeneousBasis();
            }

            std::shared_ptr<ExpList> &GetPlane(int n)
            {
                return v_GetPlane(n);
            }
             
            //expansion type
            ExpansionType m_expType;

            MULTI_REGIONS_EXPORT void CreateCollections(
                    Collections::ImplementationType ImpType
                                                    = Collections::eNoImpType);

            MULTI_REGIONS_EXPORT void ClearGlobalLinSysManager(void);

            /// Get m_coeffs to elemental value map
            MULTI_REGIONS_EXPORT inline const 
                Array<OneD, const std::pair<int, int> > 
                &GetCoeffsToElmt() const;

            MULTI_REGIONS_EXPORT inline const LocTraceToTraceMapSharedPtr 
                &GetLocTraceToTraceMap() const;
        protected:
            /// Definition of the total number of degrees of freedom and
            /// quadrature points and offsets to access data
            void SetCoeffPhysOffsets();

            std::shared_ptr<DNekMat> GenGlobalMatrixFull(
                const GlobalLinSysKey &mkey,
                const std::shared_ptr<AssemblyMapCG> &locToGloMap);

            /// Communicator
            LibUtilities::CommSharedPtr m_comm;

            /// Session
            LibUtilities::SessionReaderSharedPtr m_session;

            /// Mesh associated with this expansion list.
            SpatialDomains::MeshGraphSharedPtr m_graph;

            /// The total number of local degrees of freedom. #m_ncoeffs
            /// \f$=N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_l\f$
            int m_ncoeffs;

            /** The total number of quadrature points. #m_npoints
             *\f$=Q_{\mathrm{tot}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_Q\f$
             **/
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
            std::shared_ptr<LocalRegions::ExpansionVector> m_exp;

            Collections::CollectionVector m_collections;

            /// Offset of elemental data into the array #m_coeffs
            std::vector<int>  m_coll_coeff_offset;

            /// Offset of elemental data into the array #m_phys
            std::vector<int>  m_coll_phys_offset;

            /// Offset of elemental data into the array #m_coeffs
            Array<OneD, int>  m_coeff_offset;

            /// Offset of elemental data into the array #m_phys
            Array<OneD, int>  m_phys_offset;

            /// m_coeffs to elemental value map
            Array<OneD, std::pair<int, int> >  m_coeffsToElmt;

            NekOptimize::GlobalOptParamSharedPtr m_globalOptParam;

            BlockMatrixMapShPtr  m_blockMat;

            //@todo should this be in ExpList or ExpListHomogeneous1D.cpp
            // it's a bool which determine if the expansion is in the wave space (coefficient space)
            // or not
            bool m_WaveSpace;

            /// Mapping from geometry ID of element to index inside #m_exp
            std::unordered_map<int, int> m_elmtToExpId;

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

            /// Generates a global matrix from the given key and map.
            std::shared_ptr<GlobalMatrix>  GenGlobalMatrix(
                const GlobalMatrixKey &mkey,
                const std::shared_ptr<AssemblyMapCG> &locToGloMap);


            void GlobalEigenSystem(
                const std::shared_ptr<DNekMat> &Gmat,
                Array<OneD, NekDouble> &EigValsReal,
                Array<OneD, NekDouble> &EigValsImag,
                Array<OneD, NekDouble> &EigVecs
                = NullNekDouble1DArray);


            /// This operation constructs the global linear system of type \a
            /// mkey.
            std::shared_ptr<GlobalLinSys>  GenGlobalLinSys(
                const GlobalLinSysKey &mkey,
                const std::shared_ptr<AssemblyMapCG> &locToGloMap);

            /// Generate a GlobalLinSys from information provided by the key
            /// "mkey" and the mapping provided in LocToGloBaseMap.
            std::shared_ptr<GlobalLinSys> GenGlobalBndLinSys(
                const GlobalLinSysKey     &mkey,
                const AssemblyMapSharedPtr &locToGloMap);

            // Virtual prototypes

            virtual int v_GetNumElmts(void)
            {
                return (*m_exp).size();
            }

            virtual const Array<OneD,const std::shared_ptr<ExpList> >
                &v_GetBndCondExpansions(void);

            virtual const Array<OneD, const NekDouble>
                &v_GetBndCondBwdWeight();

            virtual void v_SetBndCondBwdWeight(
                const int index, 
                const NekDouble value);

            virtual std::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i);

            virtual void v_Upwind(
                const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
                const Array<OneD,                   const NekDouble>   &Fwd,
                const Array<OneD,                   const NekDouble>   &Bwd,
                      Array<OneD,                         NekDouble>   &Upwind);

            virtual void v_Upwind(
                const Array<OneD, const NekDouble> &Vn,
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &Upwind);

            virtual std::shared_ptr<ExpList> &v_GetTrace();

            virtual std::shared_ptr<AssemblyMapDG> &v_GetTraceMap();

            virtual const Array<OneD, const int> &v_GetTraceBndMap();

            virtual void v_GetNormals(
                Array<OneD, Array<OneD, NekDouble> > &normals);

            virtual void v_GetElmtNormalLength(
                Array<OneD, NekDouble>  &lengthsFwd,
                Array<OneD, NekDouble>  &lengthsBwd);

            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fx,
                const Array<OneD, const NekDouble> &Fy,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_AddFwdBwdTraceIntegral(
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

            virtual void v_GetFwdBwdTracePhysDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble>  &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);
            
            virtual void v_GetFwdBwdTracePhysDerivSerial(
                const int                          Dir,
                const Array<OneD, const NekDouble>  &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);

            virtual void v_GetFwdBwdTracePhysNoBndFill(
                const Array<OneD, const NekDouble>  &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);

            virtual void v_GetFwdBwdTracePhysSerial(
                const Array<OneD, const NekDouble>  &field,
                      Array<OneD, NekDouble> &Fwd,
                      Array<OneD, NekDouble> &Bwd);
            virtual void v_GetFwdBwdTracePhysInterior(
                const Array<OneD,const NekDouble>  &field,
                      Array<OneD,NekDouble> &Fwd,
                      Array<OneD,NekDouble> &Bwd);
            
            virtual void v_AddTraceQuadPhysToField(
                const Array<OneD, const NekDouble>  &Fwd,
                const Array<OneD, const NekDouble>  &Bwd,
                Array<OneD,       NekDouble>        &field);
                      
            virtual void v_FillBwdWithBound(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            
            virtual void v_FillBwdWithBoundDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            
            virtual void v_FillBwdWithBwdWeight(
                Array<OneD,       NekDouble> &weightave,
                Array<OneD,       NekDouble> &weightjmp);

            virtual void v_PeriodicBwdCopy(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);

            virtual const std::vector<bool> &v_GetLeftAdjacentFaces(void) const;

            virtual void v_ExtractTracePhys(
                Array<OneD,NekDouble> &outarray);

            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            virtual void v_MultiplyByInvMassMatrix(
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);

            virtual void v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing);

            virtual void v_LinearAdvectionDiffusionReactionSolve(
                const Array<OneD, Array<OneD, NekDouble> > &velocity,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                const NekDouble lambda,
                const Array<OneD, const NekDouble>&
                      dirForcing = NullNekDouble1DArray);

            virtual void v_LinearAdvectionReactionSolve(
                const Array<OneD, Array<OneD, NekDouble> > &velocity,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                const NekDouble lambda,
                const Array<OneD, const NekDouble>&
                      dirForcing = NullNekDouble1DArray);

            // wrapper functions about virtual functions
            virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            virtual void v_FillBndCondFromField();

            virtual void v_FillBndCondFromField(const int nreg);

            virtual void v_Reset();

            virtual void v_LocalToGlobal(bool UseComm);

            virtual void v_LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool UseComm);

            virtual void v_GlobalToLocal(void);

            virtual void v_GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);

            virtual void v_BwdTrans(
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);
            
            virtual void v_BwdTrans_IterPerExp(
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            virtual void v_FwdTrans(
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);

            virtual void v_FwdTrans_IterPerExp(
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            virtual void v_FwdTrans_BndConstrained(
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray);

            virtual void v_SmoothField(Array<OneD,NekDouble> &field);

            virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,       NekDouble> &outarray);
            
            virtual void v_IProductWRTBase_IterPerExp(
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,      NekDouble> &outarray);

            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);
            
            virtual void v_GetCoords(
                Array<OneD, NekDouble> &coord_0,
                Array<OneD, NekDouble> &coord_1,
                Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

            virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2);

            virtual void v_PhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &out_d);

            virtual void v_PhysDeriv(
                Direction edir,
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &out_d);

            virtual void v_CurlCurl(
                Array<OneD, Array<OneD, NekDouble> > &Vel,
                Array<OneD, Array<OneD, NekDouble> > &Q);

            virtual void v_PhysDirectionalDeriv(
                const Array<OneD, const NekDouble> &direction,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_GetMovingFrames(
                const SpatialDomains::GeomMMF MMFdir,
                const Array<OneD, const NekDouble> &CircCentre,
                      Array<OneD, Array<OneD, NekDouble> > &outarray);

            virtual void v_HomogeneousFwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                bool Shuff = true,
                bool UnShuff = true);

            virtual void v_HomogeneousBwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                bool Shuff = true,
                bool UnShuff = true);

            virtual void v_DealiasedProd(
                const Array<OneD, NekDouble> &inarray1,
                const Array<OneD, NekDouble> &inarray2,
                Array<OneD, NekDouble> &outarray);

            virtual void v_DealiasedDotProd(
                const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                Array<OneD, Array<OneD, NekDouble> > &outarray);

            virtual void v_GetBCValues(
                      Array<OneD, NekDouble> &BndVals,
                const Array<OneD, NekDouble> &TotField,
                int BndID);

            virtual void v_NormVectorIProductWRTBase(
                Array<OneD, const NekDouble> &V1,
                Array<OneD, const NekDouble> &V2,
                Array<OneD, NekDouble> &outarray,
                int BndID);

            virtual void v_NormVectorIProductWRTBase(
                Array<OneD, Array<OneD, NekDouble> > &V,
                Array<OneD, NekDouble> &outarray);

            virtual void v_SetUpPhysNormals();

            virtual void v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                                Array<OneD,int> &EdgeID);

            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);

            virtual void v_ExtractElmtToBndPhys( const int                      i,
                                                 const Array<OneD, NekDouble> & elmt,
                                                       Array<OneD, NekDouble> & boundary );

            virtual void v_ExtractPhysToBndElmt( const int                            i,
                                                 const Array<OneD, const NekDouble> & phys,
                                                       Array<OneD, NekDouble>       & bndElmt );

            virtual void v_ExtractPhysToBnd( const int                            i,
                                             const Array<OneD, const NekDouble> & phys,
                                                   Array<OneD, NekDouble>       & bnd );

            virtual void v_GetBoundaryNormals(int i,
                            Array<OneD, Array<OneD, NekDouble> > &normals);

            virtual std::vector<LibUtilities::FieldDefinitionsSharedPtr>
                v_GetFieldDefinitions(void);

            virtual void  v_GetFieldDefinitions(
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef);

            virtual void v_AppendFieldData(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata);

            virtual void v_AppendFieldData(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata,
                Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractDataToCoeffs(
                LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                std::vector<NekDouble> &fielddata, std::string &field,
                Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractCoeffsToCoeffs(const std::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs);

            virtual void v_WriteTecplotHeader(std::ostream &outfile,
                                              std::string var = "");
            virtual void v_WriteTecplotZone(std::ostream &outfile,
                                            int expansion);
            virtual void v_WriteTecplotField(std::ostream &outfile,
                                             int expansion);
            virtual void v_WriteTecplotConnectivity(std::ostream &outfile,
                                                    int expansion);
            virtual void v_WriteVtkPieceHeader(
                std::ostream &outfile,
                int expansion,
                int istrip);

            virtual void v_WriteVtkPieceData(
                std::ostream &outfile,
                int expansion,
                std::string var);

            virtual NekDouble v_L2(
                const Array<OneD, const NekDouble> &phys,
                const Array<OneD, const NekDouble> &soln = NullNekDouble1DArray);

            virtual NekDouble v_Integral (
                const Array<OneD, const NekDouble> &inarray);
            virtual NekDouble v_VectorFlux (
                const Array<OneD, Array<OneD, NekDouble> > &inarray);

            virtual Array<OneD, const NekDouble> v_HomogeneousEnergy(void);
            virtual LibUtilities::TranspositionSharedPtr v_GetTransposition(void);
            virtual NekDouble v_GetHomoLen(void);
            virtual void      v_SetHomoLen(const NekDouble lhom);
            virtual Array<OneD, const unsigned int> v_GetZIDs(void);
            virtual Array<OneD, const unsigned int> v_GetYIDs(void);

            // 1D Scaling and projection
            virtual void v_PhysInterp1DScaled(
                const NekDouble scale, const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_PhysGalerkinProjection1DScaled(
                const NekDouble scale,
                const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_ClearGlobalLinSysManager(void);

            void ExtractFileBCs(const std::string                &fileName,
                                LibUtilities::CommSharedPtr       comm,
                                const std::string                &varName,
                                const std::shared_ptr<ExpList>  locExpList);

            // Utility function for a common case of retrieving a
            // BoundaryCondition from a boundary condition collection.
            MULTI_REGIONS_EXPORT
                static SpatialDomains::BoundaryConditionShPtr
                    GetBoundaryCondition(const SpatialDomains::
                            BoundaryConditionCollection& collection,
                            unsigned int index, const std::string& variable);


        private:

            virtual const LocTraceToTraceMapSharedPtr 
                &v_GetLocTraceToTraceMap() const;

            virtual const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &v_GetBndConditions();

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                &v_UpdateBndConditions();

            virtual void v_EvaluateBoundaryConditions(
                const NekDouble   time    = 0.0,
                const std::string varName = "",
                const NekDouble   x2_in   = NekConstants::kNekUnsetDouble,
                const NekDouble   x3_in   = NekConstants::kNekUnsetDouble);

            virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo(void);


            virtual void v_GetPeriodicEntities(
                PeriodicMap &periodicVerts,
                PeriodicMap &periodicEdges,
                PeriodicMap &periodicFaces);

            // Homogeneous direction wrapper functions.
            virtual LibUtilities::BasisSharedPtr  v_GetHomogeneousBasis(void)
            {
                NEKERROR(ErrorUtil::efatal,
                    "This method is not defined or valid for this class type");
                return LibUtilities::NullBasisSharedPtr;
            }

            // wrapper function to set viscosity for Homo1D expansion
            virtual void v_SetHomo1DSpecVanVisc(Array<OneD, NekDouble> visc)
            {
                boost::ignore_unused(visc);
                NEKERROR(ErrorUtil::efatal,
                    "This method is not defined or valid for this class type");
            }


            virtual std::shared_ptr<ExpList> &v_GetPlane(int n);
        };


        /// Shared pointer to an ExpList object.
        typedef std::shared_ptr<ExpList>      ExpListSharedPtr;
        /// An empty ExpList object.
        static ExpList NullExpList;
        static ExpListSharedPtr NullExpListSharedPtr;

        // Inline routines follow.

        /**
         * Returns the total number of local degrees of freedom
         * \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
         */
        inline int ExpList::GetNcoeffs() const
        {
            return m_ncoeffs;
        }

        inline int ExpList::GetNcoeffs(const int eid) const
        {
            return (*m_exp)[eid]->GetNcoeffs();
        }

        /**
         * Evaulates the maximum number of modes in the elemental basis
         * order over all elements
         */
        inline int ExpList::EvalBasisNumModesMax() const
        {
            unsigned int i;
            int returnval = 0;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval = (std::max)(returnval,
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
            unsigned int i;
            Array<OneD,int> returnval((*m_exp).size(),0);

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval[i]
                    = (std::max)(returnval[i],(*m_exp)[i]->EvalBasisNumModesMax());
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


        inline int ExpList::Get1DScaledTotPoints(const NekDouble scale) const
        {
            int returnval = 0;
            int cnt;
            int nbase = (*m_exp)[0]->GetNumBases();

            for(int i = 0; i < (*m_exp).size(); ++i)
            {
                cnt = 1;
                for(int j = 0; j < nbase; ++j)
                {
                    cnt *= (int)(scale*((*m_exp)[i]->GetNumPoints(j)));
                }
                returnval += cnt;
            }
            return returnval;
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
        inline void ExpList::SetWaveSpace(const bool wavespace)
        {
            m_WaveSpace = wavespace;
        }

        /**
         *
         */
        inline bool ExpList::GetWaveSpace() const
        {
            return m_WaveSpace;
        }

        /// Set the \a i th value of\a m_phys to value \a val
        inline void ExpList::SetPhys(int i, NekDouble val)
        {
            m_phys[i] = val;
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
            ASSERTL0(inarray.size() == m_npoints,
                     "Input array does not have correct number of elements.");

            Vmath::Vcopy(m_npoints,&inarray[0],1,&m_phys[0],1);
            m_physState = true;
        }


        inline void ExpList::SetPhysArray(Array<OneD, NekDouble> &inarray)
        {
            m_phys = inarray;
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
            Array<OneD, NekDouble> &outarray)
        {
            v_IProductWRTBase(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::IProductWRTBase_IterPerExp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase_IterPerExp(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::FwdTrans(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,       NekDouble> &outarray)
        {
            v_FwdTrans(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::FwdTrans_IterPerExp (
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray)
        {
            v_FwdTrans_IterPerExp(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::FwdTrans_BndConstrained (
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray)
        {
            v_FwdTrans_BndConstrained(inarray,outarray);
        }


        /**
         *
         */
        inline void ExpList::SmoothField(Array<OneD,NekDouble> &field)
        {
            v_SmoothField(field);
        }

        /**
         *
         */
        inline void ExpList::BwdTrans (
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,       NekDouble> &outarray)
        {
            v_BwdTrans(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::BwdTrans_IterPerExp (
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_BwdTrans_IterPerExp(inarray,outarray);
        }


        /**
         *
         */
        inline void ExpList::MultiplyByInvMassMatrix(
            const Array<OneD,const NekDouble> &inarray,
            Array<OneD,      NekDouble> &outarray)
        {
            v_MultiplyByInvMassMatrix(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::HelmSolve(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::ConstFactorMap &factors,
            const StdRegions::VarCoeffMap &varcoeff,
            const MultiRegions::VarFactorsMap &varfactors,
            const Array<OneD, const NekDouble> &dirForcing,
            const bool PhysSpaceForcing)

        {
            v_HelmSolve(inarray, outarray, factors, varcoeff,
                        varfactors, dirForcing, PhysSpaceForcing);
        }


        /**
         *
         */
        inline void ExpList::LinearAdvectionDiffusionReactionSolve(
            const Array<OneD, Array<OneD, NekDouble> > &velocity,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray,
            const NekDouble lambda,
            const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionDiffusionReactionSolve(velocity,inarray, outarray,
                                                lambda, dirForcing);
        }

        inline void ExpList::LinearAdvectionReactionSolve(
            const Array<OneD, Array<OneD, NekDouble> > &velocity,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray,
            const NekDouble lambda,
            const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionReactionSolve(velocity,inarray, outarray,
                                           lambda, dirForcing);
        }

        /**
         *
         */
        inline void ExpList::GetCoords(Array<OneD, NekDouble> &coord_0,
                                       Array<OneD, NekDouble> &coord_1,
                                       Array<OneD, NekDouble> &coord_2)

        {
            v_GetCoords(coord_0,coord_1,coord_2);
        }


        /**
         *
         */
        inline void ExpList::GetMovingFrames(
            const SpatialDomains::GeomMMF MMFdir,
            const Array<OneD, const NekDouble> &CircCentre,
                  Array<OneD, Array<OneD, NekDouble> > &outarray)
        {
             v_GetMovingFrames(MMFdir,CircCentre,outarray);
        }


        /**
         *
         */
        inline void ExpList::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD, NekDouble> &out_d0,
                                       Array<OneD, NekDouble> &out_d1,
                                       Array<OneD, NekDouble> &out_d2)
        {
            v_PhysDeriv(inarray,out_d0,out_d1,out_d2);
        }

        /**
         *
         */
        inline void ExpList::PhysDeriv(
            const int dir,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(dir,inarray,out_d);
        }

        inline void ExpList::PhysDeriv(
            Direction edir,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(edir, inarray,out_d);
        }


        /**
         *
         */
        inline void ExpList::PhysDirectionalDeriv(
            const Array<OneD, const NekDouble> &direction,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            v_PhysDirectionalDeriv(direction, inarray, outarray);
        }


        /**
         *
         */

        inline void ExpList::CurlCurl(
                Array<OneD, Array<OneD, NekDouble> > &Vel,
                Array<OneD, Array<OneD, NekDouble> > &Q)
        {
            v_CurlCurl(Vel, Q);
        }

        /**
         *
         */
        inline void ExpList::HomogeneousFwdTrans(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray,
            bool Shuff,
            bool UnShuff)
        {
            v_HomogeneousFwdTrans(inarray,outarray,Shuff,UnShuff);
        }

        /**
         *
         */
        inline void ExpList::HomogeneousBwdTrans(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray,
            bool Shuff,
            bool UnShuff)
        {
            v_HomogeneousBwdTrans(inarray,outarray,Shuff,UnShuff);
        }

        /**
         *
         */
        inline void ExpList::DealiasedProd(
            const Array<OneD, NekDouble> &inarray1,
            const Array<OneD, NekDouble> &inarray2,
            Array<OneD, NekDouble> &outarray)
        {
            v_DealiasedProd(inarray1,inarray2,outarray);
        }

        /**
         *
         */
        inline void ExpList::DealiasedDotProd(
                const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                Array<OneD, Array<OneD, NekDouble> > &outarray)
        {
            v_DealiasedDotProd(inarray1,inarray2,outarray);
        }

        /**
         *
         */
        inline void ExpList::GetBCValues(
                  Array<OneD, NekDouble> &BndVals,
            const Array<OneD, NekDouble> &TotField,
            int BndID)
        {
            v_GetBCValues(BndVals,TotField,BndID);
        }

        /**
         *
         */
        inline void ExpList::NormVectorIProductWRTBase(
            Array<OneD, const NekDouble> &V1,
            Array<OneD, const NekDouble> &V2,
            Array<OneD, NekDouble> &outarray,
            int BndID)
        {
            v_NormVectorIProductWRTBase(V1,V2,outarray,BndID);
        }

        inline void ExpList::NormVectorIProductWRTBase(
            Array<OneD, Array<OneD, NekDouble> > &V,
            Array<OneD, NekDouble> &outarray)
        {
            v_NormVectorIProductWRTBase(V, outarray);
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
         * @param   eid         The index of the element to be checked.
         * @return  The dimension of the shape of the specific element.
         */
        inline int ExpList::GetShapeDimension()
        {
            return (*m_exp)[0]->GetShapeDimension();
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


        inline void ExpList::SetCoeffsArray(Array<OneD, NekDouble> &inarray)
        {
            m_coeffs = inarray;
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

        inline void ExpList::ImposeDirichletConditions(
            Array<OneD,NekDouble>& outarray)
        {
            v_ImposeDirichletConditions(outarray);
        }

        inline void ExpList::FillBndCondFromField(void)
        {
            v_FillBndCondFromField();
        }

        inline void ExpList::FillBndCondFromField(const int nreg)
        {
            v_FillBndCondFromField(nreg);
        }

        inline void ExpList::LocalToGlobal(bool useComm)
        {
            v_LocalToGlobal(useComm);
        }

        inline void ExpList::LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool useComm)
        {
            v_LocalToGlobal(inarray, outarray,useComm);
        }

        inline void ExpList::GlobalToLocal(void)
        {
            v_GlobalToLocal();
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
         * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where
         * \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @param   inarray     An array of size \f$N_\mathrm{dof}\f$
         *                      containing the global degrees of freedom
         *                      \f$\boldsymbol{x}_g\f$.
         * @param   outarray    The resulting local degrees of freedom
         *                      \f$\boldsymbol{x}_l\f$ will be stored in this
         *                      array of size \f$N_\mathrm{eof}\f$.
         */
        inline void ExpList::GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
        {
            v_GlobalToLocal(inarray, outarray);
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
        inline LocalRegions::ExpansionSharedPtr& ExpList::GetExp(int n) const
        {
            return (*m_exp)[n];
        }

        /**
         * @return  (A const shared pointer to) the local expansion vector #m_exp
         */
        inline const std::shared_ptr<LocalRegions::ExpansionVector>
            ExpList::GetExp(void) const
        {
            return m_exp;
        }


        /**
         *
         */
        inline int ExpList::GetCoeff_Offset(int n) const
        {
            return m_coeff_offset[n];
        }

        /**
         *
         */
        inline int ExpList::GetPhys_Offset(int n) const
        {
            return m_phys_offset[n];
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them, rather use the function #GetCoeffs instead.
         *
         * @return  (A reference to) the array #m_coeffs.
         */
        inline Array<OneD, NekDouble> &ExpList::UpdateCoeffs()
        {
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
        inline const Array<OneD, const  std::shared_ptr<ExpList> >
            &ExpList::GetBndCondExpansions()
        {
            return v_GetBndCondExpansions();
        }

        /// Get m_coeffs to elemental value map
        MULTI_REGIONS_EXPORT inline const 
            Array<OneD, const std::pair<int, int> > 
            &ExpList::GetCoeffsToElmt() const
        {
            return m_coeffsToElmt;
        }

        MULTI_REGIONS_EXPORT inline const LocTraceToTraceMapSharedPtr 
            &ExpList::GetLocTraceToTraceMap() const
        {
            return v_GetLocTraceToTraceMap();
        }

        inline const Array<OneD, const  NekDouble >
            &ExpList::GetBndCondBwdWeight()
        {
            return v_GetBndCondBwdWeight();
        }

        inline void ExpList::SetBndCondBwdWeight(
            const int index, 
            const NekDouble value)
        {
            v_SetBndCondBwdWeight(index, value);
        }
        
        inline std::shared_ptr<ExpList>  &ExpList::UpdateBndCondExpansion(int i)
        {
            return v_UpdateBndCondExpansion(i);
        }

        inline void ExpList::Upwind(
            const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
            const Array<OneD,                   const NekDouble>   &Fwd,
            const Array<OneD,                   const NekDouble>   &Bwd,
                  Array<OneD,                         NekDouble>   &Upwind)
        {
            v_Upwind(Vec, Fwd, Bwd, Upwind);
        }

        inline void ExpList::Upwind(
            const Array<OneD, const NekDouble> &Vn,
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &Upwind)
        {
            v_Upwind(Vn, Fwd, Bwd, Upwind);
        }

        inline std::shared_ptr<ExpList> &ExpList::GetTrace()
        {
            return v_GetTrace();
        }

        inline std::shared_ptr<AssemblyMapDG> &ExpList::GetTraceMap()
        {
            return v_GetTraceMap();
        }

        inline const Array<OneD, const int> &ExpList::GetTraceBndMap()
        {
            return v_GetTraceBndMap();
        }

        inline void ExpList::GetNormals(
            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            v_GetNormals(normals);
        }

        inline void ExpList::GetElmtNormalLength(
            Array<OneD, NekDouble>  &lengthsFwd,
            Array<OneD, NekDouble>  &lengthsBwd)
        {
            v_GetElmtNormalLength(lengthsFwd, lengthsBwd);
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

        inline void ExpList::AddFwdBwdTraceIntegral(
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD, NekDouble> &outarray)
        {
            v_AddFwdBwdTraceIntegral(Fwd,Bwd,outarray);
        }

        inline void ExpList::GetFwdBwdTracePhys(
            Array<OneD,NekDouble> &Fwd,
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

        inline void ExpList::GetFwdBwdTracePhysSerial(
            const Array<OneD, const NekDouble>  &field,
                  Array<OneD, NekDouble> &Fwd,
                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysSerial(field, Fwd, Bwd);
        }

        inline void ExpList::GetFwdBwdTracePhysNoBndFill(
            const Array<OneD, const NekDouble>  &field,
                  Array<OneD, NekDouble> &Fwd,
                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysNoBndFill(field, Fwd, Bwd);
        }

        inline void ExpList::GetFwdBwdTracePhysDeriv(
            const int                          Dir,
            const Array<OneD, const NekDouble>  &field,
                  Array<OneD, NekDouble> &Fwd,
                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysDeriv(Dir, field, Fwd, Bwd);
        }

        inline void ExpList::GetFwdBwdTracePhysDerivSerial(
            const int                          Dir,
            const Array<OneD, const NekDouble>  &field,
                  Array<OneD, NekDouble> &Fwd,
                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysDerivSerial(Dir, field, Fwd, Bwd);
        }

        inline void ExpList::AddTraceQuadPhysToField(
                const Array<OneD, const NekDouble>  &Fwd,
                const Array<OneD, const NekDouble>  &Bwd,
                Array<OneD,       NekDouble>        &field)
        {
            v_AddTraceQuadPhysToField(Fwd, Bwd, field);
        }

        inline void ExpList::FillBwdWithBound(
            const Array<OneD, const NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            v_FillBwdWithBound(Fwd, Bwd);
        }

        inline void ExpList::FillBwdWithBoundDeriv(
            const int                          Dir,
            const Array<OneD, const NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            v_FillBwdWithBoundDeriv(Dir, Fwd, Bwd);
        }

        inline void ExpList::FillBwdWithBwdWeight(
            Array<OneD,       NekDouble> &weightave,
            Array<OneD,       NekDouble> &weightjmp)
        {
            v_FillBwdWithBwdWeight(weightave, weightjmp);
        }

        inline void ExpList::PeriodicBwdCopy(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd)
        {
            v_PeriodicBwdCopy(Fwd, Bwd);
        }

        inline const std::vector<bool> &ExpList::GetLeftAdjacentFaces(void) const
        {
            return v_GetLeftAdjacentFaces();
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

        inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>
            &ExpList::UpdateBndConditions()
        {
            return v_UpdateBndConditions();
        }

        inline void ExpList::EvaluateBoundaryConditions(
            const NekDouble   time,
            const std::string varName,
            const NekDouble   x2_in,
            const NekDouble   x3_in)
        {
            v_EvaluateBoundaryConditions(time, varName, x2_in, x3_in);
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
                                Array<OneD,      NekDouble> &outarray)
        {
            v_GeneralMatrixOp(gkey,inarray,outarray);
        }


        inline void ExpList::SetUpPhysNormals()
        {
            v_SetUpPhysNormals();
        }

        inline void ExpList::GetBoundaryToElmtMap( Array<OneD, int> &ElmtID,
                                            Array<OneD,int> &EdgeID)
        {
            v_GetBoundaryToElmtMap(ElmtID,EdgeID);
        }

        inline void ExpList::GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays)
        {
            v_GetBndElmtExpansion(i, result, DeclareCoeffPhysArrays);
        }

        inline void ExpList::ExtractElmtToBndPhys(int i,
                                                  const Array<OneD, NekDouble> &elmt,
                                                  Array<OneD, NekDouble> &boundary)
        {
            v_ExtractElmtToBndPhys(i, elmt, boundary);
        }

        inline void ExpList::ExtractPhysToBndElmt(int i,
                            const Array<OneD, const NekDouble> &phys,
                            Array<OneD, NekDouble> &bndElmt)
        {
            v_ExtractPhysToBndElmt(i, phys, bndElmt);
        }

        inline void ExpList::ExtractPhysToBnd(int i,
                            const Array<OneD, const NekDouble> &phys,
                            Array<OneD, NekDouble> &bnd)
        {
            v_ExtractPhysToBnd(i, phys, bnd);
        }

        inline void ExpList::GetBoundaryNormals(int i,
                            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            v_GetBoundaryNormals(i, normals);
        }

        const static Array<OneD, ExpListSharedPtr> NullExpListSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif // EXPLIST_H

