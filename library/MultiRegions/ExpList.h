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

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/MultiRegions.hpp>
#include <StdRegions/StdExpansion.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/GlobalOptimizationParameters.h>
#include <boost/enable_shared_from_this.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>

#include <LibUtilities/Communication/Transposition.h>

#include <tinyxml/tinyxml.h>

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
        
        MultiRegions::Direction const DirCartesianMap[] =
            {
                eX,
                eY,
                eZ
            }; 
    
        /// A map between global matrix keys and their associated block
        /// matrices.
        typedef map<GlobalMatrixKey,DNekScalBlkMatSharedPtr> BlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef boost::shared_ptr<BlockMatrixMap> BlockMatrixMapShPtr;
			       

        /// Base class for all multi-elemental spectral/hp expansions.
        class ExpList: public boost::enable_shared_from_this<ExpList>
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

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList(const ExpList &in, const bool DeclareCoeffPhysArrays = true);

            /// The default destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList();
			
            ////
            //virtual boost::shared_ptr<ExpList> do_clone(void);
            
            /// Copy coefficients from concatenated list to expansion list.
            MULTI_REGIONS_EXPORT void PutCoeffsInToElmtExp(void);

            /// Copy coefficients from expansion list to concatenated list.
            MULTI_REGIONS_EXPORT void PutElmtExpInToCoeffs(void);

            /// Copy one elements coefficients from the concatenated list
            /// to the expansion list.
            MULTI_REGIONS_EXPORT void PutCoeffsInToElmtExp(int eid);

            /// Copy one elements coefficients from the expansion list to
            /// the concatenated list.
            MULTI_REGIONS_EXPORT void PutElmtExpInToCoeffs(int eid);

            /// Copy physical data from \a m_phys to expansion list.
            MULTI_REGIONS_EXPORT void PutPhysInToElmtExp(void);

            /// Copy physical data from given array to expansion list.
            MULTI_REGIONS_EXPORT void PutPhysInToElmtExp(Array<OneD, const NekDouble> &in);

            /// Copy expansion list physical data to given array.
            MULTI_REGIONS_EXPORT void PutElmtExpInToPhys(Array<OneD,NekDouble> &out);

            /// Copy expansion list physical data from one element to array.
            MULTI_REGIONS_EXPORT void PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out);

            /// Returns the total number of local degrees of freedom
            /// \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
            inline int GetNcoeffs(void) const;

            // Returns the total number of local degrees of freedom
            // for element eid
            MULTI_REGIONS_EXPORT int GetNcoeffs(const int eid) const;
			

            /// Evaulates the maximum number of modes in the elemental basis
            /// order over all elements
            inline int EvalBasisNumModesMax(void) const;

            /// Returns the vector of the number of modes in the elemental
            /// basis order over all elements.
            MULTI_REGIONS_EXPORT const Array<OneD,int> EvalBasisNumModesMaxPerExp(void) const;

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
            MULTI_REGIONS_EXPORT NekDouble PhysIntegral(const Array<OneD, const NekDouble> &inarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all \emph{local}
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            inline void   IProductWRTBase_IterPerExp(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray,
                                CoeffState coeffstate = eLocal);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to the derivative (in
            /// direction \param dir) of all \emph{local} expansion modes
            /// \f$\phi_n^e(\boldsymbol{x})\f$.
            MULTI_REGIONS_EXPORT void   IProductWRTDerivBase(const int dir,
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            /// This function elementally evaluates the forward transformation
            /// of a function \f$u(\boldsymbol{x})\f$ onto the global
            /// spectral/hp expansion.
            inline void  FwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray,
											 Array<OneD,NekDouble> &outarray);

            ///
            inline void FwdTrans(const Array<OneD, const NekDouble> &inarray,
								Array<OneD,       NekDouble> &outarray,
                                CoeffState coeffstate = eLocal);

            /// This function elementally mulplies the coefficient space of
            /// Sin my the elemental inverse of the mass matrix.
            MULTI_REGIONS_EXPORT void  MultiplyByElmtInvMass (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void MultiplyByInvMassMatrix(
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                CoeffState coeffstate = eLocal);

            /// Smooth a field across elements
            inline void SmoothField(Array<OneD,NekDouble> &field);

            /// Solve helmholtz problem
            inline void HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff = StdRegions::NullVarCoeffMap,
                    const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);

            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionDiffusionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       CoeffState coeffstate = eLocal, 
                       const Array<OneD, const NekDouble>&
                       dirForcing = NullNekDouble1DArray);


            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       CoeffState coeffstate = eLocal,
                       const Array<OneD, const NekDouble>&
                       dirForcing = NullNekDouble1DArray);

            ///
            MULTI_REGIONS_EXPORT void FwdTrans_BndConstrained(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray);


            /// This function elementally evaluates the backward transformation
            /// of the global spectral/hp element expansion.
            inline void BwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray,
											 Array<OneD,NekDouble> &outarray);

            ///
            inline void BwdTrans (const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray,
                                  CoeffState coeffstate = eLocal);

            /// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(Array<OneD, NekDouble> &coord_0,
                                  Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                                  Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);
			
			/// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(NekDouble &x, NekDouble &y, NekDouble &z);
			
			/// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoord(Array<OneD, NekDouble> &coords);
			
			// Homogeneous transforms
            inline void HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                            Array<OneD, NekDouble> &outarray, 
                                            
                                            CoeffState coeffstate = eLocal,
                                            bool Shuff = true,
                                            bool UnShuff = true);
            
            inline void HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                            Array<OneD, NekDouble> &outarray, 
                                            CoeffState coeffstate = eLocal,
                                            bool Shuff = true,
                                            bool UnShuff = true);
            
            inline void DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                      const Array<OneD, NekDouble> &inarray2,
                                      Array<OneD, NekDouble> &outarray, 
                                      CoeffState coeffstate = eLocal);
			
            inline void GetBCValues(Array<OneD, NekDouble> &BndVals, 
                                    const Array<OneD, NekDouble> &TotField, 
                                    int BndID);
            
            inline void NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
                                                  Array<OneD, const NekDouble> &V2,
                                                  Array<OneD, NekDouble> &outarray,
                                                  int BndID);
			
            /// This function calculates Surface Normal vector of a smooth
            /// manifold.
            MULTI_REGIONS_EXPORT void GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,
                                  const int k);

            /// Populate tangents vector with tangents from each element.
            MULTI_REGIONS_EXPORT void GetTangents(
                             Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tangents);

            /// Apply geometry information to each expansion.
            MULTI_REGIONS_EXPORT void ApplyGeomInfo();

            /// This function writes the spectral/hp element solution to the
            /// file \a out.
            MULTI_REGIONS_EXPORT void WriteToFile(std::ofstream &out,
                             OutputFormat format = eTecplot,
                             std::string var = "v");

            void WriteTecplotHeader(std::ofstream &outfile,
                                    std::string var = "v")
            {
                v_WriteTecplotHeader(outfile,var);
            }

            void WriteTecplotZone(std::ofstream &outfile, int expansion)
            {
                v_WriteTecplotZone(outfile,expansion);
            }

            void WriteTecplotField(std::ofstream &outfile, int expansion)
            {
                v_WriteTecplotField(outfile,expansion);
            }

            MULTI_REGIONS_EXPORT void WriteVtkHeader(std::ofstream &outfile);
            MULTI_REGIONS_EXPORT void WriteVtkFooter(std::ofstream &outfile);

            void WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
            {
                v_WriteVtkPieceHeader(outfile, expansion);
            }

            MULTI_REGIONS_EXPORT void WriteVtkPieceFooter(std::ofstream &outfile, int expansion);

            void WriteVtkPieceData  (std::ofstream &outfile, int expansion,
                                     std::string var = "v")
            {
                v_WriteVtkPieceData(outfile, expansion, var);
            }

            MULTI_REGIONS_EXPORT void ReadFromFile(std::ifstream &in,
                              OutputFormat format = eTecplot);

            /// This function returns the dimension of the coordinates of the
            /// element \a eid.
            // inline
            MULTI_REGIONS_EXPORT int GetCoordim(int eid);

            /// Set the \a i th coefficiient in \a m_coeffs to value \a val
            inline void SetCoeff(int i, NekDouble val);
            
            /// Set the coefficiient in \a m_coeffs to value \a val (0D Exapnsion)
            inline void SetCoeff(NekDouble val);
            
            /// Set the physical value in \a m_coeffs to value \a val (0D Exapnsion)
            inline void SetPhys(NekDouble val);
            
            inline const SpatialDomains::VertexComponentSharedPtr &GetGeom(void) const;
            
            inline const SpatialDomains::VertexComponentSharedPtr &GetVertex(void) const;

            /// Set the \a i th coefficiient in  #m_coeffs to value \a val
            inline void SetCoeffs(int i, NekDouble val);

            /// Set the  #m_coeffs array to inarray
            inline void SetCoeffsArray(Array<OneD, NekDouble> &inarray);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline const Array<OneD, const NekDouble> &GetCoeffs() const;

            /// Impose Dirichlet Boundary Conditions onto Array
            inline void ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            /// Put the coefficients into global ordering using m_coeffs 
            inline void LocalToGlobal(void);

            /// Put the coefficients into local ordering and place in m_coeffs
            inline void GlobalToLocal(void);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeff(int i);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeffs(int i);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD, const NekDouble> &GetPhys()  const;

            /// This function calculates the \f$L_\infty\f$ error of the global
            /// spectral/hp element approximation.
            MULTI_REGIONS_EXPORT NekDouble Linf (const Array<OneD, const NekDouble> &soln);

            /// This function calculates the \f$L_\infty\f$ error of the global
            /// spectral/hp element approximation.
            MULTI_REGIONS_EXPORT NekDouble Linf (void);

            /// This function calculates the \f$L_2\f$ error with
            /// respect to soln of the global
            /// spectral/hp element approximation.
            NekDouble L2 (const Array<OneD, const NekDouble> &soln)
            {
                return v_L2(soln);
            }

            /// This function calculates the \f$L_2\f$ measure of the global
            /// spectral/hp element approximation.
            NekDouble L2 (void)
            {
                return v_L2();
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
            /// associaed with the homgeneous expansion.
            LibUtilities::TranspositionSharedPtr GetTransposition(void)
            {
                return v_GetTransposition();
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
            void PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble>  &outarray)
            {
                v_PhysInterp1DScaled(scale, inarray,outarray);                
            }

            /// This function Galerkin projects the physical space points in
            /// \a inarray to \a outarray where inarray is assumed to
            /// be defined in the expansion but where the number of
            /// points are rescaled by \a 1DScale
            void PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
            {
                v_PhysGalerkinProjection1DScaled(scale, inarray, outarray);
            } 

            /// Calculates the \f$H^1\f$ error of the global spectral/hp
            /// element approximation.
            MULTI_REGIONS_EXPORT NekDouble H1 (const Array<OneD, const NekDouble> &soln);
            
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
            inline const boost::shared_ptr<StdRegions::StdExpansionVector> GetExp() const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion of the \f$n^{\mathrm{th}}\f$ element.
            inline StdRegions::StdExpansionSharedPtr& GetExp(int n) const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion containing the arbitrary point given by \a gloCoord.
            MULTI_REGIONS_EXPORT StdRegions::StdExpansionSharedPtr& GetExp(
                                                                           const Array<OneD, const NekDouble> &gloCoord);

            /// This function returns the index of the local elemental
            /// expansion containing the arbitrary point given by \a gloCoord.
            MULTI_REGIONS_EXPORT int GetExpIndex(const Array<OneD, const NekDouble> &gloCoord, NekDouble tol = 0.0);

            /// Get the start offset position for a global list of #m_coeffs
            /// correspoinding to element n.
            inline int GetCoeff_Offset(int n) const;

            /// Get the start offset position for a global list of m_phys
            /// correspoinding to element n.
            inline int GetPhys_Offset(int n) const;

            /// Get the element id associated with the n th
            /// consecutive block of data in  #m_phys and #m_coeffs
            inline int GetOffset_Elmt_Id(int n) const;

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline Array<OneD, NekDouble> &UpdateCoeffs();

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            inline Array<OneD, NekDouble> &UpdatePhys();

            inline void PhysDeriv(Direction edir, 
                                  const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out_d);	
	    
            /// This function discretely evaluates the derivative of a function
            /// \f$f(\boldsymbol{x})\f$ on the domain consisting of all
            /// elements of the expansion.
            inline void PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out_d0,
                                  Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                                  Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);
            
            inline void PhysDeriv(const int dir,
                                  const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out_d);
            
            
            // functions associated with DisContField
            inline const Array<OneD, const  boost::shared_ptr<ExpList> > &GetBndCondExpansions();
            
            inline boost::shared_ptr<ExpList> &UpdateBndCondExpansion(int i);
            
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
            inline boost::shared_ptr<ExpList> &GetTrace();
            
            inline boost::shared_ptr<ExpList> &GetTrace(int i);
            
            inline boost::shared_ptr<AssemblyMapDG> &GetTraceMap(void);
            
            inline void GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals);

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

            inline void GetFwdBwdTracePhys( Array<OneD,NekDouble> &Fwd,
                                            Array<OneD,NekDouble> &Bwd);

            inline void GetFwdBwdTracePhys(
                                           const Array<OneD,const NekDouble> &field,
                                           Array<OneD,NekDouble> &Fwd,
                                           Array<OneD,NekDouble> &Bwd);

            inline void ExtractTracePhys(
                                         Array<OneD,NekDouble> &outarray);

            inline void ExtractTracePhys(
                                         const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray);

            inline const Array<OneD, const SpatialDomains
                ::BoundaryConditionShPtr>& GetBndConditions();

            inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>& UpdateBndConditions();

            inline void EvaluateBoundaryConditions(
                                                   const NekDouble time = 0.0, 
                                                   const NekDouble = NekConstants::kNekUnsetDouble, 
                                                   const NekDouble = NekConstants::kNekUnsetDouble);


            // Routines for continous matrix solution
            /// This function calculates the result of the multiplication of a
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            inline void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                        const Array<OneD,const NekDouble> &inarray,
                                        Array<OneD,      NekDouble> &outarray,
                                        CoeffState coeffstate = eLocal);

            MULTI_REGIONS_EXPORT void GeneralMatrixOp_IterPerExp(
                                                                 const GlobalMatrixKey      &gkey,
                                                                 const Array<OneD,const NekDouble> &inarray,
                                                                 Array<OneD,      NekDouble> &outarray);

            inline void SetUpPhysNormals();

            inline void SetUpPhysTangents(const StdRegions::StdExpansionVector &locexp);
 	                

            inline void SetUpTangents();

            inline void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                             Array<OneD,int> &EdgeID);

            MULTI_REGIONS_EXPORT void  GeneralGetFieldDefinitions(
                                 std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef, 
                                 int NumHomoDir = 0, 
                                 Array<OneD, LibUtilities::BasisSharedPtr> &HomoBasis = LibUtilities::NullBasisSharedPtr1DArray, 
                                 std::vector<NekDouble> &HomoLen = LibUtilities::NullNekDoubleVector,
                                 std::vector<unsigned int> &HomoZIDs = LibUtilities::NullUnsignedIntVector,
                                 std::vector<unsigned int> &HomoYIDs = LibUtilities::NullUnsignedIntVector);
            
            const NekOptimize::GlobalOptParamSharedPtr &GetGlobalOptParam(void)
            {
                return m_globalOptParam;
            }

            map<int, RobinBCInfoSharedPtr> GetRobinBCInfo()
            {
                return v_GetRobinBCInfo();
            }

            void GetPeriodicEdges(
                                  vector<map<int,int> > &periodicVertices,
                                  map<int,int>          &periodicEdges)
            {
                v_GetPeriodicEdges(periodicVertices, periodicEdges);
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
            void AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                                 std::vector<NekDouble> &fielddata)
            {
                v_AppendFieldData(fielddef,fielddata);
            }

            
            /// Append the data in coeffs listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef,
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
            MULTI_REGIONS_EXPORT  void ExtractCoeffsToCoeffs(const boost::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs);
			
			
            //Extract data in fielddata into the m_coeffs_list for the 3D stability analysis (base flow is 2D)
            MULTI_REGIONS_EXPORT void ExtractDataToCoeffs(
                                       LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                                       std::vector<NekDouble> &fielddata,
                                       std::string &field,
                                       Array<OneD, NekDouble> &coeffs);
			

            /// Returns a shared pointer to the current object.
            boost::shared_ptr<ExpList> GetSharedThisPtr()
            {
                return shared_from_this();
            }

            /// Returns the session object
            boost::shared_ptr<LibUtilities::SessionReader> GetSession()
            {
                return m_session;
            }

            /// Returns the comm object
            boost::shared_ptr<LibUtilities::Comm> GetComm()
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

            boost::shared_ptr<ExpList> &GetPlane(int n)
            {
                return v_GetPlane(n);
            }
            
        protected:
            boost::shared_ptr<DNekMat> GenGlobalMatrixFull(
                                                           const GlobalLinSysKey &mkey,
                                                           const boost::shared_ptr<AssemblyMapCG> &locToGloMap);

            /// Communicator
            LibUtilities::CommSharedPtr m_comm;

            /// Session
            LibUtilities::SessionReaderSharedPtr m_session;

            /// Mesh associated with this expansion list.
            SpatialDomains::MeshGraphSharedPtr m_graph;

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

            /// Array containing the element id #m_offset_elmt_id[n]
            /// that the n^th consecutive block of data in #m_coeffs
            /// and #m_phys is associated, i.e. for an array of
            /// constant expansion size and single shape elements
            /// m_phys[n*m_npoints] is the data related to
            /// m_exp[m_offset_elmt_id[n]];
            Array<OneD, int>  m_offset_elmt_id;

            NekOptimize::GlobalOptParamSharedPtr m_globalOptParam;

            BlockMatrixMapShPtr  m_blockMat;
			
            //@todo should this be in ExpList or ExpListHomogeneous1D.cpp
            // it's a bool which determine if the expansion is in the wave space (coefficient space)
            // or not
            bool m_WaveSpace;
			
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
            boost::shared_ptr<GlobalMatrix>  GenGlobalMatrix(
                                                             const GlobalMatrixKey &mkey,
                                                             const boost::shared_ptr<AssemblyMapCG> &locToGloMap);


            void GlobalEigenSystem(const boost::shared_ptr<DNekMat> &Gmat,
                                   Array<OneD, NekDouble> &EigValsReal,
                                   Array<OneD, NekDouble> &EigValsImag,
                                   Array<OneD, NekDouble> &EigVecs
                                   = NullNekDouble1DArray);


            /// This operation constructs the global linear system of type \a
            /// mkey.
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSys(
                                                             const GlobalLinSysKey &mkey,
                                                             const boost::shared_ptr<AssemblyMapCG> &locToGloMap);

            /// Generate a GlobalLinSys from information provided by the key
            /// "mkey" and the mapping provided in LocToGloBaseMap.
            boost::shared_ptr<GlobalLinSys> GenGlobalBndLinSys(
                                                               const GlobalLinSysKey     &mkey,
                                                               const AssemblyMapSharedPtr &locToGloMap);

            void ReadGlobalOptimizationParameters()
            {
                v_ReadGlobalOptimizationParameters();
            }

            // Virtual prototypes

            virtual int v_GetNumElmts(void)
            {
                return (*m_exp).size();
            }

            virtual const Array<OneD,const boost::shared_ptr<ExpList> > &v_GetBndCondExpansions(void);

            virtual boost::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i);
            
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

            virtual boost::shared_ptr<ExpList> &v_GetTrace();
			
            virtual boost::shared_ptr<ExpList> &v_GetTrace(int i);

            virtual boost::shared_ptr<AssemblyMapDG> &v_GetTraceMap();

            virtual void v_GetNormals(
                                      Array<OneD, Array<OneD, NekDouble> > &normals);

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

            virtual void v_ExtractTracePhys(
                                            Array<OneD,NekDouble> &outarray);

            virtual void v_ExtractTracePhys(
                                            const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray);

            virtual void v_MultiplyByInvMassMatrix(
                                                   const Array<OneD,const NekDouble> &inarray,
                                                   Array<OneD,      NekDouble> &outarray,
                                                   CoeffState coeffstate);

            virtual void v_HelmSolve(
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,       NekDouble> &outarray,
                                     const FlagList &flags,
                                     const StdRegions::ConstFactorMap &factors,
                                     const StdRegions::VarCoeffMap &varcoeff,
                                     const Array<OneD, const NekDouble> &dirForcing);

            virtual void v_LinearAdvectionDiffusionReactionSolve(
                                                                 const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                                 const Array<OneD, const NekDouble> &inarray,
                                                                 Array<OneD, NekDouble> &outarray,
                                                                 const NekDouble lambda,
                                                                 CoeffState coeffstate = eLocal, 
                                                                 const Array<OneD, const NekDouble>&
                                                                 dirForcing = NullNekDouble1DArray);

            virtual void v_LinearAdvectionReactionSolve(
                                                        const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                        const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD, NekDouble> &outarray,
                                                        const NekDouble lambda,
                                                        CoeffState coeffstate = eLocal, 
                                                        const Array<OneD, const NekDouble>&
                                                        dirForcing = NullNekDouble1DArray);

            // wrapper functions about virtual functions
            virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            virtual void v_LocalToGlobal(void);

            virtual void v_GlobalToLocal(void);

            virtual void v_BwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    CoeffState coeffstate);
			
            virtual void v_BwdTrans_IterPerExp(const Array<OneD,const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray);
	    
            virtual void v_FwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    CoeffState coeffstate);

            virtual void v_FwdTrans_IterPerExp(
                                    const Array<OneD,const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray);

            virtual void v_SmoothField(Array<OneD,NekDouble> &field);

            virtual void v_IProductWRTBase(const Array<OneD,const NekDouble> &inarray,Array<OneD,      NekDouble> &outarray, CoeffState coeffstate);
			
            virtual void v_IProductWRTBase_IterPerExp(const Array<OneD,const NekDouble> &inarray,  Array<OneD,      NekDouble> &outarray);
			
            virtual void v_SetUpPhysTangents(const StdRegions::StdExpansionVector &locexp);
            
            virtual void v_GeneralMatrixOp(
                                           const GlobalMatrixKey             &gkey,
                                           const Array<OneD,const NekDouble> &inarray,
                                           Array<OneD,      NekDouble> &outarray,
                                           CoeffState coeffstate);
            
            virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                                     Array<OneD, NekDouble> &coord_1,
                                     Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);
			
            virtual void v_GetCoords(NekDouble &x,NekDouble &y,NekDouble &z);
            
            virtual void v_GetCoord(Array<OneD, NekDouble> &coords);

            virtual void v_SetCoeff(NekDouble val);
            
            virtual void v_SetPhys(NekDouble val);
            
            virtual const SpatialDomains::VertexComponentSharedPtr &v_GetGeom(void) const;
            
            virtual const SpatialDomains::VertexComponentSharedPtr &v_GetVertex(void) const;
            
            virtual void v_PhysDeriv(
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1, 
                                     Array<OneD, NekDouble> &out_d2);
            
            virtual void v_PhysDeriv(const int dir,
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d);
            
            virtual void v_PhysDeriv(Direction edir, 
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d);
            
            virtual void v_HomogeneousFwdTrans(
                                               const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD, NekDouble> &outarray, 
                                               CoeffState coeffstate = eLocal,
                                               bool Shuff = true,
                                               bool UnShuff = true);
            
            virtual void v_HomogeneousBwdTrans(
                                               const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD, NekDouble> &outarray, 
                                               CoeffState coeffstate = eLocal,
                                               bool Shuff = true,
                                               bool UnShuff = true);
            
            virtual void v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                         const Array<OneD, NekDouble> &inarray2,
                                         Array<OneD, NekDouble> &outarray, 
                                         CoeffState coeffstate = eLocal);
            
            virtual void v_GetBCValues(Array<OneD, NekDouble> &BndVals, 
                                       const Array<OneD, NekDouble> &TotField, 
                                       int BndID);
            
            virtual void v_NormVectorIProductWRTBase(
                                                     Array<OneD, const NekDouble> &V1,
                                                     Array<OneD, const NekDouble> &V2,
                                                     Array<OneD, NekDouble> &outarray,
                                                     int BndID);
            
            virtual void v_SetUpPhysNormals();
            
            virtual void v_SetUpTangents();
            
            virtual void v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                                Array<OneD,int> &EdgeID);

            virtual void v_ReadGlobalOptimizationParameters();

            virtual std::vector<LibUtilities::FieldDefinitionsSharedPtr> v_GetFieldDefinitions(void);

            virtual void  v_GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef);


            virtual void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata);

            virtual void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractDataToCoeffs(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field,
                                               Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractCoeffsToCoeffs(const boost::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs);
			
            virtual void v_WriteTecplotHeader(std::ofstream &outfile,
                                              std::string var = "v");
            virtual void v_WriteTecplotZone(std::ofstream &outfile,
                                            int expansion);
            virtual void v_WriteTecplotField(std::ofstream &outfile,
                                             int expansion);

            virtual void v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion);
            virtual void v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                             std::string var);

            virtual NekDouble v_L2(void);
            virtual NekDouble v_L2(const Array<OneD, const NekDouble> &soln);
            
            virtual Array<OneD, const NekDouble> v_HomogeneousEnergy(void);
            virtual LibUtilities::TranspositionSharedPtr v_GetTransposition(void);
            virtual Array<OneD, const unsigned int> v_GetZIDs(void);
            virtual Array<OneD, const unsigned int> v_GetYIDs(void);
            
            // 1D Scaling and projection
            virtual void v_PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray);
            
            virtual void v_PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray);
        
            // Utility function for a common case of retrieving a
            // BoundaryCondition from a boundary condition collection.
            MULTI_REGIONS_EXPORT
                static SpatialDomains::BoundaryConditionShPtr GetBoundaryCondition(const SpatialDomains::BoundaryConditionCollection& collection,
                                                                                   unsigned int index, const std::string& variable);
        
        private:
            
            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &v_GetBndConditions();
            
            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions();

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0, 
                                                      const NekDouble x2_in = NekConstants::kNekUnsetDouble,
                                                      const NekDouble x3_in = NekConstants::kNekUnsetDouble);
            
            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo(void);
            
            
            virtual void v_GetPeriodicEdges(
                vector<map<int,int> > &periodicVertices,
                map<int,int>          &periodicEdges);

            // Homogeneous direction wrapper functions. 
            virtual LibUtilities::BasisSharedPtr  v_GetHomogeneousBasis(void)
            {
                ASSERTL0(false,
                         "This method is not defined or valid for this class type");
                return LibUtilities::NullBasisSharedPtr; 
            }

            // wrapper function to set viscosity for Homo1D expansion
            virtual void v_SetHomo1DSpecVanVisc(Array<OneD, NekDouble> visc)
            {
                ASSERTL0(false,
                         "This method is not defined or valid for this class type");
            }


            virtual boost::shared_ptr<ExpList> &v_GetPlane(int n);
        };


        /// Shared pointer to an ExpList object.
        typedef boost::shared_ptr<ExpList>      ExpListSharedPtr;
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
            unsigned int i;
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
            ASSERTL0(inarray.num_elements() == m_npoints,
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
                                             Array<OneD, NekDouble> &outarray,
                                             CoeffState coeffstate)
        {
            v_IProductWRTBase(inarray,outarray, coeffstate);
        }

		/**
         *
         */
        inline void ExpList::IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray,
														Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase_IterPerExp(inarray,outarray);
        }

        /**
         *
         */
        inline void ExpList::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                      CoeffState coeffstate)
        {
            v_FwdTrans(inarray,outarray,coeffstate);
        }
		
		/**
         *
         */
        inline void ExpList::FwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray,
												  Array<OneD,NekDouble> &outarray)
        {
            v_FwdTrans_IterPerExp(inarray,outarray);
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
        inline void ExpList::BwdTrans (const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       CoeffState coeffstate)
        {
            v_BwdTrans(inarray,outarray,coeffstate);
        }
		
		/**
         *
         */
        inline void ExpList::BwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray)
        {
            v_BwdTrans_IterPerExp(inarray,outarray);
        }


        /**
         *
         */
        inline void ExpList::MultiplyByInvMassMatrix(const Array<OneD,const NekDouble> &inarray,
                                                     Array<OneD,      NekDouble> &outarray,
                                                     CoeffState coeffstate)
        {
            v_MultiplyByInvMassMatrix(inarray,outarray,coeffstate);
        }

        /**
         *
         */
        inline void ExpList::HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            v_HelmSolve(inarray, outarray, flags, factors, varcoeff, dirForcing);
        }


        /**
         *
         */
        inline void ExpList::LinearAdvectionDiffusionReactionSolve(
                                                                   const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                                   const Array<OneD, const NekDouble> &inarray,
                                                                   Array<OneD, NekDouble> &outarray,
                                                                   const NekDouble lambda,
                                                                   CoeffState coeffstate,
                                                                   const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionDiffusionReactionSolve(velocity,inarray, outarray, lambda, coeffstate,dirForcing);
        }
        
        inline void ExpList::LinearAdvectionReactionSolve(
                                                          const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                          const Array<OneD, const NekDouble> &inarray,
                                                          Array<OneD, NekDouble> &outarray,
                                                          const NekDouble lambda,
                                                          CoeffState coeffstate,
                                                          const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionReactionSolve(velocity,inarray, outarray, lambda, coeffstate,dirForcing);
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
        inline void ExpList::SetCoeff(NekDouble val)
		
        {
            v_SetCoeff(val);
        }
		
		/**
         *
         */
        inline const SpatialDomains::VertexComponentSharedPtr &ExpList::GetGeom(void) const
        {
            return v_GetGeom();
        }
	
        /**
         *
         */
        inline const SpatialDomains::VertexComponentSharedPtr &ExpList::GetVertex(void) const
        {
            return v_GetVertex();
        }
	
		
        /**
         *
         */
        inline void ExpList::SetPhys(NekDouble val)
            
        {
            v_SetPhys(val);
        }
	
        /**
         *
         */
        inline void ExpList::GetCoords(NekDouble &x,NekDouble &y,NekDouble &z)
        {
            v_GetCoords(x,y,z);
        }
	
        inline void ExpList::GetCoord(Array<OneD, NekDouble> &coords)
        {
            v_GetCoord(coords);
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
        inline void ExpList::PhysDeriv(const int dir,
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(dir,inarray,out_d);
        }
        
        inline void ExpList::PhysDeriv(Direction edir,
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(edir, inarray,out_d);
        }		
	
        /**
         *
         */
        inline void ExpList::HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                 Array<OneD, NekDouble> &outarray, 
                                                 CoeffState coeffstate,
                                                 bool Shuff,
                                                 bool UnShuff)
        {
            v_HomogeneousFwdTrans(inarray,outarray,coeffstate,Shuff,UnShuff);
        }
	
        /**
         *
         */
        inline void ExpList::HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                 Array<OneD, NekDouble> &outarray, 
                                                 CoeffState coeffstate,
                                                 bool Shuff,
                                                 bool UnShuff)
        {
            v_HomogeneousBwdTrans(inarray,outarray,coeffstate,Shuff,UnShuff);
        }
	
        /**
         *
         */
        inline void ExpList::DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                           const Array<OneD, NekDouble> &inarray2,
                                           Array<OneD, NekDouble> &outarray, 
                                           CoeffState coeffstate)
        {
            v_DealiasedProd(inarray1,inarray2,outarray,coeffstate);
        }
		
        /**
         *
         */
        inline void ExpList::GetBCValues(Array<OneD, NekDouble> &BndVals, 
                                         const Array<OneD, NekDouble> &TotField, 
                                         int BndID)
        {
            v_GetBCValues(BndVals,TotField,BndID);
        }
	
        /**
         *
         */
        inline void ExpList::NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
                                                       Array<OneD, const NekDouble> &V2,
                                                       Array<OneD, NekDouble> &outarray,
                                                       int BndID)
        {
            v_NormVectorIProductWRTBase(V1,V2,outarray,BndID);
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
        
        inline void ExpList::ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            v_ImposeDirichletConditions(outarray);
        }

        inline void ExpList::LocalToGlobal(void)
        {
            v_LocalToGlobal();
        }
        
        inline void ExpList::GlobalToLocal(void)
        {
            v_GlobalToLocal();
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
        inline StdRegions::StdExpansionSharedPtr& ExpList::GetExp(int n) const
        {
            return (*m_exp)[n];
        }

        /**
         * @return  (A const shared pointer to) the local expansion vector #m_exp
         */
        inline const boost::shared_ptr<StdRegions::StdExpansionVector> ExpList::GetExp(void) const
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
         *
         */
        inline int ExpList::GetOffset_Elmt_Id(int n) const
        {
            return m_offset_elmt_id[n];
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
        inline const Array<OneD, const  boost::shared_ptr<ExpList> > &ExpList::GetBndCondExpansions()
        {
            return v_GetBndCondExpansions();
        }
        
        inline boost::shared_ptr<ExpList>  &ExpList::UpdateBndCondExpansion(int i)
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

        inline boost::shared_ptr<ExpList> &ExpList::GetTrace()
        {
            return v_GetTrace();
        }

        inline boost::shared_ptr<ExpList> &ExpList::GetTrace(int i)
        {
            return v_GetTrace(i);
        }
		
        inline boost::shared_ptr<AssemblyMapDG> &ExpList::GetTraceMap()
        {
            return v_GetTraceMap();
        }

        inline void ExpList::GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            v_GetNormals(normals);
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


        inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>
            &ExpList::UpdateBndConditions()
        {
            return v_UpdateBndConditions();
        }

        inline void ExpList::EvaluateBoundaryConditions(const NekDouble time,
                                                        const NekDouble x2_in,
                                                        const NekDouble x3_in)
        {
            v_EvaluateBoundaryConditions(time,x2_in,x3_in);
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
                                CoeffState coeffstate)
        {
            v_GeneralMatrixOp(gkey,inarray,outarray,coeffstate);
        }


        inline void ExpList::SetUpPhysNormals()
        {
            v_SetUpPhysNormals();
        }

        inline void ExpList::SetUpPhysTangents(
                                const StdRegions::StdExpansionVector &locexp)
        {
            v_SetUpPhysTangents(locexp);
        }
        
        inline void ExpList::SetUpTangents()
        {
            v_SetUpTangents();
        }

        inline void ExpList::GetBoundaryToElmtMap( Array<OneD, int> &ElmtID,
                                            Array<OneD,int> &EdgeID)
        {
            v_GetBoundaryToElmtMap(ElmtID,EdgeID);
        }

        const static Array<OneD, ExpListSharedPtr> NullExpListSharedPtrArray;
        
    } //end of namespace
} //end of namespace

#endif // EXPLIST_H

