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
//#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalOptimizationParameters.h>

#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/SegGeom.h>

#include <SpatialDomains/MeshGraph.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class GlobalLinSys;
        class LocalToGlobalC0ContMap;
        class LocalToGlobalBaseMap;
	class LocalToGlobalDGMap;
        class GenExpList1D;
	class ExpList1D;

        typedef map<GlobalMatrixKey,DNekScalBlkMatSharedPtr> BlockMatrixMap;
        typedef boost::shared_ptr<BlockMatrixMap> BlockMatrixMapShPtr; 

        /**
         * \brief This is the base class for all multi-elemental spectral/hp 
         * expansions.
         *
         * All multi-elemental expansions \f$u^{\delta}(\boldsymbol{x})\f$ can be 
         * considered as the assembly of the various elemental contributions. 
         * On a discrete level, this yields,
         * \f[u^{\delta}(\boldsymbol{x}_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and 
         * \f$N^{e}_m\f$ is the local elemental number of expansion modes.
         * As it is the lowest level class, it contains the definition of the 
         * common data and common routines to all multi-elemental expansions.
         */
        class ExpList
        {
        public:
            /**
             * \brief The default constructor.
             */
            ExpList();            

            /**
             * \brief The copy constructor.
             */
            ExpList(const ExpList &in);

            /**
             * \brief The default destructor.
             */
            virtual ~ExpList();
            
            void PutCoeffsInToElmtExp(void);
            void PutElmtExpInToCoeffs(void);

            void PutCoeffsInToElmtExp(int eid);
            void PutElmtExpInToCoeffs(int eid);
            
            void PutPhysInToElmtExp(void)
            {
                PutPhysInToElmtExp(m_phys);
            }

            void PutPhysInToElmtExp(Array<OneD, const NekDouble> &in);

            void PutElmtExpInToPhys(Array<OneD,NekDouble> &out);
            void PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out);
            
            /**
             * \brief This function returns the total number of local degrees of freedom 
             * \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
             */
            inline int GetNcoeffs(void) const
            {
                return m_ncoeffs;
            }
      
            /**
             * \brief This function evaulates the maximum number of
             * modes in the elemental basis order over all elements
             */
            inline int EvalBasisNumModesMax(void) const
            {
                int i;
                int returnval = 0;

                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    returnval = max(returnval,(*m_exp)[i]->EvalBasisNumModesMax());
                }

                return returnval;
            }


            /**
             * \brief This function returns the total number of quadrature points #m_npoints
             * \f$=Q_{\mathrm{tot}}\f$.
             */
            inline int GetTotPoints(void) const
            {
                return m_npoints;
            }

            /**
             * \brief This function returns the total number of quadrature points #m_npoints
             * \f$=Q_{\mathrm{tot}}\f$.
             */
            inline int GetNpoints(void) const
            {
                return m_npoints;
            }
      
            /**
             * \brief This function sets the transformed state #m_transState of the 
             * coefficient arrays.
             */
            inline void SetTransState(const TransState transState)
            {
                m_transState = transState;
            }

            /**
             * \brief This function returns the transformed state #m_transState of the 
             * coefficient arrays.
             */
            inline TransState GetTransState(void) const 
            {
                return m_transState; 
            }
      
            /**
             * \brief This function fills the array #m_phys
             *
             * This function fills the array \f$\boldsymbol{u}_l\f$, the evaluation of the 
             * expansion at the quadrature points (implemented as #m_phys), with the values 
             * of the array \a inarray.
             *
             * \param inarray The array containing the values where #m_phys should be 
             * filled with.
             */
            inline void SetPhys(const Array<OneD, const NekDouble> &inarray)
            {
                Vmath::Vcopy(m_npoints,&inarray[0],1,&m_phys[0],1);
                m_physState = true;
            }

            /**
             * \brief This function manually sets whether the array of physical values 
             * \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is filled or not.
             *
             * \param physState \a true (=filled) or \a false (=not filled).
             */
            inline void SetPhysState(const bool physState)
            {
                m_physState = physState;
            }

            /**
             * \brief This function indicates whether the array of physical values 
             * \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is filled or not.
             *
             * \return physState \a true (=filled) or \a false (=not filled).
             */
            inline bool GetPhysState(void) const
            {
                return m_physState;
            }
      
            /**
             * \brief This function integrates a function \f$f(\boldsymbol{x})\f$ over 
             * the domain consisting of all the elements of the expansion.
             *
             * The integration is evaluated locally, that is
             * \f[
             * \int f(\boldsymbol{x})d\boldsymbol{x}=\sum_{e=1}^{{N_{\mathrm{el}}}}
             * \left\{\int_{\Omega_e}f(\boldsymbol{x})d\boldsymbol{x}\right\},\f]
             * where the integration over the separate elements is done by the 
             * function StdRegions#StdExpansion#Integral, which discretely evaluates the 
             * integral using Gaussian quadrature.
             *
             * Note that the array #m_phys should be filled with the values of the 
             * function \f$f(\boldsymbol{x})\f$ at the quadrature points 
             * \f$\boldsymbol{x}_i\f$.
             *
             * \return the value of the discretely evaluated integral 
             * \f$\int f(\boldsymbol{x})d\boldsymbol{x}\f$.
             */
            NekDouble PhysIntegral (void);

            /**
             * \brief This function integrates a function \f$f(\boldsymbol{x})\f$ over the 
             * domain consisting of all the elements of the expansion.
             *
             * The integration is evaluated locally, that is
             * \f[
             * \int f(\boldsymbol{x})d\boldsymbol{x}=\sum_{e=1}^{{N_{\mathrm{el}}}}
             * \left\{\int_{\Omega_e}f(\boldsymbol{x})d\boldsymbol{x}\right\},\f]
             * where the integration over the separate elements is done by the function 
             * StdRegions#StdExpansion#Integral, which discretely evaluates the integral 
             * using Gaussian quadrature.
             * 
             * \param inarray An array of size \f$Q_{\mathrm{tot}}\f$ containing the values 
             * of the function \f$f(\boldsymbol{x})\f$ at the quadrature points 
             * \f$\boldsymbol{x}_i\f$.
             * \return the value of the discretely evaluated integral 
             * \f$\int f(\boldsymbol{x})d\boldsymbol{x}\f$.
             */
            NekDouble PhysIntegral(const Array<OneD, const NekDouble> &inarray);

            /**
             * \brief This function calculates the inner product of a function 
             * \f$f(\boldsymbol{x})\f$ with respect to all \emph{local} expansion modes 
             * \f$\phi_n^e(\boldsymbol{x})\f$.
             *
             * The operation is evaluated locally for every element by the function 
             * StdRegions#StdExpansion#IProductWRTBase.
             *
             * \param inarray An array of size \f$Q_{\mathrm{tot}}\f$
             * containing the values of the function
             * \f$f(\boldsymbol{x})\f$ at the quadrature points
             * \f$\boldsymbol{x}_i\f$.  \param outarray An array of
             * size \f$N_{\mathrm{eof}}\f$ used to store the result.
             */
            void   IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, 
                                                    Array<OneD,       NekDouble> &outarray);
            
            void   IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                         Array<OneD,       NekDouble> &outarray,
                                   bool  UseContCoeffs = false)
            {
                v_IProductWRTBase(inarray,outarray,UseContCoeffs);
            }
            
            /**
             * \brief This function calculates the inner product of a
             * function \f$f(\boldsymbol{x})\f$ with respect to the
             * derivative (in direction \param dir) of all
             * \emph{local} expansion modes
             * \f$\phi_n^e(\boldsymbol{x})\f$.
             *
             * The operation is evaluated locally for every element by the function 
             * StdRegions#StdExpansion#IProductWRTDerivBase.
             *
             * \param dir=0,1 is the direction in which the derivative
             * of the basis should be taken 
             * \param inarray An array of size \f$Q_{\mathrm{tot}}\f$
             * containing the values of the function
             * \f$f(\boldsymbol{x})\f$ at the quadrature points
             * \f$\boldsymbol{x}_i\f$.  
             * \param outarray An array of
             * size \f$N_{\mathrm{eof}}\f$ used to store the result.
             */
            void   IProductWRTDerivBase(const int dir, 
                                        const Array<OneD, const NekDouble> &inarray, 
                                              Array<OneD,       NekDouble> &outarray);

            /**
             * \brief This function elementally evaluates the forward
             * transformation of a function \f$u(\boldsymbol{x})\f$
             * onto the global spectral/hp expansion.
             *
             * Given a function \f$u(\boldsymbol{x})\f$ defined at the
             * quadrature points, this function determines the
             * transformed elemental coefficients \f$\hat{u}_n^e\f$
             * employing a discrete elemental Galerkin projection from
             * physical space to coefficient space. For each element,
             * the operation is evaluated locally by the function
             * StdRegions#StdExpansion#IproductWRTBase followed by a
             * call to #MultiRegions#MultiplyByElmtInvMass.
             *
             * \param inarray An array of size \f$Q_{\mathrm{tot}}\f$
             * containing the values of the function
             * \f$f(\boldsymbol{x})\f$ at the quadrature points
             * \f$\boldsymbol{x}_i\f$.  \param outarray The resulting
             * coefficients \f$\hat{u}_n^e\f$ will be stored in this
             * array of size \f$N_{\mathrm{eof}}\f$.
             */
            void   FwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,       NekDouble> &outarray);

            void   FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,       NekDouble> &outarray,
                            bool  UseContCoeffs = false)
            {
                v_FwdTrans(inarray,outarray,UseContCoeffs);
            }

            /**
             * \brief This function elementally mulplies the
             * coefficient space of Sin my the elemental inverse of
             * the mass matrix
             *
             * The coefficients of the function to be acted upon
             * should be contained in the \param inarray. The
             * resulting coefficients are stored in \param outarray
             *
             * \param inarray An array of size \f$N_{\mathrm{eof}}\f$
             * containing the inner product. 
             */
            void  MultiplyByElmtInvMass (const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD,       NekDouble> &outarray);

            void MultiplyByInvMassMatrix(const Array<OneD,const NekDouble> &inarray, 
                                               Array<OneD,      NekDouble> &outarray,
                                         bool  UseContCoeffs = false)
            {
                v_MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);                
            }

            void HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                           NekDouble lambda,
                           bool      UseContCoeffs = false,
                           const Array<OneD, const NekDouble>& dirForcing = NullNekDouble1DArray)
            {
                v_HelmSolve(inarray,outarray,lambda,UseContCoeffs,dirForcing);
            }

            /**
             * \brief 
             */
            void   FwdTrans_BndConstrained(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,       NekDouble> &outarray);

            /**
             * \brief This function elementally evaluates the backward
             * transformation of the global spectral/hp element
             * expansion.
             *
             * Given the elemental coefficients \f$\hat{u}_n^e\f$ of
             * an expansion, this function evaluates the spectral/hp
             * expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the
             * quadrature points \f$\boldsymbol{x}_i\f$. The operation
             * is evaluated locally by the elemental function
             * StdRegions#StdExpansion#BwdTrans.
             *
             * \param inarray An array of size \f$N_{\mathrm{eof}}\f$
             * containing the local coefficients \f$\hat{u}_n^e\f$.
             * \param outarray The resulting physical values at the
             * quadrature points \f$u^{\delta}(\boldsymbol{x}_i)\f$
             * will be stored in this array of size
             * \f$Q_{\mathrm{tot}}\f$.
             */
            void   BwdTrans_IterPerExp (const Array<OneD, const NekDouble> &inarray, 
                                              Array<OneD,       NekDouble> &outarray); 

            void   BwdTrans (const Array<OneD, const NekDouble> &inarray, 
                                   Array<OneD,       NekDouble> &outarray,
                             bool  UseContCoeffs = false)
            {
                v_BwdTrans(inarray,outarray,UseContCoeffs);
            }
       
            /**
             * \brief This function calculates the coordinates of all
             * the elemental quadrature points \f$\boldsymbol{x}_i\f$.
             *
             * The operation is evaluated locally by the elemental
             * function StdRegions#StdExpansion#GetCoords.
             *
             * \param coord_0 After calculation, the \f$x_1\f$
             * coordinate will be stored in this array.  \param
             * coord_1 After calculation, the \f$x_2\f$ coordinate
             * will be stored in this array.  \param coord_2 After
             * calculation, the \f$x_3\f$ coordinate will be stored in
             * this array.
             */
            void   GetCoords(Array<OneD, NekDouble> &coord_0,
                             Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                             Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

            /**
             * \brief This function calculates Surface Normal vector of a smooth manifold.
             *
             * The operation is evaluated locally by the elemental
             * function StdRegions#StdExpansion#GetSurfaceNormal.
             *
             */

            void  GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,const int k);
            /**
             * \brief This function calculates two tangential basis of a smooth manifold.
             *
             * The operation is evaluated locally by the elemental
             * function StdRegions#StdExpansion#GetTangentBasis.
             *
             * \param tanbasis_1 : tangential basis 1, \param tanbasis_2 : tangential basis 2.
             *  Normally tanbasis_1 /cdot tanbasis_2 = 0
             */

            void  GetTanBasis(Array<OneD,NekDouble> &tbasis1,
                             Array<OneD,NekDouble> &tbasis2,
                             const int k);
            
            /**
             * \brief This function writes the spectral/hp element
             * solution to the file \a out.
             *
             * The coordinates of the quadrature points, together with
             * the content of the array #m_phys, are written to the
             * file \a out.
             *
             * \param out The file to which the solution should be
             * written.
             */
            void   WriteToFile(std::ofstream &out, OutputFormat format = eTecplot, std::string var = "v");
    
            void ReadFromFile(std::ifstream &in, OutputFormat format = eTecplot);
            
            /**
             * \brief This function returns the dimension of the
             * coordinates of the element \a eid.
             *
             * \param eid The index of the element to be checked.
             * \return The dimension of the coordinates of the specific element.
             */
            inline int GetCoordim(int eid)
            {
                ASSERTL2(eid <= (*m_exp).size(),"eid is larger than number of elements");          
                return (*m_exp)[eid]->GetCoordim();
            }
      
            /**
             * \brief Set the \a i th coefficiient in \a m_coeffs to
             * value \a val
             *
             * \param i The index of m_coeffs to be set \param val The
             * value which m_coeffs[i] is to be set to.
             */
            inline void SetCoeff(int i, NekDouble val) 
            {
                m_coeffs[i] = val;
            }

            /**
             * \brief Set the \a i th coefficiient in  #m_coeffs to
             * value \a val
             *
             * \param i The index of #m_coeffs to be set \param val The
             * value which #m_coeffs[i] is to be set to.
             */
            inline void SetCoeffs(int i, NekDouble val) 
            {
                m_coeffs[i] = val;
            }

            /**
             * \brief This function returns (a reference to) the array
             * \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
             * containing all local expansion coefficients.
             *
             * As the function returns a constant reference to a
             * <em>const Array</em>, it is not possible to modify the
             * underlying data of the array #m_coeffs. In order to do
             * so, use the function #UpdateCoeffs instead.
             *
             * \return (A constant reference to) the array #m_coeffs.
             */
            inline const Array<OneD, const NekDouble> &GetCoeffs() const 
            {
                return m_coeffs;
            }

            
            inline const Array<OneD, const NekDouble> &GetContCoeffs() const 
            {
                return v_GetContCoeffs();
            }

            /**
             * \brief  Get the \a i th value  (coefficient) of #m_coeffs
             *
             * \param i The index of #m_coeffs to be returned
             * \return The NekDouble held in #m_coeffs[i].
             */
            inline NekDouble GetCoeff(int i) 
            {
                return m_coeffs[i];
            }

            /**
             * \brief  Get the \a i th value  (coefficient) of #m_coeffs
             *
             * \param i The index of #m_coeffs to be returned
             * \return The NekDouble held in #m_coeffs[i].
             */
            inline NekDouble GetCoeffs(int i) 
            {
                return m_coeffs[i];
            }

            /**
             * \brief This function returns (a reference to) the array
             * \f$\boldsymbol{u}_l\f$ (implemented as #m_phys)
             * containing the function
             * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
             * quadrature points.
             *
             * As the function returns a constant reference to a
             * <em>const Array</em> it is not possible to modify the
             * underlying data of the array #m_phys. In order to do
             * so, use the function #UpdatePhys instead.
             *
             * \return (A constant reference to) the array #m_phys.
             */
            inline const Array<OneD, const NekDouble> &GetPhys()  const
            {
                return m_phys;
            }

            /**
             * \brief This function calculates the \f$L_\infty\f$
             * error of the global spectral/hp element approximation.
             *
             * Given a spectral/hp approximation
             * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
             * quadrature points (which should be contained in
             * #m_phys), this function calculates the \f$L_\infty\f$
             * error of this approximation with respect to an exact
             * solution. The local distribution of the quadrature
             * points allows an elemental evaluation of this operation
             * through the functions StdRegions#StdExpansion#Linf.
             *
             * The exact solution, also evaluated at the quadrature
             * points, should be contained in the variable #m_phys of
             * the ExpList object \a Sol.
             *
             * \param Sol An ExpList, containing the discrete
             * evaluation of the exact solution at the quadrature
             * points in its array #m_phys.  \return The
             * \f$L_\infty\f$ error of the approximation.
             */
            NekDouble Linf (const Array<OneD, const NekDouble> &soln);

            /**
             * \brief This function calculates the \f$L_2\f$ error of
             * the global spectral/hp element approximation.
             *
             * Given a spectral/hp approximation
             * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
             * quadrature points (which should be contained in
             * #m_phys), this function calculates the \f$L_2\f$ error
             * of this approximation with respect to an exact
             * solution. The local distribution of the quadrature
             * points allows an elemental evaluation of this operation
             * through the functions StdRegions#StdExpansion#L2.
             *
             * The exact solution, also evaluated at the quadrature
             * points, should be contained in the variable #m_phys of
             * the ExpList object \a Sol.
             *
             * \param Sol An ExpList, containing the discrete
             * evaluation of the exact solution at the quadrature
             * points in its array #m_phys.
             *
             * \return The \f$L_2\f$ error of the approximation.
             */
            NekDouble L2 (const Array<OneD, const NekDouble> &soln);

            /**
             * \brief This function returns the number of elements in
             * the expansion.
             *
             * \return \f$N_{\mathrm{el}}\f$, the number of elements
             * in the expansion.
             */
            inline int GetExpSize(void)
            {
                return (*m_exp).size();
            }
      
            /**
             * \brief This function returns (a shared pointer to) the
             * local elemental expansion of the \f$n^{\mathrm{th}}\f$
             * element.
             *
             * \param n The index of the element concerned.
             *
             * \return (A shared pointer to) the local expansion of
             * the \f$n^{\mathrm{th}}\f$ element.
             */
            inline StdRegions::StdExpansionSharedPtr& GetExp(int n)
            {
                return (*m_exp)[n];
            }

            /**
             * \brief Get the start offset position for a global list
             * of #m_coeffs correspoinding to element n
             */
            inline int GetCoeff_Offset(int n)
            {
                return m_coeff_offset[n];
            }


            /**
             * \brief Get the start offset position for a global list
             * of m_phys correspoinding to  element  n
             */
            inline int GetPhys_Offset(int n)
            {
                return m_phys_offset[n];
            }
            
            /**
             * \brief This function returns (a reference to) the array 
             * \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs) containing all 
             * local expansion coefficients.
             *
             * If one wants to get hold of the underlying data without modifying them, 
             * rather use the function #GetCoeffs instead.
             *
             * \return (A reference to) the array #m_coeffs.
             */
            inline Array<OneD, NekDouble> &UpdateCoeffs()
            {
                m_transState = eLocal;
                return m_coeffs;
            }

            /**
             * \brief This function returns (a reference to) the array 
             * \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the function 
             * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the quadrature points.
             *
             * If one wants to get hold of the underlying data without modifying them, 
             * rather use the function #GetPhys instead.
             *
             * \return (A reference to) the array #m_phys.
             */
            inline Array<OneD, NekDouble> &UpdatePhys()
            {
                m_physState = true;
                return m_phys;
            }

            /**
             * \brief This function discretely evaluates the
             * derivative of a function \f$f(\boldsymbol{x})\f$ on the
             * domain consisting of all elements of the expansion.
             *
             * Given a function \f$f(\boldsymbol{x})\f$ evaluated at
             * the quadrature points, this function calculates the
             * derivatives \f$\frac{d}{dx_1}\f$, \f$\frac{d}{dx_2}\f$
             * and \f$\frac{d}{dx_3}\f$ of the function
             * \f$f(\boldsymbol{x})\f$ at the same quadrature
             * points. The local distribution of the quadrature points
             * allows an elemental evaluation of the derivative. This
             * is done by a call to the function
             * StdRegions#StdExpansion#PhysDeriv.
             *
             * \param inarray An array of size \f$Q_{\mathrm{tot}}\f$
             * containing the values of the function
             * \f$f(\boldsymbol{x})\f$ at the quadrature points
             * \f$\boldsymbol{x}_i\f$.  \param out_d0 The discrete
             * evaluation of the derivative\f$\frac{d}{dx_1}\f$ will
             * be stored in this array of size \f$Q_{\mathrm{tot}}\f$.
             * \param out_d1 The discrete evaluation of the
             * derivative\f$\frac{d}{dx_2}\f$ will be stored in this
             * array of size \f$Q_{\mathrm{tot}}\f$. Note that if no
             * memory is allocated for \a out_d1, the derivative
             * \f$\frac{d}{dx_2}\f$ will not be calculated.  \param
             * out_d2 The discrete evaluation of the
             * derivative\f$\frac{d}{dx_3}\f$ will be stored in this
             * array of size \f$Q_{\mathrm{tot}}\f$. Note that if no
             * memory is allocated for \a out_d2, the derivative
             * \f$\frac{d}{dx_3}\f$ will not be calculated.
             */
            void   PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &out_d0, 
                             Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                             Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);
            
            void PhysDeriv(const int dir, 
                           const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &out_d);


            // functions associated with DisContField
	    
	    const Array<OneD, const  boost::shared_ptr<ExpList1D> > &GetBndCondExpansions(void)
            {
	      return v_GetBndCondExpansions();
            }


            boost::shared_ptr<GenExpList1D> &GetTrace(void)
            {
                return v_GetTrace();
            }
	      
	    boost::shared_ptr<LocalToGlobalDGMap> &GetTraceMap(void) 
	    { 
		return v_GetTraceMap(); 
	    } 
	    

            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fx, 
                                          const Array<OneD, const NekDouble> &Fy, 
                                          Array<OneD, NekDouble> &outarray)
            {
                v_AddTraceIntegral(Fx,Fy,outarray);
            }

            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fn,   Array<OneD, NekDouble> &outarray)
            {
                v_AddTraceIntegral(Fn,outarray);
            }

            void AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd, 
                                          const Array<OneD, const NekDouble> &Bwd, 
                                          Array<OneD, NekDouble> &outarray)
            {
                v_AddTraceBiIntegral(Fwd,Bwd,outarray);
            }

            void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd)
            {
                v_GetFwdBwdTracePhys(Fwd,Bwd);
            }

            void GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field, 
                                    Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd)
            {
                v_GetFwdBwdTracePhys(field,Fwd,Bwd);
            }

            virtual void ExtractTracePhys(Array<OneD,NekDouble> &outarray)
            {
                v_ExtractTracePhys(outarray);
            }


            void ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                v_ExtractTracePhys(inarray,outarray);
            }

            inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& GetBndConditions()
            {
                return v_GetBndConditions();       
            }

            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                v_EvaluateBoundaryConditions(time);
            }

            
            // Routines for continous matrix solution 
            /**
             * \brief This function calculates the result of the
             * multiplication of a matrix of type specified by \a mkey
             * with a vector given by \a inarray.
             *
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
             * where \f$\boldsymbol{M}^e\f$ are the local matrices of type specified by the 
             * key \a mkey. The decoupling of the local matrices allows for a local 
             * evaluation of the operation. However, rather than a local matrix-vector 
             * multiplication, the local operations are evaluated as implemented in the 
             * function StdRegions#StdExpansion#GeneralMatrixOp.
             *
             * \param mkey This key uniquely defines the type matrix required for 
             * the operation.
             * \param inarray The vector \f$\boldsymbol{\hat{u}}_l\f$ of size 
             * \f$N_{\mathrm{eof}}\f$.
             * \param outarray The resulting vector of size \f$N_{\mathrm{eof}}\f$.
             */      
            void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                 const Array<OneD,const NekDouble> &inarray, 
                                       Array<OneD,      NekDouble> &outarray,
                                 bool  UseContCoeffs = false)
            {
                v_GeneralMatrixOp(gkey,inarray,outarray,UseContCoeffs);
            }

            void GeneralMatrixOp_IterPerExp(const GlobalMatrixKey             &gkey,
                                            const Array<OneD,const NekDouble> &inarray, 
                                            Array<OneD,      NekDouble> &outarray);


            std::vector<SpatialDomains::FieldDefinitionsSharedPtr> GetFieldDefinitions(void);
            // \brief Append the element data listed in elements
            // fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata);

            /// \brief Extract the data in fielddata into the m_coeff list 
            void ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field);

            // load global optimisation parameters
            inline void ReadGlobalOptimizationParameters(const std::string &infilename)
            {
                m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>::AllocateSharedPtr(infilename);
            }

        protected:
            
            boost::shared_ptr<DNekMat> GenGlobalMatrixFull(const GlobalLinSysKey &mkey, const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            void SetCoeffPhys(void);
                

            /**
             * \brief The total number of local degrees of freedom. 
             * 
             * #m_ncoeffs \f$=N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_l\f$
             */
            int m_ncoeffs; 
            
            /** \brief The total number of quadrature points. 
             *
             * #m_npoints \f$=Q_{\mathrm{tot}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_Q\f$
             */
            int m_npoints; 
            
            /** 
             * \brief Concatenation of all local expansion coefficients.
             *
             * The array of length #m_ncoeffs\f$=N_{\mathrm{eof}}\f$ which is the 
             * concatenation of the local expansion coefficients \f$\hat{u}_n^e\f$ over 
             * all \f$N_{\mathrm{el}}\f$ elements 
             * \f[\mathrm{\texttt{m\_coeffs}}=\boldsymbol{\hat{u}}_{l} = 
             * \underline{\boldsymbol{\hat{u}}}^e = \left [ \begin{array}{c}
             * \boldsymbol{\hat{u}}^{1} \       \
             * \boldsymbol{\hat{u}}^{2} \       \
             * \vdots \                                                 \
             * \boldsymbol{\hat{u}}^{{{N_{\mathrm{el}}}}} \end{array} \right ],\quad
             * \mathrm{where}\quad \boldsymbol{\hat{u}}^{e}[n]=\hat{u}_n^{e}\f]
             */
            Array<OneD, NekDouble> m_coeffs;
            
            /**
             * \brief The global expansion evaluated at the quadrature points
             *
             * The array of length #m_npoints\f$=Q_{\mathrm{tot}}\f$ containing the evaluation 
             * of \f$u^{\delta}(\boldsymbol{x})\f$ at the quadrature points 
             * \f$\boldsymbol{x}_i\f$. 
             * \f[\mathrm{\texttt{m\_phys}}=\boldsymbol{u}_{l} = 
             * \underline{\boldsymbol{u}}^e = \left [ \begin{array}{c}
             * \boldsymbol{u}^{1} \             \
             * \boldsymbol{u}^{2} \             \
             * \vdots \                                                 \
             * \boldsymbol{u}^{{{N_{\mathrm{el}}}}} \end{array} \right ],\quad
             * \mathrm{where}\quad \boldsymbol{u}^{e}[i]=u^{\delta}(\boldsymbol{x}_i)\f]
             */
            Array<OneD, NekDouble> m_phys; 
            
            /**
             * \brief The transformed state of the array of coefficients of the expansion.
             *
             * where #TransState is the enumeration which can attain the following values:
             * - <em>eNotSet</em>: The coefficients are not set.
             * - <em>eLocal</em>: The array #m_coeffs is filled with the proper local coefficients.
             * - <em>eContinuous</em>: The array #m_contCoeffs is filled with the proper 
             *   global coefficients.
             * - <em>eLocalCont</em>: Both the arrays #m_coeffs and #m_contCoeffs are filled with 
             *   the proper coefficients.
             */
            TransState m_transState; 
            
            /**
             * \brief The state of the array #m_phys.
             *
             * Indicates whether the array #m_phys, created to contain the evaluation of 
             * \f$u^{\delta}(\boldsymbol{x})\f$ at the quadrature points 
             * \f$\boldsymbol{x}_i\f$, is filled with these values.
             */
            bool       m_physState;  
            
            /**
             * \brief The list of local expansions.
             *
             * The (shared pointer to the) vector containing (shared pointers to) all local
             * expansions. The fact that the local expansions are all stored as a (pointer to a) 
             * #StdExpansion, the abstract base class for all local expansions, allows a general 
             * implementation where most of the routines for the derived classes are defined in 
             * the #ExpList base class.
             */
            boost::shared_ptr<StdRegions::StdExpansionVector> m_exp; 
            
            /**
             * \brief Offset of elemental data into the array #m_coeffs
             */
            Array<OneD, int>  m_coeff_offset;

            /**
             * \brief Offset of elemental data into the array #m_phys
             */
            Array<OneD, int>  m_phys_offset;

            NekOptimize::GlobalOptParamSharedPtr m_globalOptParam;

            BlockMatrixMapShPtr  m_blockMat;

            /**
             * \brief This function assembles the block diagonal
             * matrix of local matrices of the type \a mtype.
             *
             * This function assembles the block diagonal matrix
             * \f$\underline{\boldsymbol{M}}^e\f$, which is the
             * concatenation of the local matrices
             * \f$\boldsymbol{M}^e\f$ of the type \a mtype, that is
             *
             * \f[
             * \underline{\boldsymbol{M}}^e = \left[
             * \begin{array}{cccc}
             * \boldsymbol{M}^1 & 0 & \hspace{3mm}0 \hspace{3mm}& 0 \\
             *  0 & \boldsymbol{M}^2 & 0 & 0 \\
             *  0 &  0 & \ddots &  0 \\
             *  0 &  0 & 0 & \boldsymbol{M}^{N_{\mathrm{el}}} \end{array} \right].\f]
             *
             * \param mtype the type of matrix to be assembled
             * \param scalar an optional parameter 
             * \param constant an optional parameter 
             */
            const DNekScalBlkMatSharedPtr  GenBlockMatrix(const GlobalMatrixKey &gkey);
            const DNekScalBlkMatSharedPtr& GetBlockMatrix(const GlobalMatrixKey &gkey);
            void MultiplyByBlockMatrix(const GlobalMatrixKey             &gkey,
                                       const Array<OneD,const NekDouble> &inarray, 
                                             Array<OneD,      NekDouble> &outarray);
            
            boost::shared_ptr<GlobalMatrix>  GenGlobalMatrix(const GlobalMatrixKey &mkey, 
                                                             const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            
            /**
             * \brief This operation constructs the global linear system of type \a mkey.
             *
             * Consider a linear system \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{\hat{f}}\f$ 
             * to be solved, where \f$\boldsymbol{M}\f$ is a matrix of type specified by 
             * \a mkey. This function assembles the global system matrix 
             * \f$\boldsymbol{M}\f$ out of the elemental submatrices \f$\boldsymbol{M}^e\f$.
             * This is equivalent to:
             * \f[ \boldsymbol{M}=\boldsymbol{\mathcal{A}}^T
             * \underline{\boldsymbol{M}}^e\boldsymbol{\mathcal{A}}.\f]
             * where the matrix \f$\boldsymbol{\mathcal{A}}\f$ is a sparse permutation 
             * matrix of size \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$. However, due 
             * to the size and sparsity of the matrix \f$\boldsymbol{\mathcal{A}}\f$, it is 
             * more efficient to assemble the global matrix using the mapping array 
             * \a map\f$[e][i]\f$ contained in the input argument \a locToGloMap. 
             * The global assembly is then evaluated as:
             * \f[ \boldsymbol{M}\left[\mathrm{\texttt{map}}[e][i]\right]
             * \left[\mathrm{\texttt{map}}[e][j]\right]=\mathrm{\texttt{sign}}[e][i]\cdot 
             * \mathrm{\texttt{sign}}[e][j] \cdot\boldsymbol{M}^e[i][j]\f]
             * where the values \a sign\f$[e][i]\f$ ensure the correct connectivity.
             *
             * \param mkey A key which uniquely defines the global matrix to be constructed.
             * \param locToGloMap Contains the mapping array and required information for 
             * the transformation from local to global degrees of freedom.
             * \return (A shared pointer to) the global linear system formed by the global 
             * matrix \f$\boldsymbol{M}\f$.
             */
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSysFullDirect(const GlobalLinSysKey &mkey, const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            
            void GlobalEigenSystem(const boost::shared_ptr<DNekMat> &Gmat,
                                   Array<OneD, NekDouble> &EigValsReal, 
                                   Array<OneD, NekDouble> &EigValsImag, 
                                   Array<OneD, NekDouble> &EigVecs = NullNekDouble1DArray);
            
            /**
             * \brief This function constructs the necessary global matrices required for 
             * solving the linear system of type \a mkey by static condensation.
             *
             * Consider the linear system \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{\hat{f}}\f$.
             * Distinguishing between the boundary and interior components of 
             * \f$\boldsymbol{\hat{u}}_g\f$ and \f$\boldsymbol{\hat{f}}\f$ using 
             * \f$\boldsymbol{\hat{u}}_b\f$,\f$\boldsymbol{\hat{u}}_i\f$ and 
             * \f$\boldsymbol{\hat{f}}_b\f$,\f$\boldsymbol{\hat{f}}_i\f$ respectively, 
             * this system can be split into its constituent parts as
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
             * \f$\boldsymbol{M}\f$ resulting from boundary-boundary mode interactions, 
             * \f$\boldsymbol{M}_{c1}\f$ and \f$\boldsymbol{M}_{c2}\f$ represent the 
             * components resulting from coupling between the boundary-interior modes, 
             * and \f$\boldsymbol{M}_i\f$ represents the components of \f$\boldsymbol{M}\f$ 
             * resulting from interior-interior mode interactions. 
             * 
             * The solution of the linear system can now be determined in two steps:
             * \f{eqnarray*}
             * \mathrm{step 1:}&\quad&(\boldsymbol{M}_b-\boldsymbol{M}_{c1}
             * \boldsymbol{M}_i^{-1}\boldsymbol{M}_{c2}) \boldsymbol{\hat{u}_b} = 
             * \boldsymbol{\hat{f}}_b - \boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
             * \boldsymbol{\hat{f}}_i,\nonumber \\
             * \mathrm{step 2:}&\quad&\boldsymbol{\hat{u}_i}=\boldsymbol{M}_i^{-1}
             * \left( \boldsymbol{\hat{f}}_i  - \boldsymbol{M}_{c2}\boldsymbol{\hat{u}_b}
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
             * where \f$\boldsymbol{\mathcal{A}}_b \f$ is the permutation matrix which 
             * scatters from global to local degrees of freedom, only the following four 
             * matrices should be constructed:
             * - \f$\underline{[\boldsymbol{M}^e_i]}^{-1}\f$
             * - \f$\underline{\boldsymbol{M}^e_{c1}}\underline{[\boldsymbol{M}^e_i]}^{-1}\f$
             * - \f$\underline{\boldsymbol{M}^e_{c2}}\f$
             * - The Schur complement: \f$\boldsymbol{M}_{\mathrm{Schur}}=
             *   \quad\boldsymbol{M}_b-\boldsymbol{M}_{c1}\boldsymbol{M}_i^{-1}
             *   \boldsymbol{M}_{c2}\f$
             *
             * The first three matrices are just a concatenation of the corresponding local 
             * matrices and they can be created as such. They also allow for an elemental 
             * evaluation of the operations concerned.
             *
             *
             * The global Schur complement however should be assembled from the 
             * concatenation of the local elemental Schur complements, that is,
             * \f[ \boldsymbol{M}_{\mathrm{Schur}}=\boldsymbol{M}_b - \boldsymbol{M}_{c1} 
             * \boldsymbol{M}_i^{-1} \boldsymbol{M}_{c2} =
             * \boldsymbol{\mathcal{A}}_b^T \left [\  \underline{\boldsymbol{M}^e_b  - 
             * \boldsymbol{M}^e_{c1} [\boldsymbol{M}^e_i]^{-1} (\boldsymbol{M}^e_{c2})}\ 
             * \right ] \boldsymbol{\mathcal{A}}_b \f]
             * and it is the only matrix operation that need to be evaluated on a global 
             * level when using static condensation.
             * However, due to the size and sparsity of the matrix 
             * \f$\boldsymbol{\mathcal{A}}_b\f$, it is more efficient to assemble the 
             * global Schur matrix using the mapping array bmap\f$[e][i]\f$ contained in 
             * the input argument \a locToGloMap. The global Schur complement is then 
             * constructed as:
             * \f[\boldsymbol{M}_{\mathrm{Schur}}\left[\mathrm{\a bmap}[e][i]\right]
             * \left[\mathrm{\a bmap}[e][j]\right]=\mathrm{\a bsign}[e][i]\cdot 
             * \mathrm{\a bsign}[e][j] \cdot\boldsymbol{M}^e_{\mathrm{Schur}}[i][j]\f]
             * All four matrices are stored in the \a GlobalLinSys returned by this function.
             *
             * \param mkey A key which uniquely defines the global matrix to be constructed.
             * \param locToGloMap Contains the mapping array and required information for 
             * the transformation from local to global degrees of freedom.
             * \return (A shared pointer to) the statically condensed global linear system.
             */
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSysStaticCond(const GlobalLinSysKey &mkey, 
                                                                       const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /**
             * \brief This operation constructs the global linear system of type \a mkey.
             * 
             * Consider a linear system \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{f}\f$ to be 
             * solved. Dependent on the solution method, this function constructs
             * - <b>The full linear system</b><BR>
             *   A call to the function #GenGlobalLinSysFullDirect
             * - <b>The statically condensed linear system</b><BR>
             *   A call to the function #GenGlobalLinSysStaticCond
             *
             * \param mkey A key which uniquely defines the global matrix to be constructed.
             * \param locToGloMap Contains the mapping array and required information for
             *  the transformation from local to global degrees of freedom.
             * \return (A shared pointer to) the global linear system in required format.
             */
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSys(const GlobalLinSysKey &mkey, 
                                                             const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);
            
            /**
             * \brief Generate a GlobalLinSys from information
             * provided by the key "mkey" and the mapping provided in
             * LocToGloBaseMap
             */      
            boost::shared_ptr<GlobalLinSys> GenGlobalBndLinSys(const GlobalLinSysKey     &mkey, const LocalToGlobalBaseMap &LocToGloBaseMap);


            // functions associated with DisContField
	    inline virtual const Array<OneD,const boost::shared_ptr<ExpList1D> > &v_GetBndCondExpansions(void)
            {
	      ASSERTL0(false,"This method is not defined or valid for this class type");                
          static Array<OneD,const boost::shared_ptr<ExpList1D> > result;
          return result;
            }

            inline virtual boost::shared_ptr<GenExpList1D> &v_GetTrace(void)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");
                static boost::shared_ptr<GenExpList1D> returnVal;
                return returnVal;
            }
            
	    inline virtual boost::shared_ptr<LocalToGlobalDGMap> &v_GetTraceMap(void) 
	    { 
		ASSERTL0(false,"This method is not defined or valid for this class type"); 
                static boost::shared_ptr<LocalToGlobalDGMap> result;
                return result;
	    } 

            virtual void v_AddTraceIntegral(const Array<OneD, const NekDouble> &Fx, 
                                          const Array<OneD, const NekDouble> &Fy, 
                                          Array<OneD, NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_AddTraceIntegral(const Array<OneD, const NekDouble> &Fn,   Array<OneD, NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd, 
                                          const Array<OneD, const NekDouble> &Bwd, 
                                          Array<OneD, NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                              Array<OneD,NekDouble> &Bwd)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field, 
                                              Array<OneD,NekDouble> &Fwd, 
                                              Array<OneD,NekDouble> &Bwd)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_ExtractTracePhys(Array<OneD,NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }

            virtual void v_MultiplyByInvMassMatrix(const Array<OneD,const NekDouble> &inarray, 
                                                         Array<OneD,      NekDouble> &outarray, 
                                                   bool  UseContCoeffs)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
            }


            virtual void v_HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,       NekDouble> &outarray,
                                     NekDouble lambda,
                                     bool      UseContCoeffs,
                                     const Array<OneD, const NekDouble>& dirForcing)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type"); 
            }
	  
            // wrapper functions about virtual functions
            virtual  const Array<OneD, const NekDouble> &v_GetContCoeffs() const 
            {
                ASSERTL0(false,"This method is not valid in this class type");
                return NullNekDouble1DArray;
            }

            virtual void v_BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                          Array<OneD,       NekDouble> &outarray, 
                                    bool  UseContCoeffs)
            {
                BwdTrans_IterPerExp(inarray,outarray);
            }
            
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                          Array<OneD,       NekDouble> &outarray,
                                    bool  UseContCoeffs)
            {
                FwdTrans_IterPerExp(inarray,outarray);
            }
            
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                                 Array<OneD,       NekDouble> &outarray,
                                           bool  UseContCoeffs)
            {
                IProductWRTBase_IterPerExp(inarray,outarray);
            }

            virtual void v_GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                           const Array<OneD,const NekDouble> &inarray, 
                                                 Array<OneD,      NekDouble> &outarray,
                                           bool  UseContCoeffs)
            {
                GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
            }

        private:

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");                
                static Array<OneD,const SpatialDomains::BoundaryConditionShPtr> result;
                return result;
            }

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");            
            }

            
        };

        typedef boost::shared_ptr<ExpList>      ExpListSharedPtr;

    static ExpList NullExpList;
    
  } //end of namespace
} //end of namespace

#endif // EXPLIST_H

/**
* $Log: ExpList.h,v $
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
