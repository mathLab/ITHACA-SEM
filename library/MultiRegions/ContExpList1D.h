///////////////////////////////////////////////////////////////////////////////
//
// File ContExpList1D.cpp
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
// Description: Continusou Expansion list definition in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_CONTEXPLIST1D_H
#define NEKTAR_LIB_MULTIREGIONS_CONTEXPLIST1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/LocalToGlobalMap1D.h>
#include <MultiRegions/GlobalLinSys.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * \brief This class is the abstraction of a global \f$C^0\f$ continuous 
         * one-dimensional spectral/hp element expansion.
         *
         * All local elemental expansions are now connected to form a global spectral/hp 
         * representation. In this case, by imposing \f$C^0\f$ continuity across the 
         * element interfaces, the expansion is chosen to be globally continuous. This type 
         * of global expansion can be defined as:
         * \f[u^{\delta}(x_i)=\sum_{n=0}^{N_{\mathrm{dof}}-1}\hat{u}_n
         * \Phi_n(x_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}\sum_{n=0}^{N^{e}_m-1}
         * \hat{u}_n^e\phi_n^e(x_i) \f]
         * where \f$N_{\mathrm{dof}}\f$ refers to the number of global modes 
         * \f$\Phi_n(\boldsymbol{x})\f$. The class #ContExpList1D is derived from the 
         * classes #ExpList and #ExpList1D. In addition to the member variables defined in 
         * these classes, the class #ContExpList1D also contains the member data 
         * #m_contNcoeffs and #m_contCoeffs, both related to the global coefficients 
         * \f$\hat{u}_n\f$. Furthermore, a mapping array #m_locToGloMap relating the global 
         * degrees of freedom \f$\hat{u}_n\f$ and local degrees of freedom 
         * \f$\hat{u}_n^e\f$ is provided within this class.
         */  
    class ContExpList1D: 
        public ExpList1D 
        {
        public:
            /**
             * \brief The default constructor.
             */  
            ContExpList1D();
             
            /**
             * \brief 
             */  
            ContExpList1D(const LibUtilities::BasisKey &Ba, 
                          const SpatialDomains::MeshGraph1D &graph1D,
                          const bool constructMap = true);
 
            /**
             * \brief This constructor sets up a global continuous expansions based on an 
             * input mesh.
             *
             * Given a mesh \a graph1D, containing information about the domain and the 
             * spectral/hp element expansion, this constructor fills the list of local 
             * expansions #m_exp with the proper expansions, calculates the total number of 
             * quadrature points \f$x_i\f$ and local expansion coefficients 
             * \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs and #m_phys. 
             * Furthermore, if required, it constructs the mapping array (contained in
             * #m_locToGloMap) for the transformation between local elemental level and 
             * global level, it calculates the total number global expansion coefficients 
             * \f$\hat{u}_n\f$ and allocates memory for the array #m_contCoeffs.
             *
             * \param graph1D A mesh, containing information about the domain and the 
             * spectral/hp element expansion.
             * \param constructMap An optional parameter to indicate whether the mapping 
             * array should be constructed or not. 
             */  
            ContExpList1D(SpatialDomains::MeshGraph1D &graph1D,
                          const bool constructMap = true);
             
            /**
             * \brief The copy constructor.
             */  
            ContExpList1D(const ContExpList1D &In);
             
            /**
             * \brief The default destructor.
             */  
            ~ContExpList1D();

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
            inline void ContToLocal()
            {
                m_locToGloMap->ContToLocal(m_contCoeffs,m_coeffs);
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
            inline void ContToLocal(const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->ContToLocal(inarray,outarray);
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
            inline void LocalToCont()
            {
                m_locToGloMap->LocalToCont(m_coeffs,m_contCoeffs);
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
            inline void Assemble(const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->Assemble(inarray,outarray);
            }
 
            /**
             * \brief This function returns the map from local to global level.
             */  
            inline const LocalToGlobalMapSharedPtr& GetLocalToGlobalMap() const
            {
                return  m_locToGloMap;
            }

            /**
             * \brief This function calculates the inner product of a function 
             * \f$f(x)\f$ with respect to all <em>global</em> expansion modes 
             * \f$\phi_n^e(x)\f$.
             * 
             * The operation is evaluated locally (i.e. with respect to all local expansion 
             * modes) by the function ExpList#IProductWRTBase. The inner product with 
             * respect to the global expansion modes is than obtained by a global assembly 
             * operation.
             *
             * The values of the function \f$f(x)\f$ evaluated at the 
             * quadrature points \f$x_i\f$ should be contained in the variable 
             * #m_phys of the ExpList object \a in. The result is stored in the array 
             * #m_contCoeffs.
             *
             * \param In An ExpList, containing the discrete evaluation of 
             * \f$f(x)\f$ at the quadrature points in its array #m_phys.
             */  
            void IProductWRTBase(const ExpList &In);
       
            /**
             * \brief This function performs the global forward transformation of a 
             * function \f$f(x)\f$.
             * 
             * Given a function \f$f(x)\f$ defined at the quadrature points, 
             * this function determines the global coefficients \f$\hat{u}_n\f$ employing a 
             * discrete elemental Galerkin projection from physical space to coefficient 
             * space. The operation is evaluated as
             * \f[\boldsymbol{\hat{u}}_g=\boldsymbol{M}^{-1}\boldsymbol{\hat{f}}\f]
             * where \f$\boldsymbol{\hat{f}}\f$ is the array of entries 
             * \f$\boldsymbol{\hat{f}}[n]=\int\Phi_n(x)f(x)dx\f$ 
             * and \f$\boldsymbol{M}\f$ is the global mass matrix.
             *
             * The values of the function \f$f(x)\f$ evaluated at the 
             * quadrature points \f$x_i\f$ should be contained in the 
             * variable #m_phys of the ExpList object \a In. The resulting global 
             * coefficients \f$\hat{u}_g\f$ are stored in the array #m_contCoeffs.
             *
             * \param In An ExpList, containing the discrete evaluation of 
             * \f$u(x)\f$ at the quadrature points in its array #m_phys.
             */ 
            void FwdTrans(const ExpList &In);
       
            /**
             * \brief This function performs the backward transformation of the spectral/hp 
             * element expansion.
             *
             * Given the coefficients of an expansion, this function evaluates the 
             * spectral/hp expansion \f$u^{\delta}(x)\f$ at the quadrature 
             * points \f$x_i\f$. This operation is evaluated locally by the 
             * function ExpList#BwdTrans.
             *
             * The coefficients of the expansion should be contained in the variable 
             * #m_coeffs of the ExpList object \a In. The resulting physical values at the 
             * quadrature points \f$u^{\delta}(x_i)\f$ are stored in the array 
             * #m_phys.
             *
             * \param In An ExpList, containing the local coefficients \f$\hat{u}_n^e\f$ 
             * in its array #m_coeffs.
             */  
            void BwdTrans(const ExpList &In);
 
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
            void GeneralMatrixOp(const GlobalLinSysKey            &gkey,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble>           &outarray);
	    
	protected:
            /**
             * \brief (A shared pointer to) the object which contains all the required 
             * information for the transformation from local to global degrees of freedom.
             */  
	    LocalToGlobalMapSharedPtr m_locToGloMap;
 
            /**
             * \brief The total number of global degrees of freedom. 
             * #m_contNcoeffs\f$=N_{\mathrm{dof}}\f$
             */  
      	    int                       m_contNcoeffs;
 
            /**
             * \brief The array of length #m_ncoeffs\f$=N_{\mathrm{dof}}\f$ containing the 
             * global expansion coefficients. 
             */  
	    Array<OneD, NekDouble>    m_contCoeffs;
 
            /**
             * \brief (A shared pointer to) a list which collects all the global matrices 
             * being assembled, such that they should be constructed only once.
             */  	
            GlobalLinSysMapShPtr      m_globalMat;
            
        private:

        };


        typedef boost::shared_ptr<ContExpList1D>      ContExpList1DSharedPtr;
        typedef std::vector<ContExpList1DSharedPtr>   ContExpList1DVector;
        typedef std::vector<ContExpList1DSharedPtr>::iterator ContExpList1DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

/**
* $Log: ContExpList1D.h,v $
* Revision 1.31  2008/05/29 21:35:03  pvos
* Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
*
* Revision 1.30  2008/04/06 06:00:06  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.29  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.28  2007/12/17 13:05:04  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.27  2007/12/06 22:52:29  pvos
* 2D Helmholtz solver updates
*
* Revision 1.26  2007/11/07 20:29:52  jfrazier
* Modified to use new expansion list contained in meshgraph.
*
* Revision 1.25  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.24  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.23  2007/09/25 14:25:29  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.22  2007/09/03 19:58:30  jfrazier
* Formatting.
*
* Revision 1.21  2007/07/26 00:07:50  bnelson
* Fixed linux compiler errors.
*
* Revision 1.20  2007/07/23 16:06:30  sherwin
* Put a std::map to hold global matrix systems
*
* Revision 1.19  2007/07/22 23:04:19  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.18  2007/07/20 02:04:10  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.17  2007/07/19 20:02:25  sherwin
* Generalised global matrix solver
*
* Revision 1.16  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.15  2007/07/13 15:22:12  sherwin
* Update for Helmholtz (working without bcs )
*
* Revision 1.14  2007/07/13 09:02:23  sherwin
* Mods for Helmholtz solver
*
* Revision 1.13  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.12  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.11  2007/06/08 12:58:26  sherwin
* Added ContField1D and remove previous structure using Fields
*
* Revision 1.10  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
