///////////////////////////////////////////////////////////////////////////////
//
// File Stdexpansion.h
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
// Description: Class definition StdExpansion which is the base class
// to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STANDARDEXPANSION_H
#define NEKTAR_LIB_STDREGIONS_STANDARDEXPANSION_H

#include <fstream>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/SpatialDomainsDeclarations.hpp>
#include <StdRegions/LocalRegionsDeclarations.hpp>

#include <StdRegions/StdExpMap.h>
#include <StdRegions/StdMatrixKey.h>
#include <StdRegions/StdLinSysKey.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        class StdSegExp;

        /// This is a macro which is used to encapsulating the points
        /// properties at which a given expansion basis is evaluated.
        /// The i index refers to the ith basis definition of a given expansion.
        ///
        /// This should/could be a typdef expression when such a thing
        /// comes into the convention.
#define ExpPointsProperties(i) LibUtilities::PointsManager()[m_base[i]->GetPointsKey()] 

        /** \brief The base class for all shapes
        *   
        *  This is the lowest level basic class for all shapes and so
        *  contains the definition of common data and common routine to all
        *  elements
        */
        class StdExpansion
        {

        public:

            /** \brief Default Constructor */
            StdExpansion();

            /** \brief Constructor */
            StdExpansion(const int numbases, const LibUtilities::BasisKey &Ba, 
                const LibUtilities::BasisKey &Bb,
                const LibUtilities::BasisKey &Bc, int numcoeffs);

            /** \brief Copy Constructor */
            StdExpansion(const StdExpansion &T);

            /** \brief Destructor */
            virtual ~StdExpansion();

            // Standard Expansion Routines Applicable Regardless of Region

            /** \brief This function returns the number of 1D bases used in 
            *  the expansion 
            *
            *  \return returns the number of 1D bases used in the expansion, 
            *  which is equal to number dimension of the expansion
            */
            inline int GetNumBases() const
            {
                return m_numbases;
            }

            /** \brief This function gets the shared point to basis 
            *  
            *  \return returns the shared pointer to the bases  
            */
            inline const boost::shared_array<LibUtilities::BasisSharedPtr> GetBase() const
            {
                return(m_base);
            }

            /** \brief This function gets the shared point to basis in
            *  the \a dir direction
            *  
            *  \return returns the shared pointer to the basis in
            *  directin \a dir
            */
            inline const LibUtilities::BasisSharedPtr GetBasis(int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(m_base[dir]);
            }



            /** \brief This function returns the total number of coefficients 
            *  used in the expansion 
            *  
            *  \return returns the total number of coefficients (which is 
            *  equivalent to the total number of modes) used in the expansion
            */
            inline int GetNcoeffs(void)
            {
                return(m_ncoeffs);
            }

            /** \brief This function returns a pointer to the coefficient array
            *  \f$ \mathbf{\hat{u}}\f$ 
            *
            *  The coefficient array \f$ \mathbf{\hat{u}}\f$ corresponds to the 
            *  class attribute #m_coeffs (which is in coefficient space)
            *
            *  \return returns a pointer to the coefficient array 
            *  \f$ \mathbf{\hat{u}}\f$ 
            */
            inline BstShrDArray GetCoeffs(void)
            {
                return(m_coeffs);
            }

            /** \brief This function sets the coefficient array 
            *  \f$ \mathbf{\hat{u}}\f$ (implemented as the class attribute 
            *  #m_coeffs) to the values given by \a coeffs
            *
            *  Using this function actually determines the entire expansion
            *
            *  \param coeffs the array of values to which #m_coeffs should 
            *  be set
            */
            inline void SetCoeffs(double *coeffs)
            {
                Vmath::Vcopy(m_ncoeffs, coeffs, 1, &m_coeffs[0], 1);
            }

            /** \brief This function sets the i th coefficient  
            *  \f$ \mathbf{\hat{u}}[i]\f$ to the value given by \a coeff
            *
            *  #m_coeffs[i] will be set to the  value given by \a coeff
            *
            *  \param i the index of the coefficient to be set
            *  \param coeff the value of the coefficient to be set
            */
            inline void SetCoeff(const int i, double coeff)
            {
                m_coeffs[i] = coeff;
            }

            /** \brief This function returns a pointer to the array
            *  \f$\mathbf{u}\f$ (which is in physical space)
            *
            *  The array \f$ \mathbf{u}\f$ corresponds to the 
            *  class attribute #m_phys and contains the values of a function
            *  evaluates at the quadrature points, 
            *  i.e. \f$\mathbf{u}[m]=u(\mathbf{\xi}_m)\f$
            *
            *  \return returns a pointer to the array \f$\mathbf{u}\f$ 
            */
            inline BstShrDArray  GetPhys(void)
            {
                return(m_phys);
            }

            /** \brief This function returns the total number of quadrature
            *  points used in the element 
            *  
            *  \return returns the total number of quadrature points
            */
            inline int GetTotPoints()
            {
                int i;
                int nqtot = 1;

                for(i=0; i<m_numbases; ++i)
                {
                    nqtot *= ExpPointsProperties(i)->GetNumPoints();
                }

                return  nqtot;
            }

            /** \brief This function sets the array 
            *  \f$ \mathbf{u}\f$ (implemented as the class attribute 
            *  #m_phys) to the values given by \a phys
            *
            *  Using this function corresponds to storing a function \f$u\f$
            *  (evaluated at the quadrature points) in the class attribute
            *  #m_phys 
            *
            *  \param phys the array of values to which #m_phys should be set
            */
            inline void  SetPhys(const double *phys)
            {
                int nqtot = GetTotPoints();

                Vmath::Vcopy(nqtot, phys, 1, &m_phys[0], 1);
            }

            /** \brief This function returns the type of basis used in the \a dir
            *  direction
            *  
            *  The different types of bases implemented in the code are defined 
            *  in the LibUtilities::BasisType enumeration list. As a result, the
            *  funcion will return one of the types of this enumeration list.
            *  
            *  \param dir the direction
            *  \return returns the type of basis used in the \a dir direction
            */
            inline LibUtilities::BasisType GetBasisType(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(m_base[dir]->GetBasisType());
            }

            /** \brief This function returns the number of expansion modes 
            *  in the \a dir direction
            *  
            *  \param dir the direction
            *  \return returns the number of expansion modes in the \a dir 
            *  direction
            */
            inline int GetBasisNumModes(const int dir) const
            {
                ASSERTL1(dir < m_numbases,"dir is larger than m_numbases");
                return(m_base[dir]->GetNumModes());
            }

            /** \brief This function returns the type of quadrature points used 
            *  in the \a dir direction
            *  
            *  The different types of quadrature points implemented in the code
            *  are defined in the LibUtilities::PointsType enumeration list. 
            *  As a result, the funcion will return one of the types of this 
            *  enumeration list.
            *  
            *  \param dir the direction
            *  \return returns the type of quadrature points  used in the \a dir
            *  direction
            */
            inline LibUtilities::PointsType GetPointsType(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(m_base[dir]->GetPointsType());
            }

            /** \brief This function returns the number of quadrature points 
            *  in the \a dir direction
            *  
            *  \param dir the direction
            *  \return returns the number of quadrature points in the \a dir 
            *  direction
            */
            inline int GetNumPoints(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(ExpPointsProperties(dir)->GetNumPoints());
            }

            /** \brief This function returns a pointer to the array containing
            *  the quadrature points in \a dir direction
            *
            *  \param dir the direction
            *  \return returns a pointer to the array containing
            *  the quadrature points in \a dir direction 
            */
            inline double *GetPoints(const int dir) const
            {
                return ExpPointsProperties(dir)->GetZ();
            }


            double operator[] (const int i) const
            {
                ASSERTL1((i >= 0) && (i < m_ncoeffs),
                    "Invalid Index used in [] operator");
                return m_coeffs[i];
            }

            double& operator[](const int i)
            {
                ASSERTL1((i >= 0) && (i < m_ncoeffs),
                    "Invalid Index used in [] operator");
                return m_coeffs[i];
            }

            // Wrappers around virtual Functions

            /** \brief This function returns the number of vertices of the 
            *  expansion domain
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetNverts() 
            * 
            *  \return returns the number of vertices of the expansion domain
            */
            int GetNverts()
            {
                return v_GetNverts();
            }

            /** \brief This function returns the number of edges of the 
            *  expansion domain
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetNedges() 
            * 
            *  \return returns the number of edges of the expansion domain
            */
            int GetNedges()
            {
                return v_GetNedges();
            }

            /** \brief This function returns the number of expansion coefficients
            *  belonging to the \a i-th edge  
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetEdgeNcoeffs() 
            * 
            *  \param i specifies which edge
            *  \return returns the number of expansion coefficients belonging to
            *  the \a i-th edge
            */
            int GetEdgeNcoeffs(const int i)
            {
                return v_GetEdgeNcoeffs(i);
            }

            /** \brief This function returns the type of expansion basis on the
            *  \a i-th edge  
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetEdgeBasisType() 
            *
            *  The different types of bases implemented in the code are defined 
            *  in the LibUtilities::BasisType enumeration list. As a result, the
            *  funcion will return one of the types of this enumeration list.
            * 
            *  \param i specifies which edge
            *  \return returns the expansion basis on the \a i-th edge
            */
            LibUtilities::BasisType GetEdgeBasisType(const int i)
            {
                return v_GetEdgeBasisType(i);
            }

            /** \brief This function returns the number of faces of the 
            *  expansion domain
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetNFaces() 
            * 
            *  \return returns the number of faces of the expansion domain
            */
            int GetNfaces()
            {
                return v_GetNfaces();
            }

            /** \brief This function returns the shape of the expansion domain 
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_DetShapeType() 
            *
            *  The different shape types implemented in the code are defined 
            *  in the ::ShapeType enumeration list. As a result, the
            *  funcion will return one of the types of this enumeration list.
            *  
            *  \return returns the shape of the expansion domain
            */	    
            ShapeType DetShapeType()
            {
                return v_DetShapeType();
            }

            /** \brief This function performs the Backward transformation from 
            *  coefficient space to physical space
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_BwdTrans() 
            *
            *  Based on the expansion coefficients, this function evaluates the
            *  expansion at the quadrature points. This is equivalent to the 
            *  operation \f[ u(\xi_{1i}) =
            *  \sum_{p=0}^{P-1} \hat{u}_p \phi_p(\xi_{1i}) \f] which can be 
            *  evaluated as \f$ {\bf u} = {\bf B}^T {\bf \hat{u}} \f$ with 
            *  \f${\bf B}[i][j] = \phi_i(\xi_{j})\f$
            *
            *  This function requires that the coefficient array 
            *  \f$\mathbf{\hat{u}}\f$ (implemented as the attribute #m_coeffs) 
            *  is set.
            *
            *  The resulting array \f$\mathbf{u}[m]=u(\mathbf{\xi}_m)\f$
            *  containing the expansion 
            *  evaluated at the quadrature points, is stored in the 
            *  \a outarray. (Note that the class attribute #m_phys
            *  provides a suitable location to store this result) 
            *
            *  \param outarray contains the values of the expansion evaluated
            *  at the quadrature points (output of the function)
            */
            void  BwdTrans (double *outarray)
            {
                v_BwdTrans (outarray);
            }

            /** \brief This function performs the Forward transformation from 
            *  physical space to coefficient space
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_FwdTrans() 
            *
            *  Given a function evaluated at the quadrature points, this 
            *  function calculates the expansion coefficients such that the 
            *  resulting expansion approximates the original function.
            *
            *  The calculation of the expansion coefficients is done using a 
            *  Galerkin projection. This is equivalent to the operation:
            *  \f[ \mathbf{\hat{u}} = \mathbf{M}^{-1} \mathbf{I}\f]
            *  where
            *  - \f$\mathbf{M}[p][q]= \int\phi_p(\mathbf{\xi})\phi_q(
            *  \mathbf{\xi}) d\mathbf{\xi}\f$ is the Mass matrix
            *  - \f$\mathbf{I}[p] = \int\phi_p(\mathbf{\xi}) u(\mathbf{\xi})    
            *  d\mathbf{\xi}\f$
            *
            *  This function takes the array \a inarray as the values of the 
            *  function evaluated at the quadrature points 
            *  (i.e. \f$\mathbf{u}\f$),
            *  and stores the resulting coefficients \f$\mathbf{\hat{u}}\f$  
            *  in the class attribute #m_coeffs 
            *  
            *  \param inarray array of the function discretely evaluated at the
            *  quadrature points
            */
            void  FwdTrans (const double *inarray)
            {
                v_FwdTrans(inarray);
            }

            /** \brief This function integrates the specified function over the 
            *  domain 
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_Integral()
            *
            *  Based on the values of the function evaluated at the quadrature
            *  points (which are stored in \a inarray), this function calculates
            *  the integral of this function over the domain.  This is 
            *  equivalent to the numerical evaluation of the operation
            *  \f[ I=\int u(\mathbf{\xi})d \mathbf{\xi}\f]
            *
            *  \param inarray values of the function to be integrated evaluated
            *  at the quadrature points (i.e. 
            *  \a inarray[m]=\f$u(\mathbf{\xi}_m)\f$)
            *  \return returns the value of the calculated integral
            */
            double Integral(const double *inarray )
            {
                return v_Integral(inarray);
            }
	    
	    /** \brief This function evaluates the expansion at a single
	     *  (arbitrary) point of the domain
	     *
             *  This function is a wrapper around the virtual function 
             *  \a v_Evaluate()
	     *
	     *  Based on the value of the expansion at the quadrature points,
	     *  this function calculates the value of the expansion at an 
	     *  arbitrary single points (with coordinates \f$ \mathbf{x_c}\f$ 
	     *  given by the pointer \a coords). This operation, equivalent to
	     *  \f[ u(\mathbf{x_c})  = \sum_p \phi_p(\mathbf{x_c}) \hat{u}_p \f] 
	     *  is evaluated using Lagrangian interpolants through the quadrature
	     *  points:
	     *  \f[ u(\mathbf{x_c}) = \sum_p h_p(\mathbf{x_c}) u_p\f]
	     *
             *  This function requires that the physical value array 
             *  \f$\mathbf{u}\f$ (implemented as the attribute #m_phys) 
             *  is set.
	     * 
	     *  \param coords the coordinates of the single point
	     *  \return returns the value of the expansion at the single point
	     */
            double Evaluate(const double * coords)
            {
                return v_Evaluate(coords);
            }

	    /** \brief This function fills the array \a outarray with the 
	     *  \a mode-th mode of the expansion 
	     *
             *  This function is a wrapper around the virtual function 
             *  \a v_FillMode()
	     *
	     *  The requested mode is evaluated at the quadrature points
	     *
	     *  \param mode the mode that should be filled
	     *  \param outarray contains the values of the \a mode-th mode of the
	     *  expansion evaluated at the quadrature points (output of the 
	     *  function)
	     */
            void FillMode(const int mode, double * outarray)
            {
                v_FillMode(mode, outarray);
            }

	    /** \brief this function calculates the inner product of a given 
	     *  function \a f with the different modes of the expansion
	     *  
             *  This function is a wrapper around the virtual function 
             *  \a v_IProductWRTBase()
	     * 
	     *  This is equivalent to the numerical evaluation of 
	     *  \f[ I[p] = \int \phi_p(\mathbf{x}) f(\mathbf{x}) d\mathbf{x}\f]
	     *
	     *  \param inarray contains the values of the function \a f 
	     *  evaluated at the quadrature points
	     *  \param outarray contains the values of the inner product of \a f
	     *  with the different modes, i.e. \f$ outarray[p] = I[p]\f$ 
	     *  (output of the function) 
	     */
            void IProductWRTBase(const double *inarray, double * outarray)
            {
                v_IProductWRTBase(inarray, outarray);
            }

            DNekMatSharedPtr GenNBasisTransMatrix()
            {
                return v_GenNBasisTransMatrix();
            }

            void GetCoords(double **coords)
            {
                v_GetCoords(coords);
            }

            void GetCoord(const double *Lcoords, double *coords)
            {
                v_GetCoord(Lcoords, coords);
            }

            void WriteToFile(FILE *outfile)
            {
                v_WriteToFile(outfile);
            }

            void WriteToFile(std::ofstream &outfile)
            {
                v_WriteToFile(outfile);
            }

            void WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                v_WriteToFile(outfile,dumpVar);
            }

            GeomType MinfoType(void)
            {
                return v_MinfoType();
            }


            // virtual functions related to LocalRegions

            boost::shared_ptr<LocalRegions::MetricRelatedInfo> GenMinfo(void)
            {
                return v_GenMinfo();
            }

            boost::shared_ptr<LocalRegions::MetricRelatedInfo> GetMinfo(void)
            {
                return v_GetMinfo();
            }

            void SetMinfo(boost::shared_ptr<LocalRegions::MetricRelatedInfo> minfo)
            {
                v_SetMinfo(minfo);
            }

            int GetCoordim()
            {
                return v_GetCoordim(); 
            }

            // element boundary ordering 
            // Segment mapping: Vertex to Seg
            void MapTo(EdgeOrientation dir, StdExpMap &Map)
            {
                v_MapTo(dir,Map);
            }

            // EdgeTo2D mapping 
            void  MapTo(const int edge_ncoeff, 
                const LibUtilities::BasisType Btype, 
                const int eid, const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                v_MapTo(edge_ncoeff,Btype,eid,eorient,Map);
            }

            // EdgeTo2D mapping 
            void  MapTo_ModalFormat(const int edge_ncoeff, 
                const LibUtilities::BasisType Btype, 
                const int eid, 
                const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                v_MapTo_ModalFormat(edge_ncoeff,Btype,eid,eorient,Map);
            }

            // Matrix Routines
            DNekMatSharedPtr GenerateMassMatrix();

            DNekMatSharedPtr GenMassMatrix ()
            {
                return v_GenMassMatrix();
            }

            DNekMatSharedPtr GenLapMatrix()
            {
                return v_GenLapMatrix();
            }

            void PhysDeriv (const int dim, double **outarray) 
            {
                v_PhysDeriv (dim, outarray);
            }

            void PhysDeriv (const int dim, const double *inarray, 
                double **outarray)
            {
                v_PhysDeriv (dim, inarray, outarray);
            }

            void Interp1D(const LibUtilities::BasisKey &fbasis0,
			  const double *from,
			  const LibUtilities::BasisKey &tbasis0, double *to);

            void Interp2D(const LibUtilities::BasisKey &fbasis0, 
			  const LibUtilities::BasisKey &fbasis1,
			  const double *from,   
			  const LibUtilities::BasisKey &tbasis0,
			  const LibUtilities::BasisKey &tbasis1, double *to);

            /** \brief Function to evaluate the discrete \f$ L_\infty\f$
            *  error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
            *	u_{exact}\f$ is given by the array \a sol. 
            *
            *	This function takes the physical value space array \a m_phys as
            *  approximate solution
            *
            *  \param sol array of solution function  at physical quadrature
            *  points
            *  \return returns the \f$ L_\infty \f$ error as a double. 
            */
            double Linf(const double *sol);

            /** \brief Function to evaluate the discrete \f$ L_\infty \f$ norm of
            *  the function defined at the physical points \a (this)->m_phys. 
            *
            *	This function takes the physical value space array \a m_phys as
            *  discrete function to be evaluated
            *
            *  \return returns the \f$ L_\infty \f$ norm as a double.
            */
            double Linf();

            /** \brief Function to evaluate the discrete \f$ L_2\f$ error,
            *  \f$ | \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2
            *  dx \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by 
            *  the array \a sol.
            *
            *	This function takes the physical value space array \a m_phys as
            *  approximate solution
            *
            *  \param sol array of solution function  at physical quadrature
            *  points
            *  \return returns the \f$ L_2 \f$ error as a double. 
            */
            double L2(const double *sol);

            /** \brief Function to evaluate the discrete \f$ L_2\f$ norm of the
            *  function defined at the physical points \a (this)->m_phys.
            *
            *	This function takes the physical value space array \a m_phys as
            *  discrete function to be evaluated
            *
            *  \return returns the \f$ L_2 \f$ norm as a double
            */
            double L2();

            // I/O routines
            void WriteCoeffsToFile(std::ofstream &outfile);

        protected:

            DNekMatSharedPtr    CreateStdMatrix(const StdMatrixKey &mkey);
            DNekLinSysSharedPtr CreateStdLinSys(const StdLinSysKey &mkey);

            int   m_numbases;       /**< Number of 1D basis defined in expansion */
            boost::shared_array<LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */
            //LibUtilities::BasisSharedPtr *m_base; /**< Bases needed for the expansion */

            LibUtilities::NekManager<StdMatrixKey, DNekMat,    StdMatrixKey::opLess> m_stdMatrixManager;

            LibUtilities::NekManager<StdLinSysKey, DNekLinSys, StdLinSysKey::opLess> m_stdLinSysManager;
            /** Total number of coefficients used in the expansion*/
            int  m_ncoeffs;
            BstShrDArray m_coeffs;   /**< Array containing expansion coefficients */
            /** Array containing expansion evaluated at the quad points */
            BstShrDArray m_phys;

        private:

            // Virtual functions

            virtual int v_GetNverts() = 0;
            virtual int v_GetNedges() = 0;
            virtual int v_GetNfaces() = 0;

            virtual int v_GetEdgeNcoeffs(const int i)
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i)
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return LibUtilities::eNoBasisType;
            }

            virtual ShapeType v_DetShapeType()  
            {
                ASSERTL0(false, "This expansion does not have a shape type defined");
                return eNoShapeType;
            }

            virtual void   v_BwdTrans (double *outarray)      = 0;
            virtual void   v_FwdTrans (const double *inarray) = 0;

            virtual double v_Integral(const double *inarray ) = 0;
            virtual double v_Evaluate(const double * coords)  = 0;

            virtual void   v_PhysDeriv (const int dim, double **outarray)
            {
                ASSERTL0(false, "This function is only valid for "
                    " local expansions");
            }


            virtual void   v_PhysDeriv (const int dim, const double *inarray,
                double **outarray)
            {
                ASSERTL0(false, "This function is only valid for "
                    "local expansions");
            }

            virtual void v_FillMode(const int mode, double * outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this shape");
            }

            virtual void v_IProductWRTBase(const double *inarray, double * outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this shape");
            }

            virtual DNekMatSharedPtr v_GenMassMatrix()
            {
                DNekMatSharedPtr returnval;

                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this element");

                return returnval;
            }

            virtual DNekMatSharedPtr v_GenLapMatrix()
            {
                DNekMatSharedPtr returnval;

                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this element");

                return returnval;
            }

            virtual DNekMatSharedPtr v_GenNBasisTransMatrix()
            {
                DNekMatSharedPtr returnval;

                NEKERROR(ErrorUtil::efatal, "This function is only valid "
                    "for nodal expansions");

                return returnval;
            }

            virtual void v_GetCoords(double **coords)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            virtual void v_GetCoord(const double *Lcoords, double *coords)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            virtual int v_GetCoordim(void)
            {
                NEKERROR(ErrorUtil::efatal, "Write method");		
                return -1;
            }

            // element boundary ordering 
            virtual void v_MapTo(EdgeOrientation dir, StdExpMap &Map)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );		
            }

            virtual void  v_MapTo(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype,
                const int eid, const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );		
            }

            virtual void  v_MapTo_ModalFormat(const int edge_ncoeffs, 
                const LibUtilities::BasisType Btype,
                const int eid, 
                const EdgeOrientation eorient, 
                StdExpMap &Map)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );		
            }


            virtual void v_WriteToFile(FILE *outfile)
            {
                NEKERROR(ErrorUtil::efatal, "WriteToFile: Write method");
            }

            virtual void v_WriteToFile(std::ofstream &outfile)
            {
                NEKERROR(ErrorUtil::efatal, "WriteToFile:Write method");
            }

            virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                NEKERROR(ErrorUtil::efatal, "WriteToFile: Write method");
            }

            virtual GeomType v_MinfoType()
            {
                return eRegular;
            }

            // virtual functions related to LocalRegions
            virtual boost::shared_ptr<LocalRegions::MetricRelatedInfo> v_GenMinfo()
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<LocalRegions::MetricRelatedInfo>();
            }

            virtual boost::shared_ptr<LocalRegions::MetricRelatedInfo> v_GetMinfo()
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<LocalRegions::MetricRelatedInfo>();
            }

            virtual void v_SetMinfo(boost::shared_ptr<LocalRegions::MetricRelatedInfo>
                minfo)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            }
        };

        typedef boost::shared_ptr<StdExpansion> StdExpansionSharedPtr;
        typedef std::vector< StdExpansionSharedPtr > StdExpansionVector;
        typedef std::vector< StdExpansionSharedPtr >::iterator StdExpansionVectorIter;

    } //end of namespace
} //end of namespace

#endif //STANDARDDEXPANSION_H
/**
* $Log: StdExpansion.h,v $
* Revision 1.28  2007/03/02 12:01:52  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.27  2007/03/01 17:04:07  jfrazier
* Removed extraneous basis.
*
* Revision 1.26  2007/03/01 03:52:10  jfrazier
* Added GetBasis function.
*
* Revision 1.25  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.24  2007/02/28 09:53:17  sherwin
* Update including adding GetBasis call to StdExpansion
*
* Revision 1.23  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.22  2007/02/23 19:26:08  jfrazier
* General bug fix and formatting.
*
* Revision 1.21  2007/02/22 22:02:28  sherwin
* Update with executing StdMatManager
*
* Revision 1.20  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.19  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.18  2007/02/17 04:03:23  jfrazier
* Added NekManager for holding matrices.  Need to finish the create function.
*
* Revision 1.17  2007/02/16 17:14:39  pvos
* Added documentation
*
* Revision 1.16  2007/02/07 12:51:53  sherwin
* Compiling version of Project1D
*
* Revision 1.15  2007/02/06 02:23:31  jfrazier
* Minor cleanup.
*
* Revision 1.14  2007/01/28 18:34:21  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.13  2007/01/23 23:20:21  sherwin
* New version after Jan 07 update
*
* Revision 1.12  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.11  2007/01/18 23:03:56  sherwin
* Removed for repository update in utah 07
*
* Revision 1.10  2007/01/15 11:08:40  pvos
* Updating doxygen documentation
*
* Revision 1.9  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.8  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.7  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.6  2006/06/13 18:05:02  sherwin
* Modifications to make MultiRegions demo ProjectLoc2D execute properly.
*
* Revision 1.5  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.4  2006/06/01 13:43:19  kirby
* *** empty log message ***
*
* Revision 1.3  2006/05/30 14:00:04  sherwin
* Updates to make MultiRegions and its Demos work
*
* Revision 1.2  2006/05/29 19:03:08  sherwin
* Modifications to wrap geometric information in shared_ptr
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.75  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.74  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.73  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.72  2006/03/05 23:17:53  sherwin
*
* Corrected to allow MMatrix1D and MMatrix2D to execute properly
*
* Revision 1.71  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.70  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.69  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.66  2006/03/01 22:59:12  sherwin
*
* First working version of Project1D
*
* Revision 1.65  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.64  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.63  2006/02/26 21:23:20  bnelson
* Fixed a variety of compiler errors caused by updates to the coding standard.
*
* Revision 1.62  2006/02/19 13:26:13  sherwin
*
* Coding standard revisions so that libraries compile
*
**/




