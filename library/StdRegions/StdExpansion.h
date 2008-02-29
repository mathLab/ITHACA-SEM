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
            StdExpansion(const int numcoeffs, const int numbases, 
                const LibUtilities::BasisKey &Ba, 
                const LibUtilities::BasisKey &Bb = LibUtilities::NullBasisKey,
                const LibUtilities::BasisKey &Bc = LibUtilities::NullBasisKey);

            /** \brief Copy Constructor */
            StdExpansion(const StdExpansion &T);

            /** \brief Destructor */
            virtual ~StdExpansion();


            /** \brief This function returns a pointer to the coefficient array
            *  \f$ \mathbf{\hat{u}}\f$ 
            *
            *  The coefficient array \f$ \mathbf{\hat{u}}\f$ corresponds to the 
            *  class attribute #m_coeffs (which is in coefficient space)
            *
            *  \return returns a pointer to the coefficient array 
            *  \f$ \mathbf{\hat{u}}\f$ 
            */
            inline const ConstArray<OneD, NekDouble>& GetCoeffs(void) const
            {
                return m_coeffs;
            }


            /** \brief This function returns a NekDouble to the coefficient 
            *  \f$ \mathbf{\hat{u}}[i]\f$ 
            *
            *  The coefficient \f$ \mathbf{\hat{u}}[i]\f$ corresponds
            *  to the #i th entry of the class attribute #m_coeffs
            *
            *  \return returns a NekDouble of the coefficient 
            *  \f$ \mathbf{\hat{u}}[i]\f$ 
            */
            inline NekDouble  GetCoeffs(int i) const
            {
                ASSERTL1(i < m_ncoeffs,"index out of range");

                return m_coeffs[i];
            }


            /** \brief This function returns a NekDouble to the coefficient 
            *  \f$ \mathbf{\hat{u}}[i]\f$ 
            *
            *  The coefficient \f$ \mathbf{\hat{u}}[i]\f$ corresponds
            *  to the #i th entry of the class attribute #m_coeffs
            *
            *  \return returns a NekDouble of the coefficient 
            *  \f$ \mathbf{\hat{u}}[i]\f$ 
            */
            inline NekDouble  GetCoeff(int i) const
            {
                ASSERTL1(i < m_ncoeffs,"index out of range");

                return m_coeffs[i];
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
            inline const ConstArray<OneD, NekDouble>& GetPhys(void) const
            {
                return m_phys;
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
            inline Array<OneD, NekDouble>& UpdateCoeffs(void)
            {
                return(m_coeffs);
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
            inline Array<OneD, NekDouble>& UpdatePhys(void) 
            {
                return(m_phys);
            }

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
            inline const ConstArray<OneD, LibUtilities::BasisSharedPtr>& GetBase() const
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
                ASSERTL1(dir < m_numbases, "dir is larger than number of bases");
                return(m_base[dir]);
            }

            /** \brief This function returns the total number of coefficients 
            *  used in the expansion 
            *  
            *  \return returns the total number of coefficients (which is 
            *  equivalent to the total number of modes) used in the expansion
            */
            inline int GetNcoeffs(void) const 
            {
                return(m_ncoeffs);
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
            inline void SetCoeffs(const ConstArray<OneD, NekDouble>& coeffs)
            {
                Vmath::Vcopy(m_ncoeffs, &coeffs[0], 1, &m_coeffs[0], 1);
            }

            /** \brief This function sets the i th coefficient  
            *  \f$ \mathbf{\hat{u}}[i]\f$ to the value given by \a coeff
            *
            *  #m_coeffs[i] will be set to the  value given by \a coeff
            *
            *  \param i the index of the coefficient to be set
            *  \param coeff the value of the coefficient to be set
            */
            inline void SetCoeff(const int i, const NekDouble coeff)
            {
                m_coeffs[i] = coeff;
            }

            /** \brief This function returns the total number of quadrature
            *  points used in the element 
            *  
            *  \return returns the total number of quadrature points
            */
            inline  int GetTotPoints() const
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
            inline void SetPhys(const ConstArray<OneD, NekDouble>& phys)
            {
                int nqtot = GetTotPoints();

                Vmath::Vcopy(nqtot, &phys[0], 1, &m_phys[0], 1);
            }

            /** \brief This function returns the type of basis used in the \a dir
            *  direction
            *  
            *  The different types of bases implemented in the code are defined 
            *  in the LibUtilities::BasisType enumeration list. As a result, the
            *  function will return one of the types of this enumeration list.
            *  
            *  \param dir the direction
            *  \return returns the type of basis used in the \a dir direction
            */
            inline  LibUtilities::BasisType GetBasisType(const int dir) const
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
            *  As a result, the function will return one of the types of this 
            *  enumeration list.
            *  
            *  \param dir the direction
            *  \return returns the type of quadrature points  used in the \a dir
            *  direction
            */
            inline LibUtilities::PointsType GetPointsType(const int dir)  const 
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
            inline const ConstArray<OneD, NekDouble>& GetPoints(const int dir) const
            {
                return ExpPointsProperties(dir)->GetZ();
            }


            NekDouble operator[] (const int i) const
            {
                ASSERTL1((i >= 0) && (i < m_ncoeffs),
                    "Invalid Index used in [] operator");
                return m_coeffs[i];
            }

            NekDouble& operator[](const int i)
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
            int GetNverts() const
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
            int GetNedges() const
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
            int GetEdgeNcoeffs(const int i) const
            {
                return v_GetEdgeNcoeffs(i);
            }


            int NumBndryCoeffs(void)  const
            {
                return v_NumBndryCoeffs();
            }

            int NumDGBndryCoeffs(void)  const
            {
                return v_NumDGBndryCoeffs();
            }

            /** \brief This function returns the type of expansion basis on the
            *  \a i-th edge  
            *  
            *  This function is a wrapper around the virtual function 
            *  \a v_GetEdgeBasisType() 
            *
            *  The different types of bases implemented in the code are defined 
            *  in the LibUtilities::BasisType enumeration list. As a result, the
            *  function will return one of the types of this enumeration list.
            * 
            *  \param i specifies which edge
            *  \return returns the expansion basis on the \a i-th edge
            */
            LibUtilities::BasisType GetEdgeBasisType(const int i) const
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
            int GetNfaces() const
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
            *  function will return one of the types of this enumeration list.
            *  
            *  \return returns the shape of the expansion domain
            */        
            ShapeType DetShapeType() const
            {
                return v_DetShapeType();
            }


            bool IsBoundaryInteriorExpansion() 
            {
                return v_IsBoundaryInteriorExpansion();
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
            *  \f$\mathbf{\hat{u}}\f$ provided as \a inarray. 
            *
            *  The resulting array
            *  \f$\mathbf{u}[m]=u(\mathbf{\xi}_m)\f$ containing the
            *  expansion evaluated at the quadrature points, is stored
            *  in the \a outarray. (Note that the class attribute
            *  #m_phys provides a suitable location to store this
            *  result)
            *
            *  \param inarray contains the values of the expansion
            *  coefficients (input of the function)
            *
            *  \param outarray contains the values of the expansion evaluated
            *  at the quadrature points (output of the function)
            */

            void  BwdTrans (const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
            {
                v_BwdTrans (inarray, outarray);
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
            *  in the \a outarray
            *  
            *  \param inarray array of the function discretely evaluated at the
            *  quadrature points
            *
            *  \param outarray array of the function coefficieints 
            */
            void  FwdTrans (const ConstArray<OneD, NekDouble>& inarray, 
                            Array<OneD, NekDouble> &outarray)
            {
                v_FwdTrans(inarray,outarray);
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
            NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                return v_Integral(inarray);
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
            void FillMode(const int mode, Array<OneD, NekDouble> &outarray)
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

            void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, 
                 Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTBase(inarray, outarray);
            }

            void   IProductWRTDerivBase(const int dir, 
                                         const ConstArray<OneD,NekDouble>& inarray, 
                                         Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDerivBase(dir,inarray, outarray);
            }

            /** \brief this function returns the physical coordinates of the
            *  quadrature points of the expansion
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_GetCoords()
            *
            *  \param coords an array containing the coordinates of the
            *  quadrature points (output of the function)
            */
            void GetCoords(Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray,
                Array<OneD, NekDouble> &coords_3 = NullNekDouble1DArray)
            {
                v_GetCoords(coords_1,coords_2,coords_3);
            }

            /** \brief given the coordinates of a point of the element in the 
            *  local collapsed coordinate system, this function calculates the 
            *  physical coordinates of the point
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_GetCoord()         
            *
            *  \param Lcoords the coordinates in the local collapsed 
            *  coordinate system
            *  \param coords the physical coordinates (output of the function)
            */
            void GetCoord(const ConstArray<OneD, NekDouble>& Lcoord, 
              Array<OneD, NekDouble> &coord)
            {
                v_GetCoord(Lcoord, coord);
            }

            /** \brief this function writes the solution to the file \a outfile
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_WriteToFile()
            *
            *  The expansion evaluated at the quadrature points (stored as 
            *  #m_phys), together with 
            *  the coordinates of the quadrature points, are written to the 
            *  file \a outfile
            *  
            *  \param outfile the file to which the solution is written
            */
            void WriteToFile(FILE *outfile)
            {
                v_WriteToFile(outfile);
            }

            /** \brief this function writes the solution to the file \a outfile
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_WriteToFile()
            *
            *  The expansion evaluated at the quadrature points (stored as 
            *  #m_phys), together with 
            *  the coordinates of the quadrature points, are written to the 
            *  file \a outfile
            *  
            *  \param outfile the file to which the solution is written
            */
            void WriteToFile(std::ofstream &outfile)
            {
                v_WriteToFile(outfile);
            }

            /** \brief this function writes the solution to the file \a outfile
            *
            *  This function is a wrapper around the virtual function 
            *  \a v_WriteToFile()
            *
            *  The expansion evaluated at the quadrature points (stored as 
            *  #m_phys), together with 
            *  the coordinates of the quadrature points, are written to the 
            *  file \a outfile
            *  
            *  \param outfile the file to which the solution is written
            */
            void WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                v_WriteToFile(outfile,dumpVar);
            }
            
            DNekScalMatSharedPtr GetLocMatrix(const MatrixType mtype, NekDouble lambdaval = NekUnsetDouble, NekDouble tau = NekUnsetDouble)
            {
                return v_GetLocMatrix(mtype,lambdaval,tau);
            }
            

            DNekScalMatSharedPtr GetLocMatrix(const LocalRegions::MatrixKey &mkey)
            {
                return v_GetLocMatrix(mkey);
            }


            DNekScalBlkMatSharedPtr GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
            {
                return v_GetLocStaticCondMatrix(mkey);
            }

            StdRegions::EdgeOrientation GetEorient(int edge)
            {
                return v_GetEorient(edge);
            }

            void AddUDGHelmholtzBoundaryTerms(const NekDouble tau, 
                                              const ConstArray<OneD,NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray)
            {
                v_AddUDGHelmholtzBoundaryTerms(tau,inarray,outarray);
            }
            

            void AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                           const ConstArray<OneD,NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
                
            {
                v_AddUDGHelmholtzTraceTerms(tau,inarray,outarray);
            }

            void AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                           const ConstArray<OneD,NekDouble> &inarray,
                                           Array<OneD,boost::shared_ptr<LocalRegions::SegExp> > &EdgeExp,
                                           Array<OneD,NekDouble> &outarray,
                                           bool UseLocalOrient = true)
                
            {
                v_AddUDGHelmholtzTraceTerms(tau,inarray,EdgeExp,
                                            outarray,UseLocalOrient);
            }

            virtual void AddNormBoundaryInt(const int dir,
                                            ConstArray<OneD,NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray)
            {
                v_AddNormBoundaryInt(dir,inarray,outarray);
            }


            virtual void AddNormTraceInt(const int dir,
                                         ConstArray<OneD,NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray)
            {
                v_AddNormTraceInt(dir,inarray,outarray);
            }
            // virtual functions related to LocalRegions

            int GetCoordim()
            {
                return v_GetCoordim(); 
            }

            const ConstArray<OneD, int> GetBoundaryMap(void)
            {
                return v_GetBoundaryMap();
            }

            const ConstArray<OneD, int> GetInteriorMap(void)
            {
                return v_GetInteriorMap();
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

            /** \brief this function generates the mass matrix 
            *  \f$\mathbf{M}[i][j] =
            *  \int \phi_i(\mathbf{x}) \phi_j(\mathbf{x}) d\mathbf{x}\f$
            * 
            *  \return returns the mass matrix
            */
            DNekMatSharedPtr CreateGeneralMatrix(const StdMatrixKey &mkey);

            void GeneralMatrixOp(const StdMatrixKey &mkey, 
                                 const ConstArray<OneD,NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray);

            
            void MassMatrixOp(const ConstArray<OneD,NekDouble> &inarray, 
                              Array<OneD,NekDouble> &outarray);
        
            void LaplacianMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            void LaplacianMatrixOp(const int k1, const int k2, 
                                   const ConstArray<OneD,NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            void WeakDerivMatrixOp(const int i,
                                   const ConstArray<OneD,NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);


            void BwdTransMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray);

            void HelmholtzMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const double lambda);

            void TensProdBwdTrans(const ConstArray<OneD, NekDouble>&inarray, 
                                  Array<OneD, NekDouble> &outarray);

            void TensProdBwdTrans(const StdExpansion &in)
            {
                TensProdBwdTrans(((StdExpansion &)in).GetCoeffs(), m_phys);
            }


            DNekMatSharedPtr GenMatrix (const StdMatrixKey &mkey)
            {
                return v_GenMatrix(mkey);
            }

            void PhysDeriv (const ConstArray<OneD, NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                            Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                v_PhysDeriv (inarray, out_d0, out_d1, out_d1);
            }

            virtual void PhysDeriv(const int dir, 
                                   const ConstArray<OneD, NekDouble>& inarray,
                                   Array<OneD, NekDouble> &outarray)
            {
                v_PhysDeriv (dir, inarray, outarray);
            }

            boost::shared_ptr<SpatialDomains::GeomFactors> GetMetricInfo(void) const
            {
                return v_GetMetricInfo();
            }

            virtual DNekScalMatSharedPtr v_GetLocMatrix(const LocalRegions::MatrixKey &mkey)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<DNekScalMat>();
            }

            virtual DNekScalMatSharedPtr v_GetLocMatrix(const StdRegions::MatrixType mtype, NekDouble lambdaval, NekDouble tau)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<DNekScalMat>();
            }


            virtual DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<DNekScalBlkMat>();

            }

            virtual StdRegions::EdgeOrientation v_GetEorient(int edge)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for two-dimensional  LocalRegions");                
            }


            virtual void v_AddUDGHelmholtzBoundaryTerms(const NekDouble tau, 
                                                        const ConstArray<OneD,NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

            virtual void v_AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                     const ConstArray<OneD,NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray)
            { 
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

            virtual void v_AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                    const ConstArray<OneD,NekDouble> &inarray,
                                    Array<OneD, boost::shared_ptr<LocalRegions::SegExp> > &edgeExp, 
                                    Array<OneD,NekDouble> &outarray,
                                    bool UseLocalOrient = true)
            { 
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }
            
            virtual void v_AddNormBoundaryInt(const int dir,
                                              ConstArray<OneD,NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

            virtual void v_AddNormTraceInt(const int dir,
                                           ConstArray<OneD,NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

            /** \brief this function interpolates a 1D function \f$f\f$ evaluated
            *  at the quadrature points of the basis \a fbasis0 to the 
            *  function values at the quadrature points of the basis \a tbasis0
            *
            *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
            *  of the basis \a fbasis0, this routine calculates, using 
            *  \a (Q-1)th order polynomial interpolation, the function values
            *  at the \a Q2 quadrature points of the basis \a tbasis0.
            *
            *  \param fbasis0 the basis at which's quadrature points the 
            *  function is given
            *  \param from the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a fbasis0
            *  \param tbasis0 the basis to which's quadrature points the 
            *  function should be interpolated
            *  \param to the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a tbasis0 (output of the function)
            */
            void Interp1D(const LibUtilities::BasisKey &fbasis0,
              const ConstArray<OneD, NekDouble>& from,
              const LibUtilities::BasisKey &tbasis0,
              Array<OneD, NekDouble> &to);

            void Interp1D(const LibUtilities::BasisKey &fbasis0,
                          const NekDouble *from,
                          const LibUtilities::BasisKey &tbasis0,
                          NekDouble *to);


            /** \brief this function interpolates a 2D function \f$f\f$ evaluated
            *  at the quadrature points of the 2D basis, constructed by 
            *  \a fbasis0 and \a fbasis1, to the function values at the 
            *  quadrature points of the 2D basis, constructed by \a tbasis0 and 
            *  \a tbasis1
            *
            *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
            *  of the first expansion basis, this routine calculates, using 
            *  \a (Q-1)th order polynomial interpolation, the function values
            *  at the \a Q2 quadrature points of the second basis.
            *
            *  \param fbasis0 the basis at which's quadrature points the 
            *  function is given
            *  \param from the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a fbasis0
            *  \param tbasis0 the basis to which's quadrature points the 
            *  function should be interpolated
            *  \param to the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a tbasis0 (output of the function)
            */
            void Interp2D(const LibUtilities::BasisKey &fbasis0, 
                          const LibUtilities::BasisKey &fbasis1,
                          const ConstArray<OneD, NekDouble>& from,
                          const LibUtilities::BasisKey &tbasis0,
                          const LibUtilities::BasisKey &tbasis1, 
                          Array<OneD, NekDouble> &to);

            void Interp2D(const LibUtilities::BasisKey &fbasis0, 
                          const LibUtilities::BasisKey &fbasis1,
                          const NekDouble *from,   
                          const LibUtilities::BasisKey &tbasis0,
                          const LibUtilities::BasisKey &tbasis1, 
                          NekDouble *to);


            /** \brief this function interpolates a 3D function \f$f\f$ evaluated
            *  at the quadrature points of the 3D basis, constructed by 
            *  \a fbasis0, \a fbasis1, and \a fbasis2 to the function values at the 
            *  quadrature points of the 3D basis, constructed by \a tbasis0, 
            *  \a tbasis1, and \a tbasis2.
            *
            *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
            *  of the first expansion basis, this routine calculates, using 
            *  \a (Q-1)th order polynomial interpolation, the function values
            *  at the \a Q2 quadrature points of the second basis, and the function
            *  values at the \a Q3 quadrature points of the third basis.
            *
            *  \param fbasis0 the basis at which's quadrature points the 
            *  function is given
            *  \param from the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a fbasis0
            *  \param tbasis0 the basis to which's quadrature points the 
            *  function should be interpolated
            *  \param to the array containg the function \f$ f\f$  evaluated
            *   at the quadrature points of \a tbasis0 (output of the function)
            */
           void Interp3D(const LibUtilities::BasisKey &fbasis0,
                         const LibUtilities::BasisKey &fbasis1,
                         const LibUtilities::BasisKey &fbasis2,
                         const ConstArray<OneD, NekDouble>& from,
                         const LibUtilities::BasisKey &tbasis0,
                         const LibUtilities::BasisKey &tbasis1,
                         const LibUtilities::BasisKey &tbasis2,
                         Array<OneD, NekDouble> &to);

           void Interp3D(const LibUtilities::BasisKey &fbasis0,
                         const LibUtilities::BasisKey &fbasis1,
                         const LibUtilities::BasisKey &fbasis2,
                         const NekDouble *from,
                         const LibUtilities::BasisKey &tbasis0,
                         const LibUtilities::BasisKey &tbasis1,
                         const LibUtilities::BasisKey &tbasis2,
                         NekDouble *to );


            /** \brief Function to evaluate the discrete \f$ L_\infty\f$
            *  error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
            *    u_{exact}\f$ is given by the array \a sol. 
            *
            *    This function takes the physical value space array \a m_phys as
            *  approximate solution
            *
            *  \param sol array of solution function  at physical quadrature
            *  points
            *  \return returns the \f$ L_\infty \f$ error as a NekDouble. 
            */
            NekDouble Linf(const ConstArray<OneD, NekDouble>& sol);

            /** \brief Function to evaluate the discrete \f$ L_\infty \f$ norm of
            *  the function defined at the physical points \a (this)->m_phys. 
            *
            *    This function takes the physical value space array \a m_phys as
            *  discrete function to be evaluated
            *
            *  \return returns the \f$ L_\infty \f$ norm as a double.
            */
            NekDouble Linf();

            /** \brief Function to evaluate the discrete \f$ L_2\f$ error,
            *  \f$ | \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2
            *  dx \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by 
            *  the array \a sol.
            *
            *    This function takes the physical value space array \a m_phys as
            *  approximate solution
            *
            *  \param sol array of solution function  at physical quadrature
            *  points
            *  \return returns the \f$ L_2 \f$ error as a double. 
            */
            NekDouble L2(const ConstArray<OneD, NekDouble>& sol);

            /** \brief Function to evaluate the discrete \f$ L_2\f$ norm of the
            *  function defined at the physical points \a (this)->m_phys.
            *
            *    This function takes the physical value space array \a m_phys as
            *  discrete function to be evaluated
            *
            *  \return returns the \f$ L_2 \f$ norm as a double
            */
            NekDouble L2();

            // I/O routines
            void WriteCoeffsToFile(std::ofstream &outfile);


        protected:

            int   m_numbases;   /**< Number of 1D basis defined in expansion */
            Array<OneD, LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */

            /** Total number of coefficients used in the expansion*/
            int  m_ncoeffs;
            Array<OneD, NekDouble> m_coeffs;   /**< Array containing expansion coefficients */
            /** Array containing expansion evaluated at the quad points */
            Array<OneD, NekDouble> m_phys;

            bool StdMatManagerAlreadyCreated(const StdMatrixKey &mkey)
            {
                return m_stdMatrixManager.AlreadyCreated(mkey);
            }
         
            bool StdStaticCondMatManagerAlreadyCreated(const StdMatrixKey &mkey)
            {
                return m_stdStaticCondMatrixManager.AlreadyCreated(mkey);
            }
         
            DNekMatSharedPtr CreateStdMatrix(const StdMatrixKey &mkey);

            /** \brief Create the static condensation of a matrix when
                using a boundary interior decomposition

                If a matrix system can be represented by 
                \f$ Mat = \left [ \begin{array}{cc}
                              A & B \\ 
                              C & D \end{array} \right ] \f$
               This routine creates a matrix containing the statically
               condense system of the form
               \f$ Mat = \left [ \begin{array}{cc}
                              A - B D^{-1} C & B D^{-1} \\ 
                              D^{-1} C       & D^{-1} \end{array} \right ] \f$
            **/
            DNekBlkMatSharedPtr CreateStdStaticCondMatrix(const StdMatrixKey &mkey);

            DNekMatSharedPtr GetStdMatrix(const StdMatrixKey &mkey) 
            {
                return m_stdMatrixManager[mkey];
            }

            DNekBlkMatSharedPtr GetStdStaticCondMatrix(const StdMatrixKey &mkey) 
            {
                return m_stdStaticCondMatrixManager[mkey];
            }


        private:

            LibUtilities::NekManager<StdMatrixKey, DNekMat, StdMatrixKey::opLess> m_stdMatrixManager;
            LibUtilities::NekManager<StdMatrixKey, DNekBlkMat, StdMatrixKey::opLess> m_stdStaticCondMatrixManager;

            // Virtual functions

            virtual int v_GetNverts() const = 0;
            virtual int v_GetNedges() const = 0;
            virtual int v_GetNfaces() const = 0;

            virtual int v_NumBndryCoeffs() const 
            {
                ASSERTL0(false, "This function is needs defining for this shape");
                return 0;
            }
            
            virtual int v_NumDGBndryCoeffs() const 
            {
                ASSERTL0(false, "This function is needs defining for this shape");
                return 0;
            }
            
            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return LibUtilities::eNoBasisType;
            }

            virtual ShapeType v_DetShapeType() const
            {
                ASSERTL0(false, "This expansion does not have a shape type defined");
                return eNoShapeType;
            }

            virtual bool  v_IsBoundaryInteriorExpansion() 
            {
                ASSERTL0(false,"This function has not been defined for this expansion");
                return false;
            }


            virtual void   v_BwdTrans   (const ConstArray<OneD, NekDouble>& inarray, 
                     Array<OneD, NekDouble> &outarray) = 0;
            virtual void   v_FwdTrans   (const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray) = 0;
            virtual void   v_IProductWRTBase (const ConstArray<OneD, NekDouble>& inarray, 
                          Array<OneD, NekDouble> &outarray) = 0;

            virtual void   v_IProductWRTDerivBase (const int dir, 
                                                   const ConstArray<OneD,NekDouble>& inarray, 
                                                   Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This method has not been defined");
            }

            virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                    "local expansions");
                return 0;
            }


            virtual void   v_PhysDeriv (const ConstArray<OneD, NekDouble>& inarray,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2,
                                        Array<OneD, NekDouble> &out_d3)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                    "local expansions");
            }

            virtual void v_PhysDeriv(const int dir, 
                const ConstArray<OneD, NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_d0)

            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "specific element types");
            }

            virtual void   v_StdPhysDeriv (const int dir, 
                                           const ConstArray<OneD, NekDouble>& inarray, 
                                           Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is valid in StdExpansion ");
            }
            

            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this shape");
            }

            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey)  
            {
                DNekMatSharedPtr returnval;

                NEKERROR(ErrorUtil::efatal, "This function is has not "
                    "been defined for this element");

                return returnval;
            }
            
            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            virtual void v_GetCoord(const ConstArray<OneD, NekDouble>& Lcoord, 
                Array<OneD, NekDouble> &coord)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            virtual int v_GetCoordim(void)
            {
                NEKERROR(ErrorUtil::efatal, "Write method");        
                return -1;
            }

            virtual const ConstArray<OneD, int> v_GetBoundaryMap(void)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
                return ConstArray<OneD, int>();
            }

            virtual const ConstArray<OneD, int> v_GetInteriorMap(void)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
                return ConstArray<OneD, int>();
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

            virtual boost::shared_ptr<SpatialDomains::GeomFactors> v_GetMetricInfo() const 
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<SpatialDomains::GeomFactors>();

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
* Revision 1.74  2008/02/16 06:00:37  ehan
* Added interpolation 3D.
*
* Revision 1.73  2008/01/23 09:09:46  sherwin
* Updates for Hybrized DG
*
* Revision 1.72  2008/01/03 04:15:17  bnelson
* Added a return value to some functions that abort to prevent visual studio compile errors.
*
* Revision 1.71  2007/12/17 13:03:50  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.70  2007/12/06 22:44:47  pvos
* 2D Helmholtz solver updates
*
* Revision 1.69  2007/11/29 21:40:22  sherwin
* updates for MultiRegions and DG solver
*
* Revision 1.68  2007/11/08 16:55:13  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.67  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.66  2007/08/29 23:26:49  jfrazier
* Created non-static manager that shares data across instances.
*
* Revision 1.65  2007/08/11 23:42:25  sherwin
* A few changes
*
* Revision 1.64  2007/07/27 16:56:55  jfrazier
* Changed manager to static.
*
* Revision 1.63  2007/07/22 23:04:25  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.62  2007/07/20 02:16:52  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.61  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.60  2007/07/13 09:02:25  sherwin
* Mods for Helmholtz solver
*
* Revision 1.59  2007/07/12 12:55:15  sherwin
* Simplified Matrix Generation
*
* Revision 1.57  2007/07/11 19:29:52  sherwin
* Update for ScalMat
*
* Revision 1.56  2007/07/11 13:44:08  kirby
* *** empty log message ***
*
* Revision 1.55  2007/07/11 13:40:26  kirby
* *** empty log message ***
*
* Revision 1.51  2007/07/10 20:41:46  kirby
* more fixes
*
* Revision 1.50  2007/07/10 19:27:57  kirby
* Update for new matrix structures
*
* Revision 1.49  2007/05/30 20:49:13  sherwin
* Updates to do with LocalRegions and SpatialDomains
*
* Revision 1.48  2007/05/28 16:15:01  sherwin
* Updated files in MultiRegions to make 1D demos work
*
* Revision 1.47  2007/05/17 17:59:27  sherwin
* Modification to make Demos work after introducion of Array<>
*
* Revision 1.46  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.45  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.44  2007/04/18 16:09:13  pvos
* Added some new Tensor Operations routines
*
* Revision 1.43  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.42  2007/04/08 03:36:58  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.41  2007/04/06 08:44:43  sherwin
* Update to make 2D regions work at StdRegions level
*
* Revision 1.40  2007/04/04 21:49:25  sherwin
* Update for SharedArray
*
* Revision 1.39  2007/04/04 20:48:16  sherwin
* Update to handle SharedArrays
*
* Revision 1.38  2007/03/31 00:04:03  bnelson
* *** empty log message ***
*
* Revision 1.37  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.36  2007/03/25 15:48:22  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.35  2007/03/21 20:56:43  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.34  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.33  2007/03/20 09:12:46  kirby
* update of geofac and metric info; fix style issues
*
* Revision 1.32  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.31  2007/03/08 09:34:18  pvos
* added documentation
*
* Revision 1.30  2007/03/05 08:07:12  sherwin
* Modified so that StdMatrixKey has const calling arguments in its constructor.
*
* Revision 1.29  2007/03/02 16:43:44  pvos
* Added some documentation
*
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




