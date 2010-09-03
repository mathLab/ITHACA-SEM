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
#include <StdRegions/StdMatrixKey.h>
#include <StdRegions/StdLinSysKey.hpp>
#include <StdRegions/OptimizationParametersAccess.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion1D;

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

#if 1

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
            inline const Array<OneD, const NekDouble>& GetCoeffs(void) const
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
            inline const Array<OneD, const NekDouble>& GetPhys(void) const
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
            inline const Array<OneD, const LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return(m_base);
            }

            /** \brief This function gets the shared point to basis in
             *  the \a dir direction
             *
             *  \return returns the shared pointer to the basis in
             *  directin \a dir
             */
            inline const LibUtilities::BasisSharedPtr& GetBasis(int dir) const
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
            inline void SetCoeffs(const Array<OneD, const NekDouble>& coeffs)
            {
                Vmath::Vcopy(m_ncoeffs, coeffs.get(), 1, m_coeffs.get(), 1);
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
                    nqtot *= m_base[i]->GetNumPoints();
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
            inline void SetPhys(const Array<OneD, const NekDouble>& phys)
            {
                int nqtot = GetTotPoints();

                Vmath::Vcopy(nqtot, phys.get(), 1, m_phys.get(), 1);
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


            /** \brief This function returns the maximum number of
             *  expansion modes over all local directions
             *
             *  \return returns the maximum number of expansion modes
             *  over all local directions
             */
            inline int EvalBasisNumModesMax(void) const
            {
                int i;
                int returnval = 0;

                for(i = 0; i < m_numbases; ++i)
                {
                    returnval = max(returnval, m_base[i]->GetNumModes());
                }

                return returnval;
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
                return(m_base[dir]->GetNumPoints());
            }

            /** \brief This function returns a pointer to the array containing
             *  the quadrature points in \a dir direction
             *
             *  \param dir the direction
             *  \return returns a pointer to the array containing
             *  the quadrature points in \a dir direction
             */
            inline const Array<OneD, const NekDouble>& GetPoints(const int dir) const
            {
                return m_base[dir]->GetZ();
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


            /** \brief This function returns the number of quadrature points
             *  belonging to the \a i-th edge
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetEdgeNumPoints()
             *
             *  \param i specifies which edge
             *  \return returns the number of expansion coefficients belonging to
             *  the \a i-th edge
             */
            int GetEdgeNumPoints(const int i) const
            {
                return v_GetEdgeNumPoints(i);
            }


            int DetCartesianDirOfEdge(const int edge)
            {
                return v_DetCartesianDirOfEdge(edge);
            }

            const LibUtilities::BasisKey DetEdgeBasisKey(const int i) const
            {
                return v_DetEdgeBasisKey(i);
            }

            /** \brief This function returns the number of expansion coefficients
             *  belonging to the \a i-th face
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetFaceNcoeffs()
             *
             *  \param i specifies which face
             *  \return returns the number of expansion coefficients belonging to
             *  the \a i-th face
             */
            int GetFaceNcoeffs(const int i) const
            {
                return v_GetFaceNcoeffs(i);
            }

            int GetFaceIntNcoeffs(const int i) const
            {
                return v_GetFaceIntNcoeffs(i);
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
            ExpansionType DetExpansionType() const
            {
                return v_DetExpansionType();
            }

            int GetShapeDimension() const
            {
                return v_GetShapeDimension();
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

            void  BwdTrans (const Array<OneD, const NekDouble>& inarray,
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
            void  FwdTrans (const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
            {
                v_FwdTrans(inarray,outarray);
            }

            void FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray)
            {
                v_FwdTrans_BndConstrained(inarray,outarray);
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
            NekDouble Integral(const Array<OneD, const NekDouble>& inarray )
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
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTBase(inarray, outarray);
            }

            void   IProductWRTDerivBase(const int dir,
                                        const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDerivBase(dir,inarray, outarray);
            }

            /// \brief Get the element id of this expansion when used
            /// in a list by returning value of #m_elmt_id
            inline int GetElmtId()
            {
                return m_elmt_id;
            }


            /// \brief Set the element id of this expansion when used
            /// in a list by returning value of #m_elmt_id
            inline void SetElmtId(const int id)
            {
                m_elmt_id = id;
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
            void GetCoord(const Array<OneD, const NekDouble>& Lcoord,
                          Array<OneD, NekDouble> &coord)
            {
                v_GetCoord(Lcoord, coord);
            }

            /** \brief this function returns surface normal vector of a manifold
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetSurfaceNormal()
             */
            void GetSurfaceNormal(Array<OneD, NekDouble> &SurfaceNormal,
                                  const int k)
            {
                v_GetSurfaceNormal(SurfaceNormal, k);
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
            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                v_WriteToFile(outfile,format,dumpVar,var);
            }

            void WriteTecplotZone(std::ofstream &outfile);

            void WriteTecplotField(std::ofstream &outfile);

            void ReadFromFile(std::ifstream &in, OutputFormat format, const bool dumpVar = true)
            {
                v_ReadFromFile(in,format,dumpVar);
            }

            inline DNekMatSharedPtr& GetStdMatrix(const StdMatrixKey &mkey)
            {
                return m_stdMatrixManager[mkey];
            }

            inline DNekBlkMatSharedPtr& GetStdStaticCondMatrix(const StdMatrixKey &mkey)
            {
                return m_stdStaticCondMatrixManager[mkey];
            }

            DNekScalMatSharedPtr& GetLocMatrix(const MatrixType mtype,
                                               NekDouble lambdaval = NekConstants::kNekUnsetDouble,
                                               NekDouble tau = NekConstants::kNekUnsetDouble)
            {
                return v_GetLocMatrix(mtype,lambdaval,tau);
            }


            DNekScalMatSharedPtr& GetLocMatrix(const LocalRegions::MatrixKey &mkey)
            {
                return v_GetLocMatrix(mkey);
            }


            const Array<OneD, const NekDouble>& GetPhysNormals(void)
            {
                return v_GetPhysNormals();
            }


            void SetPhysNormals(Array<OneD, const NekDouble> &normal)
            {
                v_SetPhysNormals(normal);
            }

            virtual void SetUpPhysNormals(const boost::shared_ptr<StdExpansion>  &exp2d, const int edge);


            void NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray, bool NegateNorm = false)
            {
                v_NormVectorIProductWRTBase(Fx,Fy,outarray,NegateNorm);
            }

            DNekScalBlkMatSharedPtr& GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
            {
                return v_GetLocStaticCondMatrix(mkey);
            }

            StdRegions::FaceOrientation GetFaceorient(int face)
            {
                return v_GetFaceorient(face);
            }

            StdRegions::EdgeOrientation GetEorient(int edge)
            {
                return v_GetEorient(edge);
            }


            StdRegions::EdgeOrientation GetCartesianEorient(int edge)
            {
                return v_GetCartesianEorient(edge);
            }

            void AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
            {
                v_AddHDGHelmholtzTraceTerms(tau,inarray, outarray);
            }

            void AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,boost::shared_ptr<StdExpansion1D> > &EdgeExp,
                                           Array<OneD,NekDouble> &outarray)
            {
                v_AddHDGHelmholtzTraceTerms(tau,inarray,EdgeExp, outarray);
            }

            // virtual functions related to LocalRegions

            virtual void AddEdgeNormBoundaryInt(const int edge,
                                                boost::shared_ptr<StdExpansion1D>  &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fx,
                                                const Array<OneD, const NekDouble> &Fy,
                                                Array<OneD, NekDouble> &outarray);

            virtual void AddEdgeNormBoundaryInt(const int edge,
                                                boost::shared_ptr<StdExpansion1D>  &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fn,
                                                Array<OneD, NekDouble> &outarray);

            virtual void AddEdgeNormBoundaryBiInt(const int edge,
                                                boost::shared_ptr<StdExpansion1D>  &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fwd,
                                                const Array<OneD, const NekDouble> &Bwd,
                                                Array<OneD, NekDouble> &outarray);

            virtual void AddNormTraceInt(const int dir,
                                         Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray);


            int GetCoordim()
            {
                return v_GetCoordim();
            }

            void GetBoundaryMap(Array<OneD, unsigned int> &outarray)
            {
                v_GetBoundaryMap(outarray);
            }

            void GetInteriorMap(Array<OneD, unsigned int> &outarray)
            {
                v_GetInteriorMap(outarray);
            }

            int GetVertexMap(const int localVertexId)
            {
                return v_GetVertexMap(localVertexId);
            }

            void GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                    Array<OneD, unsigned int> &maparray,
                                    Array<OneD, int> &signarray)
            {
                v_GetEdgeInteriorMap(eid,edgeOrient,maparray,signarray);
            }

            void GetFaceInteriorMap(const int fid, const FaceOrientation faceOrient,
                                    Array<OneD, unsigned int> &maparray,
                                    Array<OneD, int> &signarray)
            {
                v_GetFaceInteriorMap(fid,faceOrient,maparray,signarray);
            }

            void GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray)
            {
                v_GetEdgeToElementMap(eid,edgeOrient,maparray,signarray);
            }

            void GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray)
            {
                v_GetFaceToElementMap(fid,faceOrient,maparray,signarray);
            }

            void GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                v_GetEdgePhysVals(edge,inarray,outarray);
            }

            void GetEdgePhysVals(const int edge, const boost::shared_ptr<StdExpansion1D>   &EdgeExp, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                v_GetEdgePhysVals(edge,EdgeExp,inarray,outarray);
            }

        boost::shared_ptr<StdExpansion1D> GetEdgeExp(int edge, bool SetUpNormals=true)
            {
          return v_GetEdgeExp(edge,SetUpNormals);
            }

            // Matrix Routines

            /** \brief this function generates the mass matrix
             *  \f$\mathbf{M}[i][j] =
             *  \int \phi_i(\mathbf{x}) \phi_j(\mathbf{x}) d\mathbf{x}\f$
             *
             *  \return returns the mass matrix
             */

            DNekMatSharedPtr CreateGeneralMatrix(const StdMatrixKey &mkey);

            void GeneralMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray,
                                 const StdMatrixKey &mkey);

            void MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                              Array<OneD,NekDouble> &outarray,
                              const StdMatrixKey &mkey)
            {
                v_MassMatrixOp(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                v_LaplacianMatrixOp(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const int k1, const int k2,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                v_LaplacianMatrixOp(k1,k2,inarray,outarray,mkey);
            }

            void WeakDerivMatrixOp(const int i,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                v_WeakDerivMatrixOp(i,inarray,outarray,mkey);
            }

            void WeakDirectionalDerivMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray,
                                              const StdMatrixKey &mkey)
            {
                v_WeakDirectionalDerivMatrixOp(inarray,outarray,mkey);
            }

            void MassLevelCurvatureMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray,
                                            const StdMatrixKey &mkey)
            {
                v_MassLevelCurvatureMatrixOp(inarray,outarray,mkey);
            }

            void LinearAdvectionDiffusionReactionMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                          const StdMatrixKey &mkey,
                                                          bool addDiffusionTerm = true)
            {
                v_LinearAdvectionDiffusionReactionMatrixOp(inarray,outarray,mkey,addDiffusionTerm);
            }

            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                v_HelmholtzMatrixOp(inarray,outarray,mkey);
            }

            DNekMatSharedPtr GenMatrix (const StdMatrixKey &mkey)
            {
                return v_GenMatrix(mkey);
            }

            void PhysDeriv (const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                            Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                v_PhysDeriv (inarray, out_d0, out_d1, out_d2);
            }

            void PhysDeriv(const int dir,
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray)
            {
                v_PhysDeriv (dir, inarray, outarray);
            }

            void PhysDirectionalDeriv(const Array<OneD, const NekDouble>& inarray,
                                      const Array<OneD, const NekDouble>& direction,
                                      Array<OneD, NekDouble> &outarray)
            {
                v_PhysDirectionalDeriv (inarray, direction, outarray);
            }

            void StdPhysDeriv(const Array<OneD, const NekDouble>& inarray,
                              Array<OneD, NekDouble> &out_d0,
                              Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                              Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                v_StdPhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            void StdPhysDeriv (const int dir,
                               const Array<OneD, const NekDouble>& inarray,
                               Array<OneD, NekDouble> &outarray)
            {
                v_StdPhysDeriv(dir,inarray,outarray);
            }

            void AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
            {
                v_AddRobinMassMatrix(edgeid,primCoeffs,inoutmat);
            }

            void DGDeriv(const int dir,
                         const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, boost::shared_ptr< StdExpansion1D > > &EdgeExp,
                         Array<OneD, NekDouble> &outarray)
            {
                v_DGDeriv (dir, inarray, EdgeExp, outarray);
            }


            /** \brief This function evaluates the expansion at a single
             *  (arbitrary) point of the domain
             *
             *  This function is a wrapper around the virtual function
             *  \a v_PhysEvaluate()
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
            NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coords)
            {
                return v_PhysEvaluate(coords);
            }

            const boost::shared_ptr<SpatialDomains::GeomFactors>& GetMetricInfo(void) const
            {
                return v_GetMetricInfo();
            }

            const boost::shared_ptr<SpatialDomains::Geometry> GetGeom(void) const
            {
                return v_GetGeom();
            }

            const boost::shared_ptr<SpatialDomains::Geometry1D>& GetGeom1D(void) const
            {
                return v_GetGeom1D();
            }

            const boost::shared_ptr<SpatialDomains::Geometry2D>& GetGeom2D(void) const
            {
                return v_GetGeom2D();
            }

            const boost::shared_ptr<SpatialDomains::Geometry3D>& GetGeom3D(void) const
            {
                return v_GetGeom3D();
            }

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const LocalRegions::MatrixKey &mkey);

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype, NekDouble lambdaval, NekDouble tau);

            virtual const Array<OneD, const NekDouble>& v_GetPhysNormals(void);

            virtual void v_SetPhysNormals(Array<OneD, const NekDouble> &normal);

            virtual void v_SetUpPhysNormals(const boost::shared_ptr<StdExpansion> &exp2d, const int edge);

            virtual void v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray, bool NegateNorm = false);

            virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey);

            virtual StdRegions::FaceOrientation v_GetFaceorient(int face);

            virtual StdRegions::EdgeOrientation v_GetEorient(int edge);

            virtual StdRegions::EdgeOrientation v_GetCartesianEorient(int edge);

            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray);

            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD, boost::shared_ptr< StdExpansion1D > > &edgeExp,
                                                     Array<OneD,NekDouble> &outarray);

            virtual void v_AddEdgeNormBoundaryInt(const int edge,
                                                  boost::shared_ptr<StdExpansion1D> &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fx,
                                                  const Array<OneD, const NekDouble> &Fy,
                                                  Array<OneD, NekDouble> &outarray);

             virtual void v_AddEdgeNormBoundaryInt(const int edge,
                                                  boost::shared_ptr<StdExpansion1D> &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fn,
                                                  Array<OneD, NekDouble> &outarray);

            virtual void v_AddEdgeNormBoundaryBiInt(const int edge,
                                                    boost::shared_ptr<StdExpansion1D> &EdgeExp,
                                                    const Array<OneD, const NekDouble> &Fwd,
                                                    const Array<OneD, const NekDouble> &Bwd,
                                                    Array<OneD, NekDouble> &outarray);

           virtual void v_AddNormTraceInt(const int dir,
                                           Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray);

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
            NekDouble Linf(const Array<OneD, const NekDouble>& sol);

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
            NekDouble L2(const Array<OneD, const NekDouble>& sol);

            /** \brief Function to evaluate the discrete \f$ L_2\f$ norm of the
             *  function defined at the physical points \a (this)->m_phys.
             *
             *    This function takes the physical value space array \a m_phys as
             *  discrete function to be evaluated
             *
             *  \return returns the \f$ L_2 \f$ norm as a double
             */
            NekDouble L2();

            /** \brief Function to evaluate the discrete \f$ H^1\f$
             *  error, \f$ | \epsilon |^1_{2} = \left [ \int^1_{-1} [u -
             *  u_{exact}]^2 + \nabla(u - u_{exact})\cdot\nabla(u -
             *  u_{exact})\cdot dx \right]^{1/2} d\xi_1 \f$ where \f$
             *  u_{exact}\f$ is given by the array \a sol.
             *
             *    This function takes the physical value space array
             *  \a m_phys as approximate solution
             *
             *  \param sol array of solution function  at physical quadrature
             *  points
             *  \return returns the \f$ H_1 \f$ error as a double.
             */
            NekDouble H1(const Array<OneD, const NekDouble>& sol);

            /** \brief Function to evaluate the discrete \f$ H^1\f$ norm of the
             *  function defined at the physical points \a (this)->m_phys.
             *
             *    This function takes the physical value space array
             *  \a m_phys as discrete function to be evaluated
             *
             *  \return returns the \f$ H^1 \f$ norm as a double
             */
            NekDouble H1();

            // I/O routines
            void WriteCoeffsToFile(std::ofstream &outfile);


        protected:


            int   m_elmt_id;  ///< id of element when used in a list.
            int   m_numbases;                                 /**< Number of 1D basis defined in expansion */
            Array<OneD, LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */
            int  m_ncoeffs;                                   /**< Total number of coefficients used in the expansion */
            Array<OneD, NekDouble> m_coeffs;                  /**< Array containing expansion coefficients */
            Array<OneD, NekDouble> m_phys;                    /**< Array containing expansion evaluated at the quad points */
            LibUtilities::NekManager<StdMatrixKey, DNekMat, StdMatrixKey::opLess> m_stdMatrixManager;
            LibUtilities::NekManager<StdMatrixKey, DNekBlkMat, StdMatrixKey::opLess> m_stdStaticCondMatrixManager;

            bool StdMatManagerAlreadyCreated(const StdMatrixKey &mkey)
            {
                return m_stdMatrixManager.AlreadyCreated(mkey);
            }

            bool StdStaticCondMatManagerAlreadyCreated(const StdMatrixKey &mkey)
            {
                return m_stdStaticCondMatrixManager.AlreadyCreated(mkey);
            }

            DNekMatSharedPtr CreateStdMatrix(const StdMatrixKey &mkey)
            {
                return v_CreateStdMatrix(mkey);
            }

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

            void BwdTrans_MatOp(const Array<OneD, const NekDouble>& inarray,
                                Array<OneD, NekDouble> &outarray);

            void BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                v_BwdTrans_SumFac(inarray,outarray);
            }

            void IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTBase_SumFac(inarray,outarray);
            }

            void IProductWRTDerivBase_SumFac(const int dir,
                                             const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }

            // The term _MatFree denotes that the action of the MatrixOperation
            // is done withouth actually using the matrix (which then needs to be stored/calculated).
            // Although this does not strictly mean that no matrix operations are involved in the
            // evaluation of the operation, we use this term in the same context used as in the following
            // paper:
            // R. C. Kirby, M. G. Knepley, A. Logg, and L. R. Scott,
            // "Optimizing the evaluation of finite element matrices," SISC 27:741-758 (2005)
            void GeneralMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray,
                                               const StdMatrixKey &mkey);

            void MassMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray,
                                            const StdMatrixKey &mkey);

            void LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdMatrixKey &mkey)
            {
                v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp_MatFree_GenericImpl(const Array<OneD, const NekDouble> &inarray,
                                                             Array<OneD,NekDouble> &outarray,
                                                             const StdMatrixKey &mkey);

            void LaplacianMatrixOp_MatFree(const int k1, const int k2,
                                                 const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdMatrixKey &mkey);

            void WeakDerivMatrixOp_MatFree(const int i,
                                                 const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdMatrixKey &mkey);

            void WeakDirectionalDerivMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                      Array<OneD,NekDouble> &outarray,
                                                      const StdMatrixKey &mkey);

            void MassLevelCurvatureMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                    Array<OneD,NekDouble> &outarray,
                                                    const StdMatrixKey &mkey);

            void LinearAdvectionDiffusionReactionMatrixOp_MatFree( const Array<OneD, const NekDouble> &inarray,
                                                                   Array<OneD,NekDouble> &outarray,
                                                                   const StdMatrixKey &mkey,
                                                                   bool addDiffusionTerm = true);

            void HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdMatrixKey &mkey)
            {
                v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
            }

            void HelmholtzMatrixOp_MatFree_GenericImpl(const Array<OneD, const NekDouble> &inarray,
                                                             Array<OneD,NekDouble> &outarray,
                                                             const StdMatrixKey &mkey);

        private:
            // Virtual functions
            virtual int v_GetNverts() const = 0;
            virtual int v_GetNedges() const;

            virtual int v_GetNfaces() const;

            virtual int v_NumBndryCoeffs() const;

            virtual int v_NumDGBndryCoeffs() const;

            virtual int v_GetEdgeNcoeffs(const int i) const;

            virtual int v_GetEdgeNumPoints(const int i) const;

            virtual int v_DetCartesianDirOfEdge(const int edge);

            virtual const LibUtilities::BasisKey v_DetEdgeBasisKey(const int i) const;

            virtual int v_GetFaceNcoeffs(const int i) const;

            virtual int v_GetFaceIntNcoeffs(const int i) const;

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const;

            virtual ExpansionType v_DetExpansionType() const;

            virtual int v_GetShapeDimension() const;

            virtual bool  v_IsBoundaryInteriorExpansion();

            virtual void   v_BwdTrans   (const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray) = 0;
            virtual void   v_FwdTrans   (const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray) = 0;
            virtual void  v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray) = 0;

            virtual void  v_IProductWRTDerivBase (const int dir,
                                                   const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray);

            virtual void v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray);

            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray );

            virtual void   v_PhysDeriv (const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2,
                                        Array<OneD, NekDouble> &out_d3);

            virtual void v_PhysDeriv(const int dir,
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_d0);

            virtual void v_PhysDirectionalDeriv(const Array<OneD, const NekDouble>& inarray,
                                                const Array<OneD, const NekDouble>& direction,
                                                Array<OneD, NekDouble> &outarray);

            virtual void v_StdPhysDeriv (const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &out_d1,
                                         Array<OneD, NekDouble> &out_d2,
                                         Array<OneD, NekDouble> &out_d3);

            virtual void   v_StdPhysDeriv (const int dir,
                                           const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray);

            virtual void v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat);

            virtual void v_DGDeriv(const int dir,
                                   const Array<OneD, const NekDouble>& inarray,
                                   Array<OneD, boost::shared_ptr< StdExpansion1D > > &EdgeExp,
                                   Array<OneD, NekDouble> &outarray);

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords);

            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey);

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey);

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                     Array<OneD, NekDouble> &coords_1,
                                     Array<OneD, NekDouble> &coords_2);

            virtual void v_GetCoord(const Array<OneD, const NekDouble>& Lcoord,
                                    Array<OneD, NekDouble> &coord);

            virtual int v_GetCoordim(void);

            virtual void v_GetSurfaceNormal(Array<OneD, NekDouble> &SurfaceNormal,
                                            const int k);

            virtual void v_GetBoundaryMap(Array<OneD, unsigned int>& outarray);

            virtual void v_GetInteriorMap(Array<OneD, unsigned int>& outarray);

            virtual int v_GetVertexMap(int localVertexId);

            virtual void v_GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray);

            virtual void v_GetFaceInteriorMap(const int fid, const FaceOrientation faceOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray);

            virtual void v_GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray);

            virtual void v_GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray);

            virtual void v_GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray);

            virtual void v_GetEdgePhysVals(const int edge,  const boost::shared_ptr<StdExpansion1D>  &EdgeExp, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray);

            virtual boost::shared_ptr<StdExpansion1D> v_GetEdgeExp(const int edge, bool SetUpNormals=true);

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");

            virtual void v_ReadFromFile(std::ifstream &infile, OutputFormat format, const bool dumpVar = true);

            virtual const  boost::shared_ptr<SpatialDomains::GeomFactors>& v_GetMetricInfo() const;

            virtual const boost::shared_ptr<SpatialDomains::Geometry> v_GetGeom() const;

            virtual const boost::shared_ptr<SpatialDomains::Geometry1D>& v_GetGeom1D() const;

            virtual const boost::shared_ptr<SpatialDomains::Geometry2D>& v_GetGeom2D() const;

            virtual const boost::shared_ptr<SpatialDomains::Geometry3D>& v_GetGeom3D() const;

            virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray);

            virtual void v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                                  Array<OneD, NekDouble> &outarray);

            virtual void v_IProductWRTDerivBase_SumFac(const int dir,
                                                       const Array<OneD, const NekDouble>& inarray,
                                                       Array<OneD, NekDouble> &outarray);

            virtual void v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray,
                                        const StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp(const int k1, const int k2,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey);

            virtual void v_WeakDerivMatrixOp(const int i,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey);

            virtual void v_WeakDirectionalDerivMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray,
                                                        const StdMatrixKey &mkey);

            virtual void v_MassLevelCurvatureMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray,
                                                        const StdMatrixKey &mkey);

            virtual void v_LinearAdvectionDiffusionReactionMatrixOp(const Array<OneD,
                                                                    const NekDouble> &inarray,
                                                                    Array<OneD,NekDouble> &outarray,
                                                                    const StdMatrixKey &mkey,
                                                                    bool addDiffusionTerm=true);

            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                           Array<OneD,NekDouble> &outarray,
                                                           const StdMatrixKey &mkey);

            virtual void v_HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                           Array<OneD,NekDouble> &outarray,
                                                           const StdMatrixKey &mkey);
#endif

        };


        typedef boost::shared_ptr<StdExpansion> StdExpansionSharedPtr;
        typedef std::vector< StdExpansionSharedPtr > StdExpansionVector;
        typedef std::vector< StdExpansionSharedPtr >::iterator StdExpansionVectorIter;

    } //end of namespace
} //end of namespace

#endif //STANDARDDEXPANSION_H
/**
 * $Log: StdExpansion.h,v $
 * Revision 1.127  2010/01/03 19:39:09  cantwell
 * Added FldToVtk converter.
 * Added XmlToVtk converter.
 *
 * Revision 1.126  2009/12/14 18:02:55  cbiotto
 * Adding functions for tecplot file
 *
 * Revision 1.125  2009/11/10 19:01:37  sehunchun
 * Update related to Variable coefficients of HDG2D Solver
 *
 * Revision 1.124  2009/11/06 21:42:17  sherwin
 * Added call to DGDeriv function
 *
 * Revision 1.123  2009/11/02 19:15:43  cantwell
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
 * Revision 1.122  2009/10/25 18:53:39  sherwin
 * Added H1 norm definition
 *
 * Revision 1.121  2009/09/22 21:23:48  bnelson
 * Fixed windows compile error.
 * Committed on the Free edition of March Hare Software CVSNT Server.
 * Upgrade to CVS Suite for more features and support:
 * http://march-hare.com/cvsnt/
 *
 * Revision 1.120  2009/09/06 21:55:26  sherwin
 * Updates related to Navier Stokes Solver
 *
 * Revision 1.119  2009/08/14 09:24:08  cbiotto
 * Fixed bug in GetEdgeExp
 *
 * Revision 1.118  2009/07/20 14:15:44  cbiotto
 * Adding GetEdgeExp function
 *
 * Revision 1.117  2009/07/08 17:21:12  sehunchun
 * Deleting GetTanBasis
 *
 * Revision 1.116  2009/07/08 11:12:48  sehunchun
 * Adding GetSurfaceNormal function
 *
 * Revision 1.115  2009/07/07 16:33:14  sehunchun
 * Adding AddEdgeNormBoundaryBiInt
 *
 * Revision 1.114  2009/07/03 15:37:12  sehunchun
 * Adding GetTanBasis function
 *
 * Revision 1.113  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.112  2009/04/22 22:30:48  sherwin
 * Added ReadFromFile method to read back in .dat file
 *
 * Revision 1.111  2009/04/20 16:11:47  sherwin
 * Mods to handle output and optimise DG work
 *
 * Revision 1.110  2009/04/03 14:57:34  sherwin
 * Linear Advection matrices added, corrected unsigned int intialisation
 *
 * Revision 1.109  2009/03/04 14:17:38  pvos
 * Removed all methods that take and Expansion as argument
 *
 * Revision 1.108  2008/12/18 14:11:35  pvos
 * NekConstants Update
 *
 * Revision 1.107  2008/11/24 10:31:14  pvos
 * Changed name from _PartitionedOp to _MatFree
 *
 * Revision 1.106  2008/11/19 16:02:47  pvos
 * Added functionality for variable Laplacian coeffcients
 *
 * Revision 1.105  2008/11/05 16:08:15  pvos
 * Added elemental optimisation functionality
 *
 * Revision 1.104  2008/11/01 22:36:06  bnelson
 * Removed uneeded files.
 *
 * Revision 1.103  2008/10/19 15:55:46  sherwin
 * Added methods EvalBasisNumModesMax
 *
 * Revision 1.102  2008/10/04 19:27:41  sherwin
 * Added AddEdgeNormBoundaryInt
 *
 * Revision 1.101  2008/09/23 18:19:26  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.100  2008/09/17 13:46:02  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.99  2008/09/16 13:37:03  pvos
 * Restructured the LocalToGlobalMap classes
 *
 * Revision 1.98  2008/09/15 13:18:08  pvos
 * Added more hexahedron mapping routines
 *
 * Revision 1.97  2008/09/09 14:14:27  sherwin
 * Remove Interp1D/2D/3D and added DetCartesianDirOfEdge
 *
 * Revision 1.96  2008/08/27 16:34:53  pvos
 * Small efficiency update
 *
 * Revision 1.95  2008/08/14 22:09:50  sherwin
 * Modifications to remove HDG routines from StdRegions and removed StdExpMap
 *
 * Revision 1.94  2008/08/03 20:13:03  sherwin
 * Put return values in virtual functions
 *
 * Revision 1.93  2008/07/31 11:10:15  sherwin
 * Updates for handling EdgeBasisKey for use with DG advection. Depracated GetEdgeBasis and added DetEdgeBasisKey
 *
 * Revision 1.92  2008/07/29 22:21:15  sherwin
 * A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
 *
 * Revision 1.91  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.90  2008/07/16 22:20:27  sherwin
 * Added AddEdgeNormBoundaryInt
 *
 * Revision 1.89  2008/07/12 19:08:29  sherwin
 * Modifications for DG advection routines
 *
 * Revision 1.88  2008/07/12 16:30:07  sherwin
 * Added an new member m_elmt_id so that there is an element number for use later in lists
 *
 * Revision 1.87  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.86  2008/07/02 14:08:56  pvos
 * Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 * Revision 1.85  2008/06/16 22:45:15  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.84  2008/06/13 00:27:20  ehan
 * Added GetFaceNCoeffs() function.
 *
 * Revision 1.83  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.82  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.81  2008/05/10 18:27:33  sherwin
 * Modifications necessary for QuadExp Unified DG Solver
 *
 * Revision 1.80  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.79  2008/04/06 06:04:14  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.78  2008/04/02 22:18:10  pvos
 * Update for 2D local to global mapping
 *
 * Revision 1.77  2008/03/18 14:15:45  pvos
 * Update for nodal triangular helmholtz solver
 *
 * Revision 1.76  2008/03/12 15:25:09  pvos
 * Clean up of the code
 *
 * Revision 1.75  2008/02/29 19:15:19  sherwin
 * Update for UDG stuff
 *
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

