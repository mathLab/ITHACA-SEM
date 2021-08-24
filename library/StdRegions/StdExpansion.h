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
#include <vector>
#include <memory>

#include <boost/core/ignore_unused.hpp>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdRegionsDeclspec.h>
#include <StdRegions/StdMatrixKey.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
namespace Nektar { namespace LocalRegions { class MatrixKey; class Expansion; } }

namespace Nektar
{
    namespace StdRegions
    {

        /** \brief The base class for all shapes
         *
         *  This is the lowest level basic class for all shapes and so
         *  contains the definition of common data and common routine to all
         *  elements
         */
        class StdExpansion : public std::enable_shared_from_this<StdExpansion>
        {
        public:

            /** \brief Default Constructor */
            STD_REGIONS_EXPORT StdExpansion();

            /** \brief Constructor */
            STD_REGIONS_EXPORT StdExpansion(const int numcoeffs, const int numbases,
                         const LibUtilities::BasisKey &Ba = LibUtilities::NullBasisKey,
                         const LibUtilities::BasisKey &Bb = LibUtilities::NullBasisKey,
                         const LibUtilities::BasisKey &Bc = LibUtilities::NullBasisKey);


            /** \brief Copy Constructor */
            STD_REGIONS_EXPORT StdExpansion(const StdExpansion &T);

            /** \brief Destructor */
            STD_REGIONS_EXPORT virtual ~StdExpansion();


            // Standard Expansion Routines Applicable Regardless of Region

            /** \brief This function returns the number of 1D bases used in
             *  the expansion
             *
             *  \return returns the number of 1D bases used in the expansion,
             *  which is equal to number dimension of the expansion
             */
            inline int GetNumBases() const
            {
                return m_base.size();
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
                ASSERTL1(dir < m_base.size(),
                         "dir is larger than number of bases");
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

            /** \brief This function returns the total number of quadrature
             *  points used in the element
             *
             *  \return returns the total number of quadrature points
             */
            inline  int GetTotPoints() const
            {
                int i;
                int nqtot = 1;

                for(i=0; i < m_base.size(); ++i)
                {
                    nqtot *= m_base[i]->GetNumPoints();
                }

                return  nqtot;
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
                ASSERTL1(dir < m_base.size(), "dir is larger than m_numbases");
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
                ASSERTL1(dir < m_base.size(),"dir is larger than m_numbases");
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

                for(i = 0; i < m_base.size(); ++i)
                {
                    returnval = std::max(returnval, m_base[i]->GetNumModes());
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
                ASSERTL1(dir < m_base.size(), "dir is larger than m_numbases");
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
                ASSERTL1(dir < m_base.size() || dir == 0,
                         "dir is larger than m_numbases");
                return(m_base.size() > 0 ? m_base[dir]->GetNumPoints() : 1);
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


            /** \brief This function returns the number of expansion coefficients
             *  belonging to the \a i-th trace
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetTraceNcoeffs()
             *
             *  \param i specifies which trace
             *  \return returns the number of expansion coefficients belonging to
             *  the \a i-th trace
             */
            int GetTraceNcoeffs(const int i) const
            {
                return v_GetTraceNcoeffs(i);
            }

            int GetTraceIntNcoeffs(const int i) const
            {
                return v_GetTraceIntNcoeffs(i);
            }

            
            /** \brief This function returns the number of quadrature points
             *  belonging to the \a i-th trace
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetTraceNumPoints()
             *
             *  \param i specifies which trace id
             *  \return returns the number of quadrature points belonging to
             *   the \a i-th trace
             */
            int GetTraceNumPoints(const int i) const
            {
                return v_GetTraceNumPoints(i);
            }

            /** \brief This function returns the basis key belonging
             *   to the \a i-th trace
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetTraceBasisKey()
             *
             *  \param i specifies which trace id
             *  \param k is the direction of the basis key for 2D traces
             *
             *  \return returns the number of Basis key of the ith
             *  trace in the k th direction (when trace is a 2D
             *  object)
             */
            const LibUtilities::BasisKey GetTraceBasisKey(const int i,
                                                          int k = -1) const
            {
                return v_GetTraceBasisKey(i, k);
            }

            /** \brief This function returns the basis key belonging
             *   to the \a i-th trace
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetTracePointsKey()
             *
             *  \param i specifies which trace id
             *  \param k is the direction of the basis key for 2D traces
             *
             *  \return returns the number of Points key of the ith
             *  trace in the k th direction (when trace is a 2D
             *  object)
             */
            LibUtilities::PointsKey GetTracePointsKey(const int i,
                                                      int k = -1) const
            {
                return v_GetTracePointsKey(i, k);
            }


            int NumBndryCoeffs(void)  const
            {
                return v_NumBndryCoeffs();
            }

            int NumDGBndryCoeffs(void)  const
            {
                return v_NumDGBndryCoeffs();
            }


            /** \brief This function returns the type of expansion
             *  Nodal point type if defined
             *
             *  This function is a wrapper around the virtual function
             *  \a v_GetNodalPointsKey()
             *
             */
            const LibUtilities::PointsKey GetNodalPointsKey() const
            {
                return v_GetNodalPointsKey();
            };

            /**
             * @brief Returns the number of trace elements connected to this
             * element.
             *
             * For example, a quadrilateral has four edges, so this function
             * would return 4.
             */
            int GetNtraces() const
            {
                return v_GetNtraces(); 
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
            LibUtilities::ShapeType DetShapeType() const
            {
                return v_DetShapeType();
            }

            std::shared_ptr<StdExpansion> GetStdExp(void) const
            {
                return v_GetStdExp();
            }

            std::shared_ptr<StdExpansion> GetLinStdExp(void) const
            {
                return v_GetLinStdExp();
            }

            int GetShapeDimension() const
            {
                return v_GetShapeDimension();
            }

            bool IsBoundaryInteriorExpansion()
            {
                return v_IsBoundaryInteriorExpansion();
            }

            bool IsNodalNonTensorialExp()
            {
                return v_IsNodalNonTensorialExp();
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
             *  in the \a outarray.
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

            /**
             * @brief This function performs the Forward transformation from
             * physical space to coefficient space.
             */
            inline void FwdTrans (const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray);

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
             *
             *              Inputs:\n

            - \a inarray: definition of function to be returned at quadrature point
            of expansion.

            Outputs:\n

            - returns \f$\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2) J[i,j] d
            \xi_1 d \xi_2 \f$ where \f$inarray[i,j] =
            u(\xi_{1i},\xi_{2j}) \f$ and \f$ J[i,j] \f$ is the
            Jacobian evaluated at the quadrature point.
             *
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
             *            \f$ \begin{array}{rcl} I_{pq} = (\phi_q \phi_q, u) & = &
            \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \phi_p(\xi_{0,i})
            \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j})
            J_{i,j}\\ & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
            \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
            J_{i,j} \end{array} \f$

            where

            \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$

            which can be implemented as

            \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} =
            {\bf B_1 U}  \f$
            \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} =
            {\bf B_0 F}  \f$
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

            void IProductWRTBase(
                    const Array<OneD, const NekDouble>& base,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray,
                    int coll_check)
            {
                v_IProductWRTBase(base, inarray, outarray, coll_check);
            }


            void   IProductWRTDerivBase(
                    const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDerivBase(dir,inarray, outarray);
            }

            void   IProductWRTDirectionalDerivBase(
                    const Array<OneD, const NekDouble>& direction,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDirectionalDerivBase(direction, inarray, outarray);
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

            inline DNekMatSharedPtr GetStdMatrix(const StdMatrixKey &mkey)
            {
                return m_stdMatrixManager[mkey];
            }

            inline DNekBlkMatSharedPtr GetStdStaticCondMatrix(const StdMatrixKey &mkey)
            {
                return m_stdStaticCondMatrixManager[mkey];
            }

            void NormVectorIProductWRTBase
                (const Array<OneD, const NekDouble> &Fx,
                 Array< OneD, NekDouble> &outarray)
            {
                v_NormVectorIProductWRTBase(Fx,outarray);
            }

            void NormVectorIProductWRTBase
                (const Array<OneD, const NekDouble> &Fx,
                 const Array<OneD, NekDouble> &Fy,
                 Array< OneD, NekDouble> &outarray)
            {
                v_NormVectorIProductWRTBase(Fx,Fy,outarray);
            }

            void NormVectorIProductWRTBase
                 (const Array<OneD, const NekDouble> &Fx,
                  const Array<OneD, const NekDouble> &Fy,
                  const Array<OneD, const NekDouble> &Fz,
                  Array< OneD, NekDouble> &outarray)
            {
                v_NormVectorIProductWRTBase(Fx,Fy,Fz,outarray);
            }

            void NormVectorIProductWRTBase
                 (const Array<OneD,
                  const Array<OneD, NekDouble> > &Fvec,
                  Array< OneD, NekDouble> &outarray)
            {
                v_NormVectorIProductWRTBase(Fvec, outarray);
            }

            DNekScalBlkMatSharedPtr GetLocStaticCondMatrix
                (const LocalRegions::MatrixKey &mkey)
            {
                return v_GetLocStaticCondMatrix(mkey);
            }

            STD_REGIONS_EXPORT void DropLocStaticCondMatrix
                 (const LocalRegions::MatrixKey &mkey)
            {
                return v_DropLocStaticCondMatrix(mkey);
            }

            int CalcNumberOfCoefficients(const std::vector<unsigned int>  &nummodes,
                                         int &modes_offset)
            {
                return v_CalcNumberOfCoefficients(nummodes,modes_offset);
            }

            // virtual functions related to LocalRegions
            STD_REGIONS_EXPORT NekDouble StdPhysEvaluate(
                                            const Array<OneD, const NekDouble> &Lcoord,
                                            const Array<OneD, const NekDouble> &physvals);

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

            int GetVertexMap(const int localVertexId,
                             bool useCoeffPacking = false)
            {
                return v_GetVertexMap(localVertexId,useCoeffPacking);
            }

            void GetTraceToElementMap(
                    const int                  tid,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD,          int> &signarray,
                    Orientation                traceOrient = eForwards,
                    int                        P = -1,
                    int                        Q = -1)
            {
                v_GetTraceToElementMap(tid,maparray,signarray,traceOrient,P,Q);
            }

            void GetTraceInteriorToElementMap(
                    const int                  tid,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD,          int> &signarray,
                    const Orientation          traceOrient = eForwards)
            {
                v_GetTraceInteriorToElementMap(tid,maparray,signarray,traceOrient);
            }
            

            void GetTraceNumModes(const int tid, 
                                 int &numModes0,
                                 int &numModes1,
                                 const Orientation traceOrient
                                                     = eDir1FwdDir1_Dir2FwdDir2)
            {
                v_GetTraceNumModes(tid,numModes0,numModes1,traceOrient);
            }

            void MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                v_MultiplyByQuadratureMetric(inarray, outarray);
            }

            void MultiplyByStdQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD, NekDouble> & outarray)
            {
                v_MultiplyByStdQuadratureMetric(inarray, outarray);
            }

            // Matrix Routines

            /** \brief this function generates the mass matrix
             *  \f$\mathbf{M}[i][j] =
             *  \int \phi_i(\mathbf{x}) \phi_j(\mathbf{x}) d\mathbf{x}\f$
             *
             *  \return returns the mass matrix
             */

            STD_REGIONS_EXPORT DNekMatSharedPtr CreateGeneralMatrix(const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void GeneralMatrixOp(const Array<OneD, const NekDouble> &inarray,
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

            void ReduceOrderCoeffs(int numMin,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray)
            {
                v_ReduceOrderCoeffs(numMin,inarray,outarray);
            }

            void SVVLaplacianFilter(Array<OneD,NekDouble> &array,
                                    const StdMatrixKey &mkey)
            {
                v_SVVLaplacianFilter(array,mkey);
            }

            void ExponentialFilter( Array<OneD, NekDouble> &array,
                                    const NekDouble        alpha,
                                    const NekDouble        exponent,
                                    const NekDouble        cutoff)
            {
                v_ExponentialFilter(array, alpha, exponent, cutoff);
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

            void LinearAdvectionDiffusionReactionMatrixOp(
                                        const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray,
                                        const StdMatrixKey &mkey,
                                        bool addDiffusionTerm = true)
            {
                v_LinearAdvectionDiffusionReactionMatrixOp
                    (inarray,outarray,mkey,addDiffusionTerm);
            }

            /**
             * @param   inarray     Input array @f$ \mathbf{u} @f$.
             * @param   outarray    Output array @f$ \boldsymbol{\nabla^2u}
             *                          + \lambda \boldsymbol{u} @f$.
             * @param   mkey
             */
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

            void PhysDeriv_s(const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &out_ds)
            {
                v_PhysDeriv_s(inarray,out_ds);
            }

            void PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
            	             Array<OneD, NekDouble>& out_dn)
            {
            	 v_PhysDeriv_n(inarray,out_dn);
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

            /** \brief This function evaluates the expansion at a single
             *  (arbitrary) point of the domain
             *
             *  This function is a wrapper around the virtual function
             *  \a v_PhysEvaluate()
             *
             *  Based on the value of the expansion at the quadrature
             *  points provided in \a physvals, this function
             *  calculates the value of the expansion at an arbitrary
             *  single points (with coordinates \f$ \mathbf{x_c}\f$
             *  given by the pointer \a coords). This operation,
             *  equivalent to \f[ u(\mathbf{x_c}) = \sum_p
             *  \phi_p(\mathbf{x_c}) \hat{u}_p \f] is evaluated using
             *  Lagrangian interpolants through the quadrature points:
             *  \f[ u(\mathbf{x_c}) = \sum_p h_p(\mathbf{x_c}) u_p\f]
             *
             *  \param coords the coordinates of the single point
             *  \param physvals the interpolated field at the quadrature points
             *
             *  \return returns the value of the expansion at the
             *  single point
             */
            NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coords,
                                   const Array<OneD, const NekDouble>& physvals)
            {
                return v_PhysEvaluate(coords,physvals);
            }


            /** \brief This function evaluates the expansion at a single
             *  (arbitrary) point of the domain
             *
             *  This function is a wrapper around the virtual function
             *  \a v_PhysEvaluate()
             *
             *  Based on the value of the expansion at the quadrature
             *  points provided in \a physvals, this function
             *  calculates the value of the expansion at an arbitrary
             *  single points associated with the interpolation
             *  matrices provided in \f$ I \f$.
             *
             *  \param I an Array of lagrange interpolantes evaluated
             *  at the coordinate and going through the local physical
             *  quadrature
             *  \param physvals the interpolated field at the quadrature points
             *
             *  \return returns the value of the expansion at the
             *  single point
             */
            NekDouble PhysEvaluate(const Array<OneD, DNekMatSharedPtr>& I,
                                   const Array<OneD, const NekDouble >& physvals)
            {
                return v_PhysEvaluate(I,physvals);
            }

            /**
             * @brief This function evaluates the basis function mode @p mode at a
             * point @p coords of the domain.
             *
             * This function uses barycentric interpolation with the tensor
             * product separation of the basis function to improve performance.
             *
             * @param coord   The coordinate inside the standard region.
             * @param mode    The mode number to be evaluated.
             *
             * @return The value of the basis function @p mode at @p coords.
             */
            NekDouble PhysEvaluateBasis(
                const Array<OneD, const NekDouble>& coords,
                int mode)
            {
                return v_PhysEvaluateBasis(coords, mode);
            }

            /**
             * \brief Convert local cartesian coordinate \a xi into local
             * collapsed coordinates \a eta
             **/
            void LocCoordToLocCollapsed(const Array<OneD, const NekDouble>& xi,
                                        Array<OneD, NekDouble>& eta)
            {
                v_LocCoordToLocCollapsed(xi,eta);
            }

            /**
             * \brief Convert local collapsed coordinates \a eta into local
             * cartesian coordinate \a xi
             **/
            void LocCollapsedToLocCoord(const Array<OneD, const NekDouble>& eta,
                                        Array<OneD, NekDouble>& xi)
            {
                v_LocCollapsedToLocCoord(eta,xi);
            }

            STD_REGIONS_EXPORT virtual int v_CalcNumberOfCoefficients
              (const std::vector<unsigned int>  &nummodes, int &modes_offset);

            STD_REGIONS_EXPORT virtual void v_NormVectorIProductWRTBase
            (const Array<OneD, const NekDouble> &Fx, Array< OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_NormVectorIProductWRTBase(
                     const Array<OneD, const NekDouble> &Fx,
                     const Array<OneD, const NekDouble> &Fy,
                     Array< OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_NormVectorIProductWRTBase(
                     const Array<OneD, const NekDouble> &Fx,
                     const Array<OneD, const NekDouble> &Fy,
                     const Array<OneD, const NekDouble> &Fz,
                     Array< OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_NormVectorIProductWRTBase(
                     const Array<OneD, const Array<OneD, NekDouble> > &Fvec,
                     Array< OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual DNekScalBlkMatSharedPtr
              v_GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void
              v_DropLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey);

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
             STD_REGIONS_EXPORT NekDouble Linf(
                                      const Array<OneD,
                                      const NekDouble>& phys,
                                      const Array<OneD,
                                      const NekDouble>& sol = NullNekDouble1DArray);
            
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
             STD_REGIONS_EXPORT NekDouble L2(
                                  const Array<OneD, const NekDouble>& phys,
                                  const Array<OneD, const NekDouble>& sol =
                                  NullNekDouble1DArray);

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
             STD_REGIONS_EXPORT NekDouble H1(
                        const Array<OneD,
                        const NekDouble>& phys,
                        const Array<OneD, const NekDouble>& sol = NullNekDouble1DArray);

            // I/O routines
            const LibUtilities::PointsKeyVector GetPointsKeys() const
            {
                LibUtilities::PointsKeyVector p;
                p.reserve(m_base.size());
                for (int i = 0; i < m_base.size(); ++i)
                {
                    p.push_back(m_base[i]->GetPointsKey());
                }
                return p;
            }

            STD_REGIONS_EXPORT DNekMatSharedPtr BuildInverseTransformationMatrix(
                const DNekScalMatSharedPtr & m_transformationmatrix)
            {
                return v_BuildInverseTransformationMatrix(
                    m_transformationmatrix);
            }


            /** \brief This function performs an interpolation from
             * the physical space points provided at input into an
             * array of equispaced points which are not the collapsed
             * coordinate. So for a tetrahedron you will only get a
             * tetrahedral number of values.
             *
             * This is primarily used for output purposes to get a
             * better distribution of points more suitable for most
             * postprocessing
             */
            STD_REGIONS_EXPORT void PhysInterpToSimplexEquiSpaced(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble>       &outarray,
                int npset = -1);


            /** \brief This function provides the connectivity of
             *   local simplices (triangles or tets) to connect the
             *   equispaced data points provided by
             *   PhysInterpToSimplexEquiSpaced
             *
             *  This is a virtual call to the function
             *  \a v_GetSimplexEquiSpaceConnectivity
             */
            STD_REGIONS_EXPORT void GetSimplexEquiSpacedConnectivity(
                Array<OneD, int> &conn,
                bool              standard = true)
            {
                v_GetSimplexEquiSpacedConnectivity(conn,standard);
            }

            /** \brief This function performs a
             * projection/interpolation from the equispaced points
             * sometimes used in post-processing onto the coefficient
             * space
             *
             * This is primarily used for output purposes to use a
             * more even distribution of points more suitable for alot of
             * postprocessing
             */
             STD_REGIONS_EXPORT void EquiSpacedToCoeffs(
                           const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble>       &outarray);

            template<class T>
            std::shared_ptr<T> as()
            {
                return std::dynamic_pointer_cast<T>( shared_from_this() );
            }

            void IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &outarray,
                                        bool multiplybyweights = true)
            {
                v_IProductWRTBase_SumFac(inarray,outarray,multiplybyweights);
            }

            STD_REGIONS_EXPORT void GenStdMatBwdDeriv(
                const int dir,
                DNekMatSharedPtr &mat)
            {
                v_GenStdMatBwdDeriv(dir,mat);
            }

        protected:
            Array<OneD, LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */
            int m_elmt_id;
            int m_ncoeffs;                                   /**< Total number of coefficients used in the expansion */

            LibUtilities::NekManager<StdMatrixKey, DNekMat, StdMatrixKey::opLess>
            m_stdMatrixManager;
            LibUtilities::NekManager<StdMatrixKey, DNekBlkMat, StdMatrixKey::opLess>
            m_stdStaticCondMatrixManager;

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
            STD_REGIONS_EXPORT DNekBlkMatSharedPtr CreateStdStaticCondMatrix
                         (const StdMatrixKey &mkey);


            STD_REGIONS_EXPORT void BwdTrans_MatOp
                      (const Array<OneD, const NekDouble>& inarray,
                       Array<OneD, NekDouble> &outarray);

            void BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                v_BwdTrans_SumFac(inarray,outarray);
            }

            void IProductWRTDerivBase_SumFac(
                    const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }


            void IProductWRTDirectionalDerivBase_SumFac(
                    const Array<OneD, const NekDouble>& direction,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTDirectionalDerivBase_SumFac(direction,
                                                         inarray,outarray);
            }


            // The term _MatFree denotes that the action of the
            // MatrixOperation is done withouth actually using the
            // matrix (which then needs to be stored/calculated).
            // Although this does not strictly mean that no matrix
            // operations are involved in the evaluation of the
            // operation, we use this term in the same context used as
            // in the following paper: R. C. Kirby, M. G. Knepley,
            // A. Logg, and L. R. Scott, "Optimizing the evaluation of
            // finite element matrices," SISC 27:741-758 (2005)
            STD_REGIONS_EXPORT void GeneralMatrixOp_MatFree
                                (const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray,
                                 const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void MassMatrixOp_MatFree
                                 (const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray,
                                  const StdMatrixKey &mkey);

            void LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdMatrixKey &mkey)
            {
                v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
            }

            STD_REGIONS_EXPORT void LaplacianMatrixOp_MatFree_Kernel(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                      Array<OneD,       NekDouble> &wsp)
            {
                v_LaplacianMatrixOp_MatFree_Kernel(inarray, outarray, wsp);
            }

            STD_REGIONS_EXPORT void LaplacianMatrixOp_MatFree_GenericImpl
                       (const Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray,
                        const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void LaplacianMatrixOp_MatFree
                        (const int k1, const int k2,
                         const Array<OneD, const NekDouble> &inarray,
                         Array<OneD,NekDouble> &outarray,
                         const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void WeakDerivMatrixOp_MatFree(const int i,
                                      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray,
                                      const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void WeakDirectionalDerivMatrixOp_MatFree
                                   (const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,NekDouble> &outarray,
                                    const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void MassLevelCurvatureMatrixOp_MatFree
                                    (const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray,
                                     const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT void LinearAdvectionDiffusionReactionMatrixOp_MatFree
                                    (const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray,
                                     const StdMatrixKey &mkey,
                                     bool addDiffusionTerm = true);

            void HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           const StdMatrixKey &mkey)
            {
                v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
            }

            STD_REGIONS_EXPORT void HelmholtzMatrixOp_MatFree_GenericImpl
                                      (const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,NekDouble> &outarray,
                                       const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_SetCoeffsToOrientation
                                      (StdRegions::Orientation dir,
                                       Array<OneD, const NekDouble> &inarray,
                                       Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_SetCoeffsToOrientation
                                       (Array<OneD, NekDouble> &coeffs,
                                        StdRegions::Orientation dir);
            
            STD_REGIONS_EXPORT virtual NekDouble v_StdPhysEvaluate(
                                const Array<OneD, const NekDouble> &Lcoord,
                                const Array<OneD, const NekDouble> &physvals);

            STD_REGIONS_EXPORT virtual void v_GenStdMatBwdDeriv(
                  const int dir,
                  DNekMatSharedPtr &mat)
            {
                boost::ignore_unused(dir,mat);
                NEKERROR(ErrorUtil::efatal,"not defined");
            }

            STD_REGIONS_EXPORT virtual void v_MultiplyByStdQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            /**
             * @brief This function performs the barycentric interpolation of
             * the polynomial stored in @p coord at a point @p physvals using
             * barycentric interpolation weights in direction @tparam DIR.
             *
             * This method is intended to be used a helper function for
             * StdExpansion::PhysEvaluate and its elemental instances, so that
             * the calling method should provide @p coord for x, y and z
             * sequentially and the appropriate @p physvals and @p weights for
             * that particular direction.
             *
             * @param  coord    The coordinate of the single point.
             * @param  physvals The polynomial stored at each quadrature point.
             * @tparam DIR      The direction of evaluation.
             *
             * @return The value of @p physvals at @p coord in direction @p dir.
             */
            template<int DIR>
            inline NekDouble BaryEvaluate(
                const NekDouble &coord,
                const NekDouble *physvals)
            {
                NekDouble numer = 0.0, denom = 0.0;

                ASSERTL2(DIR < m_base.size(),
                         "Direction should be less than shape dimension.");

                const Array<OneD, const NekDouble> &z = m_base[DIR]->GetZ();
                const Array<OneD, const NekDouble> &bw =
                    m_base[DIR]->GetBaryWeights();

                const auto nquad = z.size();

                for (int i = 0; i < nquad; ++i)
                {
                    NekDouble xdiff = z[i] - coord;
                    NekDouble pval = physvals[i];

                    /*
                     * (in this specific case) you actually 
                     * want to do the comparison exactly 
                     * (believe it or not!) See chapter 7 of 
                     * the paper here:
                     *https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
                     */
                    if (xdiff == 0.0)
                    {
                        return pval;
                    }

                    NekDouble tmp = bw[i] / xdiff;
                    numer += tmp * pval;
                    denom += tmp;
                }

                return numer / denom;
            }

            /**
             * @brief This function evaluates the basis function mode @p mode at
             * a point @p coords of the domain in direction @dir.
             *
             * @param coord   The coordinate inside the standard region.
             * @param mode    The mode number to be evaluated of #m_base[dir]
             * @param dir     The direction of interpolation.
             *
             * @return The value of the basis function @p mode at @p coords in
             *         direction @p dir.
             */
            template<int DIR>
            inline NekDouble BaryEvaluateBasis(const NekDouble &coord,
                                               const int &mode)
            {
                const int nquad = m_base[DIR]->GetNumPoints();
                return BaryEvaluate<DIR>(
                    coord, &(m_base[DIR]->GetBdata())[0] + nquad * mode);
            }

        private:
            // Virtual functions
            STD_REGIONS_EXPORT virtual int v_GetNverts() const = 0;
            STD_REGIONS_EXPORT virtual int v_GetNtraces() const;

            STD_REGIONS_EXPORT virtual int v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int v_NumDGBndryCoeffs() const;

            STD_REGIONS_EXPORT virtual int v_GetTraceNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetTotalTraceIntNcoeffs() const;
            STD_REGIONS_EXPORT virtual int v_GetTraceIntNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetTraceNumPoints(const int i) const;


            STD_REGIONS_EXPORT virtual const LibUtilities::BasisKey
                v_GetTraceBasisKey(const int i, const int k) const;

            STD_REGIONS_EXPORT virtual LibUtilities::PointsKey
                v_GetTracePointsKey(const int i, const int j) const;

            STD_REGIONS_EXPORT virtual const LibUtilities::PointsKey
                v_GetNodalPointsKey() const;

            STD_REGIONS_EXPORT virtual LibUtilities::ShapeType
                v_DetShapeType() const;

            STD_REGIONS_EXPORT virtual std::shared_ptr<StdExpansion>
                v_GetStdExp(void) const;

            STD_REGIONS_EXPORT virtual std::shared_ptr<StdExpansion>
                v_GetLinStdExp(void) const;

            STD_REGIONS_EXPORT virtual int
                v_GetShapeDimension() const;

            STD_REGIONS_EXPORT virtual bool
                v_IsBoundaryInteriorExpansion();

            STD_REGIONS_EXPORT virtual bool
                v_IsNodalNonTensorialExp();

            STD_REGIONS_EXPORT virtual void   v_BwdTrans
                  (const Array<OneD, const NekDouble>& inarray,
                   Array<OneD, NekDouble> &outarray) = 0;

            /**
             * @brief Transform a given function from physical quadrature space
             * to coefficient space.
             * @see StdExpansion::FwdTrans
             */
            STD_REGIONS_EXPORT virtual void   v_FwdTrans   (
                            const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD,       NekDouble> &outarray) = 0;

            /**
             * @brief Calculates the inner product of a given function \a f
             * with the different modes of the expansion
             */
            STD_REGIONS_EXPORT virtual void  v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD,       NekDouble> &outarray) = 0;

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& base,
                            const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD,       NekDouble>& outarray,
                                  int coll_check)
            {
                boost::ignore_unused(base, inarray, outarray, coll_check);
                NEKERROR(ErrorUtil::efatal,
                         "StdExpansion::v_IProductWRTBase has no "
                         "(and should have no) implementation");
            }

            STD_REGIONS_EXPORT virtual void  v_IProductWRTDerivBase(
             const int dir,
             const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &outarray);
    
            STD_REGIONS_EXPORT virtual void  v_IProductWRTDirectionalDerivBase(
             const Array<OneD, const NekDouble>& direction,
             const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_FwdTrans_BndConstrained(
             const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
             const Array<OneD, const NekDouble>& inarray );

            STD_REGIONS_EXPORT virtual void   v_PhysDeriv
            (const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &out_d1,
             Array<OneD, NekDouble> &out_d2,
             Array<OneD, NekDouble> &out_d3);

	    STD_REGIONS_EXPORT virtual void v_PhysDeriv_s
            (const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &out_ds);

            STD_REGIONS_EXPORT virtual void v_PhysDeriv_n
            (const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble>& out_dn);

            STD_REGIONS_EXPORT virtual void v_PhysDeriv
            (const int dir,
             const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &out_d0);

            STD_REGIONS_EXPORT virtual void v_PhysDirectionalDeriv
            (const Array<OneD, const NekDouble>& inarray,
             const Array<OneD, const NekDouble>& direction,
             Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv
            (const Array<OneD,
             const NekDouble>& inarray,
             Array<OneD, NekDouble> &out_d1,
             Array<OneD, NekDouble> &out_d2,
             Array<OneD, NekDouble> &out_d3);
            
            STD_REGIONS_EXPORT virtual void   v_StdPhysDeriv
            (const int dir,
             const Array<OneD, const NekDouble>& inarray,
             Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate
            (const Array<OneD,
             const NekDouble>& coords,
             const Array<OneD, const NekDouble> & physvals);
            
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate
            (const Array<OneD, DNekMatSharedPtr >& I,
             const Array<OneD, const NekDouble> & physvals);
            
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluateBasis
            (const Array<OneD, const NekDouble>& coords, int mode);

            STD_REGIONS_EXPORT virtual void v_LocCoordToLocCollapsed(
                                        const Array<OneD, const NekDouble>& xi,
                                        Array<OneD, NekDouble>& eta);

            STD_REGIONS_EXPORT virtual void v_LocCollapsedToLocCoord(
                                        const Array<OneD, const NekDouble>& eta,
                                        Array<OneD, NekDouble>& xi);

            STD_REGIONS_EXPORT virtual void v_FillMode(const int mode,
                                                  Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                                                  const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                                                  const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_GetCoords(
                                                   Array<OneD, NekDouble> &coords_0,
                                                   Array<OneD, NekDouble> &coords_1,
                                                   Array<OneD, NekDouble> &coords_2);

            STD_REGIONS_EXPORT virtual void v_GetCoord(
                                          const Array<OneD, const NekDouble>& Lcoord,
                                          Array<OneD, NekDouble> &coord);
            
            STD_REGIONS_EXPORT virtual int v_GetCoordim(void);

            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                                          Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                                          Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual int v_GetVertexMap(int localVertexId,
                                            bool useCoeffPacking = false);

            STD_REGIONS_EXPORT virtual void v_GetTraceToElementMap(
                 const int                  tid,
                 Array<OneD, unsigned int> &maparray,
                 Array<OneD,          int> &signarray,
                 Orientation                traceOrient = eForwards,
                 int                        P = -1,
                 int                        Q = -1);

            STD_REGIONS_EXPORT  virtual void v_GetTraceInteriorToElementMap(
                const int                  eid,
                Array<OneD, unsigned int>& maparray,
                Array<OneD,          int>& signarray,
                const Orientation          traceOrient = eForwards);


            STD_REGIONS_EXPORT virtual void v_GetTraceNumModes(
                                              const int fid,
                                              int &numModes0,
                                              int &numModes1,
                                              Orientation traceOrient
                                                     = eDir1FwdDir1_Dir2FwdDir2);
            
            STD_REGIONS_EXPORT virtual void v_GetVertexPhysVals
            (const int vertex,
             const Array<OneD, const NekDouble> &inarray,
             NekDouble &outarray);

            STD_REGIONS_EXPORT virtual void v_MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac
                    (const Array<OneD, const NekDouble>& inarray,
                     Array<OneD, NekDouble> &outarray);


            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac
                    (const Array<OneD, const NekDouble>& inarray,
                     Array<OneD, NekDouble> &outarray,
                     bool multiplybyweights = true);

            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_SumFac(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void
                    v_IProductWRTDirectionalDerivBase_SumFac(
                const Array<OneD, const NekDouble>& direction,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_MassMatrixOp
                     (const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                      const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp
                     (const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                      const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_SVVLaplacianFilter
                      (Array<OneD,NekDouble> &array,
                       const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_ExponentialFilter(
                                          Array<OneD, NekDouble> &array,
                                    const NekDouble        alpha,
                                    const NekDouble        exponent,
                                    const NekDouble        cutoff);

            STD_REGIONS_EXPORT virtual void v_ReduceOrderCoeffs(
                                          int numMin,
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                                          const int k1, const int k2,
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_WeakDerivMatrixOp(
                                          const int i,
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_WeakDirectionalDerivMatrixOp(
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_MassLevelCurvatureMatrixOp(
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_LinearAdvectionDiffusionReactionMatrixOp(
                                          const Array<OneD,
                                          const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey,
                                          bool addDiffusionTerm=true);
            
            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp_MatFree(
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp_MatFree_Kernel(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                      Array<OneD,       NekDouble> &wsp);

            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp_MatFree(
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr
               v_BuildInverseTransformationMatrix(const DNekScalMatSharedPtr &
                                                  m_transformationmatrix);

            STD_REGIONS_EXPORT virtual void v_GetSimplexEquiSpacedConnectivity(
                Array<OneD, int> &conn,
                bool              standard = true);
        };

        typedef std::shared_ptr<StdExpansion> StdExpansionSharedPtr;
        typedef std::vector< StdExpansionSharedPtr > StdExpansionVector;

        /**
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
         *  @param inarray array of the function discretely evaluated at the
         *  quadrature points
         *
         *  @param outarray array of the function coefficieints
         */
        inline void StdExpansion::FwdTrans (const Array<OneD, const NekDouble>& inarray,
                        Array<OneD, NekDouble> &outarray)
        {
            v_FwdTrans(inarray,outarray);
        }

    } //end of namespace
} //end of namespace

#endif //STANDARDDEXPANSION_H
