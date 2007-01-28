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
namespace Nektar
{
    namespace StdRegions
    {
    
        class StdMatContainer;
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
         *  contains the defintiion of common data and common routine to all
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
			 const LibUtilities::BasisKey &Bc, 
			 int numcoeffs, double *coeffs, double *phys);

            /** \brief Copy Constructor */
            StdExpansion(const StdExpansion &T);

            /** \brief Destructor */
            virtual ~StdExpansion();

            // Standard Expansion Routines Applicable Regardless of Region

            /** \brief Return the number of 1D basis used in expansion */
            inline int GetNumBases() const
            {
                return m_numbases;
            }

            /** \brief Return a double pointer to the basis in the \a dir 
	     *  direction 
	     */
            inline const LibUtilities::BasisSharedPtr GetBasis(const int dir) const
            {
                ASSERTL1(dir <= m_numbases,"Base_id was larger than _numbases");
                return(m_base[dir]);
            }

            /** \brief Return the number of coefficient used in the total 
	     *  expansion 
	     */
            inline int GetNcoeffs(void)
            {
                return(m_ncoeffs);
            }

            /** \brief Return the double pointer to the coefficient array 
	     *  \a m_coeffs 
	     */
            inline BstShrDArray GetCoeffs(void)
            {
                return(m_coeffs);
            }

            /** \brief Set the coefficient array \a (this).m_coeffs to the values
	     *  given by \a coeffs
	     */
            inline void SetCoeffs(double *coeffs)
            {
                Vmath::Vcopy(m_ncoeffs, coeffs, 1, &m_coeffs[0], 1);
            }

            /** \brief Set the \a i th coefficient \a (this).m_coeffs[i] to the
	     *   values given by \a coeff
	     */
            inline void SetCoeff(const int i, double coeff)
            {
                m_coeffs[i] = coeff;
            }

            /** \brief Return the double pointer to the physical quarature points
	     *  array \a m_phys
	     */
            inline BstShrDArray  GetPhys(void)
            {
                return(m_phys);
            }

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

            /** \brief Sets the value of \a m_phys to the input argument 
	     *  \a phys 
	     */
            inline void  SetPhys(const double *phys)
            {
                int nqtot = GetTotPoints();

                Vmath::Vcopy(nqtot, phys, 1, &m_phys[0], 1);
            }


            /** Return the basis type in the \a dir direction using the enum
	     *  BasisType list
	     */
            inline LibUtilities::BasisType GetBasisType(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(m_base[dir]->GetBasisType());
            }

            /** Return the expansion order of the 1D expansion in the \a dir
	     *  direction
	     */
            inline int GetBasisNumModes(const int dir) const
            {
                ASSERTL1(dir < m_numbases,"dir is larger than m_numbases");
                return(m_base[dir]->GetNumModes());
            }

            /** Return the quadrature type in the \a dir direction using the
	     *  enum PointType list
	     */
            inline LibUtilities::PointsType GetPointsType(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(m_base[dir]->GetPointsType());
            }

            /** Return the number of quadrature points in the \a dir
	     *  direction
	     */
            inline int GetNumPoints(const int dir) const
            {
                ASSERTL1(dir < m_numbases, "dir is larger than m_numbases");
                return(ExpPointsProperties(dir)->GetNumPoints());
            }

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
	    int GetNverts()
	    {
		return v_GetNverts();
	    }

	    int GetNedges()
	    {
		return v_GetNedges();
	    }

	    int GetEdgeNcoeffs(const int i)
	    {
		return v_GetEdgeNcoeffs(i);
	    }

	    LibUtilities::BasisType GetEdgeBasisType(const int i)
	    {
		return v_GetEdgeBasisType(i);
	    }


	    int GetNfaces()
	    {
		return v_GetNfaces();
	    }
	    
            ShapeType DetShapeType()
            {
                return v_DetShapeType();
            }

            void  BwdTrans (double *outarray)
            {
                v_BwdTrans (outarray);
            }

            void  FwdTrans (const double *inarray)
            {
                v_FwdTrans(inarray);
            }

            double Integral(const double *inarray )
            {
                return v_Integral(inarray);
            }

            double Evaluate(const double * coords)
            {
                return v_Evaluate(coords);
            }

            void FillMode(const int mode, double * outarray)
            {
                v_FillMode(mode, outarray);
            }

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

            GeomType GeoFacType(void)
            {
                return v_GeoFacType();
            }

            // virtual functions related to LocalRegions
            boost::shared_ptr<LocalRegions::MetricRelatedInfo> GenGeoFac(void)
            {
		return v_GenGeoFac();
            }

            void SetGeoFac(boost::shared_ptr<LocalRegions::MetricRelatedInfo> minfo)
            {
                v_SetGeoFac(minfo);
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


            // Overloaded Operators
            friend bool operator == (const StdExpansion &x,const StdMatContainer *y);
            friend bool operator != (const StdExpansion &x,const StdMatContainer *y);
            friend bool operator == (const StdMatContainer *x,const StdExpansion &y);
            friend bool operator != (const StdMatContainer *x,const StdExpansion &y);


            void Interp1D(const LibUtilities::BasisKey *fbasis0, const double *from,
                const LibUtilities::BasisKey *tbasis0, double *to);

            void Interp2D(const LibUtilities::BasisKey *fbasis0, const LibUtilities::BasisKey *fbasis1,
                const double *from,   const LibUtilities::BasisKey *tbasis0,
                const LibUtilities::BasisKey* tbasis1, double *to);

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

            int   m_numbases;   /**< Number of 1D basis defined in expansion */
	    LibUtilities::BasisSharedPtr *m_base; /**< Bases needed for the expansion */	    

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
		return (LibUtilities::BasisType) NULL;
	    }


            virtual ShapeType v_DetShapeType()                = 0;


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
                NEKERROR(ErrorUtil::efatal, "This function is has not "
			 "been defined for this element");
            }

            virtual DNekMatSharedPtr v_GenLapMatrix()
            {
                NEKERROR(ErrorUtil::efatal, "This function is has not "
			 "been defined for this element");
            }

            virtual DNekMatSharedPtr v_GenNBasisTransMatrix()
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid "
			 "for nodal expansions");
            }

            virtual StdMatContainer *v_GetNBasisTransMatrix()
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid "
			 "for nodal expansions");
                return NULL;
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

            virtual GeomType v_GeoFacType()
            {
                return eRegular;
            }

            // virtual functions related to LocalRegions
            virtual boost::shared_ptr<LocalRegions::MetricRelatedInfo> v_GenGeoFac()
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return boost::shared_ptr<LocalRegions::MetricRelatedInfo>();
            }
	    
            virtual void v_SetGeoFac(boost::shared_ptr<LocalRegions::MetricRelatedInfo>
				     minfo)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            }

	    virtual void v_SetInvInfo(StdMatContainer *mat, MatrixType Mform)
	    {
		NEKERROR(ErrorUtil::efatal,"SetInvInfo not set for this expansion type");
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




