///////////////////////////////////////////////////////////////////////////////
//
// File LinearisedAdvection.h
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
// Description: TBA
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_LINEARISEDADVECTION_H
#define NEKTAR_SOLVERS_LINEARISEDADVECTION_H

#include <SpatialDomains/MeshComponents.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ExpListHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>
//#include <Auxiliary/EquationSystem.h>


#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>

//#define TIMING

//#ifdef TIMING
//#include <time.h>
//#include <sys/time.h>
//#endif


namespace Nektar
{     


    class LinearisedAdvection: public AdvectionTerm
    {
		enum FloquetMatType
        {
            eForwardsCoeff,
            eForwardsPhys
        };
		
		/// A map between  matrix keys and their associated block
        /// matrices.
        typedef map< FloquetMatType, DNekBlkMatSharedPtr> FloquetBlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef boost::shared_ptr<FloquetBlockMatrixMap> FloquetBlockMatrixMapShPtr;
		
    public:
        friend class MemoryManager<LinearisedAdvection>;

        /// Creates an instance of this class
        static AdvectionTermSharedPtr create(
                                const LibUtilities::SessionReaderSharedPtr& pSession,
                                const SpatialDomains::MeshGraphSharedPtr& pGraph) {
            AdvectionTermSharedPtr p = MemoryManager<LinearisedAdvection>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

	protected:
        //Storage of the base flow
        Array<OneD, MultiRegions::ExpListSharedPtr>     m_base;
		
		//Auxiliary base flow for half mode analysis
		Array<OneD, MultiRegions::ExpListSharedPtr>     m_base_aux;
		//number of slices
		int                                             m_slices;
		//period length
		NekDouble										m_period;
		//interpolation vector
		Array<OneD, Array<OneD, NekDouble> >			m_interp;
		//auxiliary variables
		LibUtilities::NektarFFTSharedPtr				m_FFT;
		Array<OneD,NekDouble>							m_tmpIN;
		Array<OneD,NekDouble>							m_tmpOUT;
		bool											    m_useFFTW;
		bool m_SingleMode;			 ///< flag to determine if use single mode or not
		bool m_HalfMode;		         ///< flag to determine if use half mode or not
		bool m_MultipleModes;		 ///< flag to determine if use multiple mode or not
		
		DNekBlkMatSharedPtr GetFloquetBlockMatrix(FloquetMatType mattype, bool UseContCoeffs = false) const;
		DNekBlkMatSharedPtr GenFloquetBlockMatrix(FloquetMatType mattype, bool UseContCoeffs = false) const;
		FloquetBlockMatrixMapShPtr       m_FloquetBlockMat;

		
        LinearisedAdvection(const LibUtilities::SessionReaderSharedPtr&        pSession,
							const SpatialDomains::MeshGraphSharedPtr&          pGraph);

        virtual ~LinearisedAdvection();

        virtual void v_InitObject();

        void SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh);
		void UpdateBase(const NekDouble m_slices,
						Array<OneD, const NekDouble> &inarray,
						Array<OneD, NekDouble> &outarray,
						const NekDouble m_time,
						const NekDouble m_period);
		void DFT(const string file, const NekDouble m_slices);
		
        /// Import Base flow
		void ImportFldBase(std::string pInfile,
						   SpatialDomains::MeshGraphSharedPtr pGraph, int cnt);
        void ImportFldBase(std::string pInfile,
                SpatialDomains::MeshGraphSharedPtr pGraph);
		
		/// Write field data to the given filename.
        void WriteFldBase(std::string &outname);
		
        /// Write input fields to the given filename.
        void WriteFldBase(
					  std::string &outname,
					  MultiRegions::ExpListSharedPtr &field,
					  Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs,
					  Array<OneD, std::string> &variables);
		
    private:
        //Function for the evaluation of the linearised advective terms
        virtual void v_ComputeAdvectionTerm(
                         Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                         const Array<OneD, Array<OneD, NekDouble> > &pV,
                         const Array<OneD, const NekDouble> &pU,
                         Array<OneD, NekDouble> &pOutarray,
                         int pVelocityComponent,
						 NekDouble m_time,
                         Array<OneD, NekDouble> &pWk);
		
		///Parameter for homogeneous expansions
        enum HomogeneousType
        {
            eHomogeneous1D,
            eHomogeneous2D,
            eHomogeneous3D,
            eNotHomogeneous
        };
		
        bool m_useFFT;               ///< flag to determine if use or not the FFT for transformations

        enum HomogeneousType m_HomogeneousType;
		
        NekDouble m_LhomX; ///< physical length in X direction (if homogeneous)
        NekDouble m_LhomY; ///< physical length in Y direction (if homogeneous)
        NekDouble m_LhomZ; ///< physical length in Z direction (if homogeneous)
		
        int m_npointsX;    ///< number of points in X direction (if homogeneous)
        int m_npointsY;    ///< number of points in Y direction (if homogeneous)
        int m_npointsZ;    ///< number of points in Z direction (if homogeneous)
		
        int m_HomoDirec;   ///< number of homogenous directions
		
		int m_NumMode;     ///< Mode to use in case of single mode analysis

		SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
    };
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
