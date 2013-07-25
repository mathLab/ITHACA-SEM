///////////////////////////////////////////////////////////////////////////////
//
// File: Forcing.h
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGSPONGE
#define NEKTAR_SOLVERUTILS_FORCINGSPONGE

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

#include <MultiRegions/ExpList2D.h>     // for ExpList2D, etc
#include <MultiRegions/ExpList3D.h>     // for ExpList3D
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

namespace Nektar
{
    namespace SolverUtils
    {
        
        class ForcingSponge : public Forcing
        {
	public:
	    
	    friend class MemoryManager<ForcingSponge>;

            /// Creates an instance of this class
            SOLVER_UTILS_EXPORT static ForcingSharedPtr create() {
                return ForcingSharedPtr(new ForcingSponge());
                //ForcingSharedPtr p = MemoryManager<ForcingSponge>::AllocateSharedPtr();
		//return p;
            }
            SOLVER_UTILS_EXPORT ForcingSponge();
	    ///Name of the class
            static std::string className;
	    //SOLVER_UTILS_EXPORT ForcingSponge();
        protected:
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
		SpatialDomains::MeshGraphSharedPtr                pGraph);
            virtual void v_Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
	protected:
	    LibUtilities::SessionReaderSharedPtr                  m_Session;
	    Array<OneD, Array<OneD, NekDouble> >                  m_Sponge;
            Array<OneD, Array<OneD, NekDouble> >                  m_Refflow;
	    Array<OneD, Array<OneD, NekDouble> >                  m_Forcing;
            enum HomogeneousType
            {
                eHomogeneous1D,
                eHomogeneous2D,
                eHomogeneous3D,
                eNotHomogeneous
            }; 
            enum HomogeneousType                                  m_HomogeneousType;
	    int                                                   m_NumVariable;
            bool                                                  m_SingleMode;
            bool                                                  m_HalfMode;
	private:
	    void EvaluateFunction(
		Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
		LibUtilities::SessionReaderSharedPtr              pSession,
		std::string 					  pFieldName,
            	Array<OneD, NekDouble>&                           pArray,
            	const std::string&                                pFunctionName);
	    void v_ReadSpongeInfo(
		LibUtilities::SessionReaderSharedPtr              pSession,
		Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
                SpatialDomains::MeshGraphSharedPtr                pGraph);
        };
        
    }
}
// Hui XU  21 Jul 2013 Created 
#endif
