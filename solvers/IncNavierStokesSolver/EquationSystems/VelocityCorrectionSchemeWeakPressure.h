///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionSchemeWeakPressure.h
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
// Description: Velocity Correction Scheme header with Weak Pressure Formulation
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEMEWEAKPRESSURE_H
#define NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEMEWEAKPRESSURE_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>

namespace Nektar
{
    class VCSWeakPressure: public VelocityCorrectionScheme
    {
    public:

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<VCSWeakPressure>::AllocateSharedPtr(
                    pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;


        /// Constructor.
        VCSWeakPressure(const LibUtilities::SessionReaderSharedPtr& pSession,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph);

        virtual ~VCSWeakPressure();

    protected:
        // Virtual functions
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        virtual void v_SetUpPressureForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        
        virtual void v_SolvePressure( const Array<OneD, NekDouble>  &Forcing);
        
        virtual std::string v_GetExtrapolateStr(void)
        {
            return "WeakPressure";
        }
        
        virtual std::string v_GetSubSteppingExtrapolateStr(const std::string &instr)
        {
            if(boost::iequals(instr,"SubStepping"))
            {
                return  "SubSteppingWeakPressure";
            }
            else
            {
                return instr;
            }
        }

    private:
        
    };

    typedef std::shared_ptr<VCSWeakPressure>
                VCSWeakPressureSharedPtr;

} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H
