///////////////////////////////////////////////////////////////////////////////
//
// File: RiemannSolver.h
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
// Description: Abstract base class for Riemann solvers with factory.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_RIEMANNSOLVER
#define NEKTAR_SOLVERUTILS_RIEMANNSOLVER

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <boost/function.hpp>
#include <boost/call_traits.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        typedef boost::function<Array<OneD, const NekDouble>& ()> RSScalarFuncType;
        typedef boost::function<Array<OneD, Array<OneD, const NekDouble> >& ()> RSVectorFuncType;
        typedef boost::function<NekDouble ()> RSParamFuncType;
        
        /*
        class RiemannSolverKey
        {
        public:
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddScalar(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                scalars[name] = boost::bind(func, obj);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddVector(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                vectors[name] = boost::bind(func, obj);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddParam(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                params[name] = boost::bind(func, obj);
            }
            
            /// Map of scalar function types.
            std::map<std::string, RSScalarFuncType> scalars;
            /// Map of vector function types.
            std::map<std::string, RSVectorFuncType> vectors;
            /// Map of parameter function types.
            std::map<std::string, RSParamFuncType > params;
        };
        */
        
        class RiemannSolver
        {
        public:
            void Solve(
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                      Array<OneD,       Array<OneD, NekDouble> > &flux);
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddScalar(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                m_scalars[name] = boost::bind(func, obj);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddVector(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                m_vectors[name] = boost::bind(func, obj);
            }
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            inline void AddParam(std::string name, FuncPointerT func, ObjectPointerT obj)
            {
                m_params[name] = boost::bind(func, obj);
            }


        protected:
            RiemannSolver();
            
            /// Map of scalar function types.
            std::map<std::string, RSScalarFuncType> m_scalars;
            /// Map of vector function types.
            std::map<std::string, RSVectorFuncType> m_vectors;
            /// Map of parameter function types.
            std::map<std::string, RSParamFuncType > m_params;
            
            virtual void v_Solve(
                const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                      Array<OneD,       Array<OneD, NekDouble> > &flux) = 0;

            bool CheckScalars(std::string name);
            bool CheckVectors(std::string name);
            bool CheckParams (std::string name);
        }; 
        
        /// A shared pointer to an EquationSystem object
        typedef boost::shared_ptr<RiemannSolver> RiemannSolverSharedPtr;
        /// Datatype of the NekFactory used to instantiate classes derived
        /// from the RiemannSolver class.
        typedef LibUtilities::NekFactory<std::string, RiemannSolver>
            RiemannSolverFactory;
        RiemannSolverFactory& GetRiemannSolverFactory();
    }
}

#endif
