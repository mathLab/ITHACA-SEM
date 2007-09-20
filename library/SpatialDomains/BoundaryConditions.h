////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/BoundaryConditions.h,v $
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include <SpatialDomains/Equation.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>

class TiXmlElement;
class TiXmlDocument;

class MeshGraph;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(BoundaryConditionType type):
                m_BoundaryConditionType(type)
            {
            };

            virtual ~BoundaryConditionBase(){};

            BoundaryConditionType GetBoundaryConditionType()
            {
                return m_BoundaryConditionType;
            };

        protected:
            BoundaryConditionType m_BoundaryConditionType;
        };

        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {
            DirichletBoundaryCondition(const std::string &eqn):
                m_DirichletCondition(eqn),
                BoundaryConditionBase(eDirichlet)
            {
            };

            Equation m_DirichletCondition;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(const std::string &eqn):
                m_NeumannCondition(eqn),
                BoundaryConditionBase(eNeumann)
            {
            };

            Equation m_NeumannCondition;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition(const std::string &a, const std::string &b):
                m_a(a), m_b(b),
                BoundaryConditionBase(eRobin)
            {
            }

            // u = a(x,y,z) + b(x,y,z)*\frac{\partial{u}}{\partial{n}}
            Equation m_a;
            Equation m_b;
        };

        typedef std::map<std::string, NekDouble> ParamMap;
        typedef std::map<std::string, std::string> FunctionMap;
        typedef std::vector<std::string> Variable;

        typedef std::vector<Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition> NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition> RobinBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr> BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap> BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        typedef Equation ForcingFunction;
        typedef boost::shared_ptr<ForcingFunction> ForcingFunctionShPtr;
        typedef boost::shared_ptr<const ForcingFunction> ConstForcingFunctionShPtr;
        typedef std::map<std::string, ForcingFunctionShPtr> ForcingFunctionsMap;

        typedef Equation ExactSolution;
        typedef boost::shared_ptr<ExactSolution> ExactSolutionShPtr;
        typedef boost::shared_ptr<const ExactSolution> ConstExactSolutionShPtr;
        typedef std::map<std::string, ExactSolutionShPtr> ExactSolutionMap;

        typedef Equation InitialCondition;
        typedef boost::shared_ptr<InitialCondition> InitialConditionShPtr;
        typedef boost::shared_ptr<const InitialCondition> ConstInitialConditionShPtr;
        typedef std::map<std::string, InitialConditionShPtr> InitialConditionsMap;

        class BoundaryConditions
        {
        public:
            BoundaryConditions(const MeshGraph *meshGraph);
            ~BoundaryConditions();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);

            static NekDouble GetParameter(const std::string &parmName);

            BoundaryRegionCollection &GetBoundaryRegions(void) 
            {
                return m_BoundaryRegions;
            }

            BoundaryConditionCollection &GetBoundaryConditions(void)
            {
                return m_BoundaryConditions;
            }

            /// Get forcing function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            ConstForcingFunctionShPtr GetForcingFunction(int indx) const;

            /// Get forcing function based on name of variable.
            ConstForcingFunctionShPtr GetForcingFunction(const std::string &var) const;

            ConstExactSolutionShPtr GetExactSolution(int indx) const;
            ConstExactSolutionShPtr GetExactSolution(const std::string &var) const;

            /// Get initial condition function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            ConstInitialConditionShPtr GetInitialCondition(int indx) const;

            /// Get initial condition function based on name of variable.
            ConstInitialConditionShPtr GetInitialCondition(const std::string &var) const;

            const std::string &GetVariable(unsigned int indx)
            {
                ASSERTL0(0 <= indx && indx < m_Variables.size(),"Variable index is out of range");
                return m_Variables[indx];
            }

            static const ParamMap &GetParameters(void)
            {
                return m_Parameters;
            }

            const std::string &GetFunction(const std::string &lhs);
            Equation GetFunctionAsEquation(const std::string &lhs);

            /// Will look for the lhs equal to str and if found
            /// will return the function in str and return true.
            /// If not found it will return false and leave str
            /// as it was coming in.
            bool SubstituteFunction(std::string &str);

        protected:
            void ReadParameters(TiXmlElement *parameters);
            void ReadVariables(TiXmlElement *variables);
            void ReadFunctions(TiXmlElement *conditions);
            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
            void ReadForcingFunctions(TiXmlElement *functions);
            void ReadInitialConditions(TiXmlElement *conditions);
            void ReadExactSolution(TiXmlElement *solution);

            // Containers to hold conditions and associated data
            static ParamMap m_Parameters;
            FunctionMap m_Functions;
            Variable m_Variables;
            BoundaryRegionCollection m_BoundaryRegions;
            BoundaryConditionCollection m_BoundaryConditions;
            ForcingFunctionsMap m_ForcingFunctions;
            InitialConditionsMap m_InitialConditions;
            ExactSolutionMap m_ExactSolution;

            /// The mesh graph to use for referencing geometry info.
            const MeshGraph *m_MeshGraph;

        private:
            BoundaryConditions();
        };
    };
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
    
