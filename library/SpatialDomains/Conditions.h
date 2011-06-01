////////////////////////////////////////////////////////////////////////////////
//
//  File:  Conditions.h
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
#include <SpatialDomains/SpatialDomainsDeclspec.h>

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
            eRobin,
            ePeriodic
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(BoundaryConditionType type):
                m_boundaryConditionType(type)
            {
            }

            BoundaryConditionBase(BoundaryConditionType type, Equation userDefined):
                m_boundaryConditionType(type),  m_userDefined(userDefined)
            {
            }

            virtual ~BoundaryConditionBase()
            {};

            BoundaryConditionType GetBoundaryConditionType() const
            {
                return m_boundaryConditionType;
            }

            void SetUserDefined(Equation equation)
            {
                m_userDefined = equation;
            }

            Equation GetUserDefined() const
            {
                return m_userDefined;
            }


        protected:
            BoundaryConditionType m_boundaryConditionType;
            Equation m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {

             DirichletBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
            BoundaryConditionBase(eDirichlet, userDefined),
                m_dirichletCondition(eqn)
                {
                }

            Equation m_dirichletCondition;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_neumannCondition(eqn)
            {
            }

            Equation m_neumannCondition;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition( const std::string &a, const std::string &b, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eRobin, userDefined),
                m_robinFunction(a), m_robinPrimitiveCoeff(b)
            {
            }
                // \frac{\partial {u}}{\partial{n}} +
                // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
                Equation m_robinFunction;
                Equation m_robinPrimitiveCoeff;
        };


        struct PeriodicBoundaryCondition : public BoundaryConditionBase
        {
            PeriodicBoundaryCondition(const unsigned int n):
                BoundaryConditionBase(ePeriodic),
                m_connectedBoundaryRegion(n)
            {
            }

            unsigned int m_connectedBoundaryRegion;
        };

        typedef std::map<std::string, NekDouble> ParamMap;
        typedef std::map<std::string, std::string> FunctionMap;
        typedef std::vector<std::string> Variable;

        typedef std::map<int, Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;

        typedef Equation ForcingFunction;
        typedef boost::shared_ptr<ForcingFunction> ForcingFunctionShPtr;
        typedef boost::shared_ptr<const ForcingFunction> ConstForcingFunctionShPtr;
        typedef std::map<std::string, ForcingFunctionShPtr> ForcingFunctionsMap;

        typedef Equation ExactSolution;
        typedef boost::shared_ptr<ExactSolution> ExactSolutionShPtr;
        typedef boost::shared_ptr<const ExactSolution> ConstExactSolutionShPtr;
        typedef std::map<std::string, ExactSolutionShPtr> ExactSolutionMap;

        typedef Equation UserDefinedEqn;
        typedef boost::shared_ptr<UserDefinedEqn> UserDefinedEqnShPtr;
        typedef boost::shared_ptr<const UserDefinedEqn> ConstUserDefinedEqnShPtr;
        typedef std::map<std::string, UserDefinedEqnShPtr> UserDefinedEqnMap;

        typedef Equation InitialCondition;
        typedef boost::shared_ptr<InitialCondition> InitialConditionShPtr;
        typedef boost::shared_ptr<const InitialCondition> ConstInitialConditionShPtr;
        typedef std::map<std::string, InitialConditionShPtr> InitialConditionsMap;

        typedef std::map<std::string, std::string> SolverInfoMap;

        class BoundaryConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT BoundaryConditions(const MeshGraph *meshGraph);
            SPATIAL_DOMAINS_EXPORT ~BoundaryConditions();

            SPATIAL_DOMAINS_EXPORT void Read(const std::string& infilename);
            SPATIAL_DOMAINS_EXPORT void Read(TiXmlDocument &doc);

            SPATIAL_DOMAINS_EXPORT bool   CheckForParameter(const std::string &paramName);
            SPATIAL_DOMAINS_EXPORT static NekDouble GetParameter(const std::string &parmName);

            BoundaryRegionCollection &GetBoundaryRegions(void)
            {
                return m_boundaryRegions;
            }

            BoundaryConditionCollection &GetBoundaryConditions(void)
            {
                return m_boundaryConditions;
            }

            /// Get forcing function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            SPATIAL_DOMAINS_EXPORT ConstForcingFunctionShPtr GetForcingFunction(int indx) const;

            /// Get forcing function based on name of variable.
            SPATIAL_DOMAINS_EXPORT ConstForcingFunctionShPtr GetForcingFunction(const std::string &var) const;

            SPATIAL_DOMAINS_EXPORT bool ExactSolutionExists(int indx) const;
            SPATIAL_DOMAINS_EXPORT ConstExactSolutionShPtr GetExactSolution(int indx) const;
            SPATIAL_DOMAINS_EXPORT ConstExactSolutionShPtr GetExactSolution(const std::string &var) const;

            SPATIAL_DOMAINS_EXPORT bool UserDefinedEqnExists(const std::string &var) const;
            SPATIAL_DOMAINS_EXPORT ConstUserDefinedEqnShPtr GetUserDefinedEqn(int indx) const;
            SPATIAL_DOMAINS_EXPORT ConstUserDefinedEqnShPtr GetUserDefinedEqn(const std::string &var) const;

            /// Get initial condition function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            SPATIAL_DOMAINS_EXPORT bool InitialConditionExists(int indx) const;
            SPATIAL_DOMAINS_EXPORT ConstInitialConditionShPtr GetInitialCondition(int indx) const;

            /// Get initial condition function based on name of variable.
            SPATIAL_DOMAINS_EXPORT ConstInitialConditionShPtr GetInitialCondition(const std::string &var) const;
            /// Check to see if initial condition exists in list.
            SPATIAL_DOMAINS_EXPORT bool FoundInitialCondition(const std::string &var);

            const std::string &GetVariable(unsigned int indx)
            {
                ASSERTL0(0 <= indx && indx < m_variables.size(),"Variable index is out of range");
                return m_variables[indx];
            }

            inline int GetNumVariables() const
            {
                return m_variables.size();
            }

            static const ParamMap &GetParameters(void)
            {
                return m_parameters;
            }

            SPATIAL_DOMAINS_EXPORT const std::string &GetSolverInfo(const std::string &lhs);
            SPATIAL_DOMAINS_EXPORT bool SolverInfoExists(const std::string &property);

            SPATIAL_DOMAINS_EXPORT const std::string &GetFunction(const std::string &lhs);
            SPATIAL_DOMAINS_EXPORT Equation GetFunctionAsEquation(const std::string &lhs);

            /// Will look for the lhs equal to str and if found
            /// will return the function in str and return true.
            /// If not found it will return false and leave str
            /// as it was coming in.
            SPATIAL_DOMAINS_EXPORT bool SubstituteFunction(std::string &str);

        protected:
            void ReadSolverInfo(TiXmlElement *functions);
            void ReadParameters(TiXmlElement *parameters);
            void ReadVariables(TiXmlElement *variables);
            void ReadFunctions(TiXmlElement *conditions);
            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
            void ReadForcingFunctions(TiXmlElement *functions);
            void ReadInitialConditions(TiXmlElement *conditions);
            void ReadExactSolution(TiXmlElement *solution);
            void ReadUserDefinedEqn(TiXmlElement *functions);

            // Containers to hold conditions and associated data
            static ParamMap m_parameters;
            FunctionMap     m_functions;
            Variable        m_variables;
            BoundaryRegionCollection    m_boundaryRegions;
            BoundaryConditionCollection m_boundaryConditions;
            ForcingFunctionsMap         m_forcingFunctions;
            InitialConditionsMap        m_initialConditions;
            ExactSolutionMap            m_exactSolution;
            UserDefinedEqnMap           m_userDefinedEqn;

            SolverInfoMap               m_solverInfo; //< Solver Information

            /// The mesh graph to use for referencing geometry info.
            const MeshGraph *m_meshGraph;

        private:
            BoundaryConditions();
        };

        typedef boost::shared_ptr<BoundaryConditions> BoundaryConditionsSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

