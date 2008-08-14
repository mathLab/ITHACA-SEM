///////////////////////////////////////////////////////////////////////////////
//
// File MultiRegsions.hpp
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
// Description: Multiregion overal header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_H
#define MULTIREGIONS_H

#include <LocalRegions/LocalRegions.hpp>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <SpatialDomains/BoundaryConditions.h>
#include <LibUtilities/BasicUtils/Metis.hpp>
//#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
//#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>


namespace Nektar
{ 
    /** \brief The namespace associated with the the MultiRegions library 
     * (\ref pageMultiRegions "MultiRegions introduction")
     */
    namespace MultiRegions
    {
        /** \page pageMultiRegions The MultiRegions library
         *
         * In the MultiRegions library, all classes and routines are related to 
         * the process of assembling a global spectral/hp expansion out of local
         * elemental contributions are bundled together. The most important entities 
         * of this library are the base class ExpList and its daughter classes. 
         * These classes all are the abstraction of a multi-elemental spectral/hp 
         * element expansion. Three different types of multi-elemental expansions 
         * can be distinguished:<BR><BR>
         * <b> A collection of local expansions.</b><BR>
         * This collection is just a list of local expansions, without any coupling
         * between the expansions on the different elements, and can be formulated 
         * as:
         * \f[ u^{\delta}(\boldsymbol{x})=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}) \f]
         * where 
         * - \f${N_{\mathrm{el}}}\f$ is the number of elements,
         * - \f$N^{e}_m\f$ is the number of local expansion modes within the 
         *   element \f$e\f$,
         * - \f$\phi_n^e(\boldsymbol{x})\f$ is the \f$n^{th}\f$ local expansion 
         *   mode within the element \f$e\f$,
         * - \f$\hat{u}_n^e\f$ is the \f$n^{th}\f$ local expansion coefficient 
         * - within the element \f$e\f$.
         *
         * These types of expansion are represented by the classes ExpList1D, 
         * ExpList2D and ExpList3D, depending on the dimension of the problem.
         *
         * <b> A multi-elemental continuous global expansion.</b><BR>
         * All local elemental expansions are now connected to form a global 
         * spectral/hp representation. In this case, by imposing \f$C^0\f$ 
         * continuity across the element interfaces, the expansion is chosen to be 
         * globally continuous. This type of global expansion can be defined as:
         * \f[u^{\delta}(\boldsymbol{x})=\sum_{n=0}^{N_{\mathrm{dof}}-1}\hat{u}_n
         * \Phi_n(\boldsymbol{x})=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x})\f]
         * where 
         * - \f$N_{\mathrm{dof}}\f$ refers to the number of global modes,
         * - \f$\Phi_n(\boldsymbol{x})\f$ is the \f$n^{th}\f$ global expansion mode,
         * - \f$\hat{u}_n\f$ is the \f$n^{th}\f$ global expansion coefficient.
         *
         * Typically, a mapping array to relate the global degrees of freedom 
         * \f$\hat{u}_n\f$ and local degrees of freedom \f$\hat{u}_n^e\f$ is 
         * required to assemble the global expansion out of the local
         * contributions. <BR>
         * These types of expansion are represented by the classes ContExpList1D, 
         * ContExpList2D and ContExpList3D, depending on the dimension of the 
         * problem.
         * - In order to solve (second-order) partial differential equations, 
         *   information about the boundary conditions should be incorporated in 
         *   the expansion. In case of a standard Galerkin implementation, the 
         *   Dirichlet boundary conditions can be enforced by lifting a known 
         *   solution satisfying these conditions, leaving a homogeneous Dirichlet 
         *   problem to be solved. If we denote the unknown solution by 
         *   \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet boundary
         *   conditions by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$ then we can 
         *   decompose the solution \f$u^{\delta}(\boldsymbol{x})\f$ into the form
         *   \f[u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
         *   u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
         *   \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
         *   \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}-1}
         *   \hat{u}_n^{\mathcal{H}}\Phi_n(\boldsymbol{x}_i).\f]
         *   Implementation-wise, the known solution can be lifted by ordering the 
         *   known degrees of freedom \f$\hat{u}_n^{\mathcal{H}}\f$ first in the 
         *   global solution array \f$\boldsymbol{\hat{u}}\f$.<BR>
         *   This type of global continuous expansion which incorporates the 
         *   boundary conditions are represented by the classes ContField1D, 
         *   ContField2D and ContField3D. Objects of these classes should be used 
         *   when solving partial differential equations using a standard Galerkin 
         *   approach.
         *
         * <b>A multi-elemental discontinuous global expansion.</b><BR>
         * The expansions are represented by the classes DiscContField1D, 
         * DiscContField2D and DiscContField3D. Objects of these classes should 
         * be used when solving partial differential equations using a 
         * discontinuous Galerkin approach.
         *
         * 
         * Furthermore, we have two more sets of classes:
         * - The class LocalToGlobalBndryMap and its daughter classes.<BR>
         *   These classes are an abstraction of the mapping from local to global 
         *   degrees of freedom and contain one or both of the following mapping 
         *   arrays:
         *   - map[\f$e\f$][\f$n\f$]<BR>
         *     This array contains the index of the global degree of freedom 
         *     corresponding to the \f$n^{th}\f$ local expansion mode within the 
         *     \f$e^{th}\f$ element.
         *   - bmap[\f$e\f$][\f$n\f$]<BR>
         *     This array contains the index of the global degree of freedom 
         *     corresponding to the \f$n^{th}\f$ local boundary mode within the 
         *     \f$e^{th}\f$ element.
         *   Next to the mapping array, these classes also contain routines to 
         *   assemble the global system from the local contributions, and other 
         *   routines to transform between local and global level.
         * - The classes GlobalLinSys and GlobalLinSysKey.<BR>
         *   The class GlobalLinSys is an abstraction of the global system matrix 
         *   resulting from the global assembly procedure. Depending of the choice 
         *   to statically condense the global matrix or not, the relevant blocks 
         *   are stored as a member of this class. Given a proper right hand side 
         *   vector, this class also contains a routine to solve the resulting matrix 
         *   system.<BR>
         *   The class GlobalLinSysKey represents a key which uniquely defines a 
         *   global matrix. This key can be used to construct or retrieve the global 
         *   matrix associated to a certain key. 
         *
         * More information about the implementation of connectivity between elements
         * in Nektar++ can be found \ref pageConnectivity "here".
         */
        
        // multiregion stuff here
        enum TransState
        {
            eNotSet,      ///< No transformed state set 
            eLocal,       ///< Local  Modal space array contains "true" expansion values
            eContinuous,  ///< Continuous Modal space array contains "true" expansion values
            eLocalCont,   ///< Both local and continuous space array contains "true" Expansion values 
        };
        
        enum GlobalSysSolnType
        {
            eNoSolnType,    ///< No Solution type specified
            eDirectFullMatrix,
            eDirectStaticCond,
        };

        // Orientation of adjacent edge for use with boundary
        // constraints
        enum AdjacentTraceOrientation
        {
            eAdjacentEdgeIsForwards,
            eAdjacentEdgeIsBackwards
        };
        
        
        const char* const GlobalSysSolnTypeMap[] = 
            {
            "No Solution Type",
            "Direct Solve: Full Matrix",
            "Direct Solve: Static Condensation"
        };

        
        typedef std::vector<SpatialDomains::BoundaryConditionType>  BndTypesVector;
        typedef std::vector<SpatialDomains::BoundaryConditionType>::iterator BndTypesVectorIter;

    }// end of namespace
}// end of namespace

#endif
