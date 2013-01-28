////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "MeshElements.h"
#include "ProcessJac.h"

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

#include <vector>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessJac::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "jac"), ProcessJac::create,
                "Process elements based on values of Jacobian.");

        ProcessJac::ProcessJac(MeshSharedPtr m) : ProcessModule(m)
        {
            
        }

        ProcessJac::~ProcessJac()
        {
            
        }
        
        void ProcessJac::Process()
        {
            if (m->verbose)
            {
                cout << "ProcessJac: Calculating Jacobians..." << endl;
            }

            vector<ElementSharedPtr> &el = m->element[m->expDim];
            
            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                // Create elemental geometry.
                SpatialDomains::GeometrySharedPtr geom = el[i]->GetGeom(m->spaceDim);
                
                // Define basis key using MeshGraph functions. Need a better
                // way of determining the number of modes!
                LibUtilities::BasisKeyVector b = 
                    SpatialDomains::MeshGraph::DefineBasisKeyFromExpansionType(
                        geom, SpatialDomains::eModified, 5);
                
                Array<OneD, LibUtilities::BasisSharedPtr> basis(m->expDim);

                // Generate basis functions.
                for (int j = 0; j < m->expDim; ++j)
                {
                    basis[j] = LibUtilities::BasisManager()[b[j]];
                }
                
                // Generate geometric factors.
                SpatialDomains::GeomFactorsSharedPtr gfac = 
                    geom->GetGeomFactors(basis);
                
                // Get the Jacobian and, if it is negative, print a warning
                // message.
                Array<OneD, NekDouble> jac = gfac->GetJac();
                NekDouble d = Vmath::Vmin(jac.num_elements(),&jac[0],1);
                if (d <= 0)
                {
                    cout << "Negative Jacobian in element " 
                         << el[i]->GetId() << " (value = "
                         << d << ")" << endl;
                }
            }
        }
    }
}
