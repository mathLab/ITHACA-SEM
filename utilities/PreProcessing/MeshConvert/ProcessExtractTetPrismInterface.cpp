////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessExtractTetPrismInterface.cpp
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
//  Description: Extract 
//
////////////////////////////////////////////////////////////////////////////////

#include "MeshElements.h"
#include "ProcessExtractTetPrismInterface.h"

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

#include <vector>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessExtractTetPrismInterface::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "extracttetprism"), ProcessExtractTetPrismInterface::create,
                "Process elements to extract the faces between tets and prisms.");

        ProcessExtractTetPrismInterface::ProcessExtractTetPrismInterface(MeshSharedPtr m) : ProcessModule(m)
        {
        }

        ProcessExtractTetPrismInterface::~ProcessExtractTetPrismInterface()
        {
        }
        
        void ProcessExtractTetPrismInterface::Process()
        {
            if (m_mesh->m_verbose)
            {
                cout << "ProcessExtractTetPrismInterface: Extracting interface... " << endl;
            }

            vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

            vector<int> pfaces;
            vector<int> tfaces;
            vector<int> inter;

            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                // Create elemental geometry.
                SpatialDomains::GeometrySharedPtr geom =
                    el[i]->GetGeom(m_mesh->m_spaceDim);

                int nFaces= el[i]->GetFaceCount();

                //Get the number of nodes on a face of a prism 
                if (nFaces ==5)
                {
                    //First triangular face
                    pfaces.push_back(el[i]->GetFace(1)->m_id);

                    //Second triangular face
                    pfaces.push_back(el[i]->GetFace(3)->m_id);

                }

                if (nFaces ==4)
                {
                    for(int j=0; j<nFaces; ++j)
                    {
                        tfaces.push_back(el[i]->GetFace(j)->m_id);
                    }
                }

            }
            
            //order the vectors
            sort(tfaces.begin(), tfaces.end());
            sort(pfaces.begin(), pfaces.end());

            //Find the faces that are in both vectors
            set_intersection(tfaces.begin(), tfaces.end(),
                             pfaces.begin(), pfaces.end(),
                             back_inserter(inter));
            std::vector<int>::iterator it;
            cout << "The intersection has " << (inter.size()) << " elements:\n";
            for (it=inter.begin(); it!=inter.end(); ++it)
                std::cout << ',' << *it;
            std::cout << '\n';
        }

    }
}
