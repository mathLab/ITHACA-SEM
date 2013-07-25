////////////////////////////////////////////////////////////////////////////////
//
//  File: Field.h
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
//  Description: Field converter module base classes.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ExpList.h>

#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>


using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        struct Field {
            Field() : verbose(false) {}

            vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
            vector<vector<double> > data;
            vector<MultiRegions::ExpListSharedPtr> exp;
            
            LibUtilities::SessionReaderSharedPtr session;
            SpatialDomains::MeshGraphSharedPtr graph;
            
           

            MultiRegions::ExpListSharedPtr AppendExpList()
            {
                MultiRegions::ExpListSharedPtr tmp;
                switch (graph->GetMeshDimension())
                {
                    case 1:
                    {
                        if (fielddef[0]->m_numHomogeneousDir == 1)
                        {
                            MultiRegions::ExpList2DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList2DHomogeneous1D>(exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::
                                ExpList2DHomogeneous1D>::
                                    AllocateSharedPtr(*tmp2);

                        }
                        else if (fielddef[0]->m_numHomogeneousDir == 2)
                        {
                            MultiRegions::ExpList3DHomogeneous2DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList3DHomogeneous2D>(exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous2D>::
                                    AllocateSharedPtr(*tmp2);
                        }
                        else
                        {
                            MultiRegions::ExpList1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList1D>(exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::ExpList1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                        break;
                    case 2:
                    {   
                        if (fielddef[0]->m_numHomogeneousDir == 1)
                        {
                            
                            MultiRegions::ExpList3DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList3DHomogeneous1D>(exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous1D>::
                                    AllocateSharedPtr(*tmp2);
                        }
                        else
                        {
                            MultiRegions::ExpList2DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList2D>(exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::ExpList2D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                        break;
                    case 3:
                    {
                        MultiRegions::ExpList3DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                                ExpList3D>(exp[0]);
                        
                        tmp = MemoryManager<MultiRegions::ExpList3D>::
                            AllocateSharedPtr(*tmp2);

                    }
                        break;
                    default:
                        ASSERTL0(false, "Expansion dimension not recognised");
                        break;
                }
                
                //data.push_back(newdata);
                //exp.push_back(tmp);
                //fielddef.push_back();
                
                return tmp;
            }
            
            bool verbose;
        };

        typedef boost::shared_ptr<Field> FieldSharedPtr;
    }
}

