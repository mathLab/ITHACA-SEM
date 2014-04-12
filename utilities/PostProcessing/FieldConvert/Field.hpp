////////////////////////////////////////////////////////////////////////////////
//
//  File: Field.hpp
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
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/DisContField3D.h>


using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        struct Field {
            Field() : m_verbose(false),m_declareExpansionAsContField(false),
                      m_declareExpansionAsDisContField(false),
                      m_writeBndFld(false){}
            
            bool m_verbose;
            vector<LibUtilities::FieldDefinitionsSharedPtr> m_fielddef;
            vector<vector<double> > m_data;
            vector<MultiRegions::ExpListSharedPtr> m_exp;

            bool m_declareExpansionAsContField;
            bool m_declareExpansionAsDisContField;
            
            LibUtilities::CommSharedPtr          m_comm;
            LibUtilities::SessionReaderSharedPtr m_session;
            SpatialDomains::MeshGraphSharedPtr   m_graph;
            LibUtilities::FieldIOSharedPtr       m_fld;
            map<string, vector<string> >         m_inputfiles;

            bool                 m_writeBndFld; 
            vector<unsigned int> m_bndRegionsToWrite;
            bool                 m_fldToBnd; 


            MultiRegions::ExpListSharedPtr AppendExpList(string var = "DefaultVar")
            {
                MultiRegions::ExpListSharedPtr tmp;
                switch (m_graph->GetMeshDimension())
                {
                case 1:
                    {
                        if (m_fielddef[0]->m_numHomogeneousDir == 1)
                        {
                            if(m_declareExpansionAsContField||
                               m_declareExpansionAsDisContField)
                            {
                                ASSERTL0(false,"ContField2DHomogeneous1D or DisContField2DHomogenenous1D has not been implemented");
                            }
                            
                            MultiRegions::ExpList2DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ExpList2DHomogeneous1D>(m_exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::
                                ExpList2DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                            
                        }
                        else if (m_fielddef[0]->m_numHomogeneousDir == 2)
                        {
                            if(m_declareExpansionAsContField)
                            {
                                MultiRegions::ContField3DHomogeneous2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ContField3DHomogeneous2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    ContField3DHomogeneous2D>::
                                    AllocateSharedPtr(*tmp2);
                            }
                            else  if(m_declareExpansionAsContField)
                            {
                                MultiRegions::DisContField3DHomogeneous2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    DisContField3DHomogeneous2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    DisContField3DHomogeneous2D>::
                                    AllocateSharedPtr(*tmp2);
                            }
                            else 
                            {
                                MultiRegions::ExpList3DHomogeneous2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList3DHomogeneous2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    ExpList3DHomogeneous2D>::
                                    AllocateSharedPtr(*tmp2);
                            }

                            
                        }
                        else
                        {
                            if(m_declareExpansionAsContField)
                            {
                                MultiRegions::ContField1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ContField1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::ContField1D>::
                                    AllocateSharedPtr(m_session,m_graph,var);
                            }
                            else if(m_declareExpansionAsDisContField)
                            {
                                MultiRegions::DisContField1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    DisContField1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::DisContField1D>::
                                    AllocateSharedPtr(m_session,m_graph,var);
                            }
                            else
                            {
                                MultiRegions::ExpList1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::ExpList1D>::
                                    AllocateSharedPtr(*tmp2);
                            }

                        }
                    }
                    break;
                case 2:
                    {   
                        if (m_fielddef[0]->m_numHomogeneousDir == 1)
                        {                            
                            if(m_declareExpansionAsContField)
                            {
                                MultiRegions::ContField3DHomogeneous1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ContField3DHomogeneous1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    ContField3DHomogeneous1D>::
                                    AllocateSharedPtr(*tmp2);

                                //WARNINGL0(false,"ContField is not copying BCs");
                            }
                            else  if(m_declareExpansionAsDisContField)
                            {
                                MultiRegions::DisContField3DHomogeneous1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    DisContField3DHomogeneous1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    DisContField3DHomogeneous1D>::
                                    AllocateSharedPtr(*tmp2);
                                //  WARNINGL0(false,"DisContField is not copying BCs");
                            }
                            else  
                            {
                                MultiRegions::ExpList3DHomogeneous1DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList3DHomogeneous1D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::
                                    ExpList3DHomogeneous1D>::
                                    AllocateSharedPtr(*tmp2);
                            }

                        }
                        else
                        {
                            if(m_declareExpansionAsContField)
                            {
                                MultiRegions::ContField2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ContField2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::ContField2D>::
                                    AllocateSharedPtr(*tmp2,m_graph,var);
                            }
                            else if(m_declareExpansionAsDisContField)
                            {
                                MultiRegions::DisContField2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    DisContField2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::DisContField2D>::
                                    AllocateSharedPtr(*tmp2,m_graph,var);
                            }
                            else
                            {
                                MultiRegions::ExpList2DSharedPtr tmp2 =
                                    boost::dynamic_pointer_cast<MultiRegions::
                                    ExpList2D>(m_exp[0]);
                                
                                tmp = MemoryManager<MultiRegions::ExpList2D>::
                                    AllocateSharedPtr(*tmp2);
                            }       
                        }
                    }
                    break;
                case 3:
                    {
                        if(m_declareExpansionAsContField)
                        {
                            MultiRegions::ContField3DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ContField3D>(m_exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::ContField3D>::
                                AllocateSharedPtr(*tmp2,m_graph,var);
                        }
                        else if(m_declareExpansionAsDisContField)
                        {
                            MultiRegions::DisContField3DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                DisContField3D>(m_exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::DisContField3D>::
                                AllocateSharedPtr(*tmp2,m_graph,var);
                        }
                        else
                        {
                            MultiRegions::ExpList3DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ExpList3D>(m_exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::ExpList3D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                    break;
                default:
                    ASSERTL0(false, "Expansion dimension not recognised");
                    break;
                }
                
                return tmp;
            }
            
        };

        typedef boost::shared_ptr<Field> FieldSharedPtr;
    }
}

