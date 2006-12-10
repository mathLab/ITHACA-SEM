///////////////////////////////////////////////////////////////////////////////
//
// File NodalBasisManager.cpp
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
// Description: Manager for Nodal basis
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/NodalBasisManager.h>

#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        char * NodalBasisFile [SIZE_NodalBasisType] = 
        {
            "NodalTriElec.dat",
                "NodalTriFekete.dat",
                "NodalTetElec.dat"
        };


        NodalBasisManager::NodalBasisManager()
        {
            int i,j,dim,nrecords;
            FILE * fp;
            char buf[NekConstants::kBufSize];
            BaryNodeContainer * rc;

            m_nodecontainer = new std::vector<BaryNodeContainer*>[SIZE_NodalBasisType];

            for(i=0;i<SIZE_NodalBasisType;i++)
            {
                sprintf(buf,"%s\0",NodalBasisFile[i]);
                if((fp = fopen(buf,"r")) == (FILE *) NULL)
                {
                    ASSERTL0(false, "Nodal Data File Opening Error");
                }

                fgets(buf,NekConstants::kBufSize,fp); // Dummy line
                fgets(buf,NekConstants::kBufSize,fp); // Dummy line
                fgets(buf,NekConstants::kBufSize,fp); // Dummy line
                sscanf(buf,"%d %d",&dim, &nrecords);

                switch(dim)
                {
                case 2:
                    for(j=0;j<nrecords;j++)
                    {
                        rc = new BaryNodeContainer();
                        rc->m_dim = dim;
                        NodalFileRead2D(fp,rc);
                        rc->set_ntype((NodalBasisType) i);
                        m_nodecontainer[i].push_back(rc);
                    }
                    break;
                case 3:
                    for(j=0;j<nrecords;j++)
                    {
                        rc = new BaryNodeContainer();
                        rc->m_dim = dim;
                        NodalFileRead3D(fp,rc);
                        rc->set_ntype((NodalBasisType) i);
                        m_nodecontainer[i].push_back(rc);
                    }
                    break;

                default:
                    ASSERTL0(false, "Dimension Check Violated");
                }
                fclose(fp);
            }
        }//end constructor

        NodalBasisManager::~NodalBasisManager()
        {
            std::vector<BaryNodeContainer*>::iterator def;

            for(int i=0;i<SIZE_NodalBasisType;i++)
            {
                for(def = m_nodecontainer[i].begin(); def != m_nodecontainer[i].end(); 
                    ++def)
                {
                    delete def[0];
                }
            }

            delete[] m_nodecontainer;

        }//end destructor

	
        int NodalBasisManager::GetNodePoints(NodalBasisType ntype, 
					     const int order, 
					     const double* &x, 
					     const double* &y, 
					     const double* &z)
        {
            std::vector<BaryNodeContainer*>::iterator def;
            BaryNodeContainer key;
            int  id_ntype = (int)ntype;

            key.m_ntype  = ntype;
            key.m_porder = order-1; // constructor order is 1D expansion order

            def = find(m_nodecontainer[id_ntype].begin(),
		       m_nodecontainer[id_ntype].end(),key);

            if(def != m_nodecontainer[id_ntype].end())
            {
                if(def[0]->m_xi[0]==NULL)
                {
                    switch(def[0]->m_dim)
                    {
                    case 2:
                        NodalPointExpand2D(def[0]);
                        break;
                    case 3:
                        NodalPointExpand3D(def[0]);
                        break;
                    default:
                        ASSERTL0(false, "Dimension Check Violated");
                    }
                }
                x = def[0]->m_xi[0];
                y = def[0]->m_xi[1];
                z = def[0]->m_xi[2];
            }
            else
	    {
                ASSERTL0(false, "Nodal order requested which is not available");
            }

            return def[0]->m_npts;
        }

        void NodalBasisManager::NodalFileRead2D(FILE *fp, BaryNodeContainer * rc)
        {
            int j;
            int sym1,sym3,sym6;
            double l1,l2,l3;
            char buf[NekConstants::kBufSize];

            fgets(buf,NekConstants::kBufSize,fp);
            sscanf(buf,"%d %d",&(rc->m_porder),&(rc->m_bnpts));

            rc->m_symm[0] = new int[rc->m_bnpts];
            rc->m_symm[1] = new int[rc->m_bnpts];
            rc->m_symm[2] = new int[rc->m_bnpts];

            rc->m_bpoints[0] = new double[rc->m_bnpts];
            rc->m_bpoints[1] = new double[rc->m_bnpts];
            rc->m_bpoints[2] = new double[rc->m_bnpts];

            for(j=0;j<rc->m_bnpts;j++)
            {
                fgets(buf,NekConstants::kBufSize,fp);
                sscanf(buf,"%d %d %d %lf %lf %lf",&sym1,&sym3,&sym6,&l1,&l2,&l3);
                rc->m_symm[0][j] = sym1;
                rc->m_symm[1][j] = sym3;
                rc->m_symm[2][j] = sym6;

                rc->m_bpoints[0][j] = l1;
                rc->m_bpoints[1][j] = l2;
                rc->m_bpoints[2][j] = l3;
            }
        }

        void NodalBasisManager::NodalFileRead3D(FILE *fp, BaryNodeContainer * rc)
        {
            int j;
            int sym1,sym3,sym6,sym12,sym24;
            double l1,l2,l3,l4;
            char buf[NekConstants::kBufSize];

            fgets(buf,NekConstants::kBufSize,fp);
            sscanf(buf,"%d %d",&(rc->m_porder),&(rc->m_bnpts));

            rc->m_symm[0] = new int[rc->m_bnpts];
            rc->m_symm[1] = new int[rc->m_bnpts];
            rc->m_symm[2] = new int[rc->m_bnpts];
            rc->m_symm[3] = new int[rc->m_bnpts];
            rc->m_symm[4] = new int[rc->m_bnpts];

            rc->m_bpoints[0] = new double[rc->m_bnpts];
            rc->m_bpoints[1] = new double[rc->m_bnpts];
            rc->m_bpoints[2] = new double[rc->m_bnpts];
            rc->m_bpoints[3] = new double[rc->m_bnpts];

            for(j=0;j<rc->m_bnpts;j++)
            {
                fgets(buf,NekConstants::kBufSize,fp);
                sscanf(buf,"%d %d %d %d %d %lf %lf %lf %lf",&sym1,&sym3,&sym6,&sym12,
                    &sym24,&l1,&l2,&l3,&l4);
                rc->m_symm[0][j] = sym1;
                rc->m_symm[1][j] = sym3;
                rc->m_symm[2][j] = sym6;
                rc->m_symm[3][j] = sym12;
                rc->m_symm[4][j] = sym24;

                rc->m_bpoints[0][j] = l1;
                rc->m_bpoints[1][j] = l2;
                rc->m_bpoints[2][j] = l3;
                rc->m_bpoints[3][j] = l4;
            }
        }

        void NodalBasisManager::NodalPointExpand2D(BaryNodeContainer * rc)
        {
            int i,j,sum=0;
            double a,b,c;

            rc->m_npts = (rc->m_porder+1)*(rc->m_porder+2)/2;
            rc->m_xi[0] = new double[rc->m_npts];
            rc->m_xi[1] = new double[rc->m_npts];

            for(i=0;i<rc->m_bnpts;i++)
            {
                if(rc->m_symm[0][i])
                {  
                    a = rc->m_bpoints[0][i]; 
                    b = rc->m_bpoints[1][i]; 
                    c = rc->m_bpoints[2][i];
                    rc->m_xi[0][sum] = 2.0*b - 1.0;
                    rc->m_xi[1][sum] = 2.0*c - 1.0;
                    sum++;
                }//end symmetry1

                if(rc->m_symm[1][i] == 1)
                {
                    for(j=0;j<3;j++)
                    {
                        a = rc->m_bpoints[perm3A_2d[j][0]][i];
                        b = rc->m_bpoints[perm3A_2d[j][1]][i];
                        c = rc->m_bpoints[perm3A_2d[j][2]][i];
                        rc->m_xi[0][sum] = 2.0*b - 1.0;
                        rc->m_xi[1][sum] = 2.0*c - 1.0;
                        sum++;
                    }
                }//end symmetry3a

                if(rc->m_symm[1][i] == 2)
                {
                    for(j=0;j<3;j++)
                    {
                        a = rc->m_bpoints[perm3B_2d[j][0]][i];
                        b = rc->m_bpoints[perm3B_2d[j][1]][i];
                        c = rc->m_bpoints[perm3B_2d[j][2]][i];
                        rc->m_xi[0][sum] = 2.0*b - 1.0;
                        rc->m_xi[1][sum] = 2.0*c - 1.0;
                        sum++;	  
                    }
                }//end symmetry3b


                if(rc->m_symm[2][i])
                {
                    for(j=0;j<6;j++)
                    {
                        a = rc->m_bpoints[perm6_2d[j][0]][i];
                        b = rc->m_bpoints[perm6_2d[j][1]][i];
                        c = rc->m_bpoints[perm6_2d[j][2]][i];
                        rc->m_xi[0][sum] = 2.0*b - 1.0;
                        rc->m_xi[1][sum] = 2.0*c - 1.0;
                        sum++;	  
                    }
                }//end symmetry6
            }//end npts

            NodalPointReorder2d(rc);

            ASSERTL1((sum==rc->m_npts),"sum not equal to npts");

        }//end NodalPointExpand2D


        void NodalBasisManager::NodalPointReorder2d(BaryNodeContainer * rc)
        {
            int i,j,istart,iend,sum=0;
            const int numvert = 3;
            const int numepoints = rc->m_porder-1;
            double tmpdouble;

            if(numepoints==0)
            {
                return;
            }

            // bubble sort for first edge
            istart = numvert + sum;
            iend = istart + numepoints;
            for(i=istart;i<iend;++i)
            {
                for(j=istart;j<iend-1;++j)
                {
                    if(rc->m_xi[0][j+1]<rc->m_xi[0][j])
                    {
                        tmpdouble = rc->m_xi[0][j+1];
                        rc->m_xi[0][j+1] = rc->m_xi[0][j];
                        rc->m_xi[0][j] = tmpdouble;

                        tmpdouble = rc->m_xi[1][j+1];
                        rc->m_xi[1][j+1] = rc->m_xi[0][j];
                        rc->m_xi[1][j] = tmpdouble;
                    }
                }
            }
            sum += numepoints;

            // bubble sort for second edge
            istart = numvert + sum;
            iend = istart + numepoints;
            for(i=istart;i<iend;++i)
            {
                for(j=istart;j<iend-1;++j)
                {
                    if(rc->m_xi[0][j+1]>rc->m_xi[0][j])
                    {
                        tmpdouble = rc->m_xi[0][j+1];
                        rc->m_xi[0][j+1] = rc->m_xi[0][j];
                        rc->m_xi[0][j] = tmpdouble;

                        tmpdouble = rc->m_xi[1][j+1];
                        rc->m_xi[1][j+1] = rc->m_xi[0][j];
                        rc->m_xi[1][j] = tmpdouble;
                    }
                }
            }
            sum += numepoints;

            // bubble sort for third edge
            istart = numvert + sum;
            iend = istart + numepoints;
            for(i=istart;i<iend;++i)
            {
                for(j=istart;j<iend-1;++j)
                {
                    if(rc->m_xi[1][j+1]>rc->m_xi[1][j])
                    {
                        tmpdouble = rc->m_xi[0][j+1];
                        rc->m_xi[0][j+1] = rc->m_xi[0][j];
                        rc->m_xi[0][j] = tmpdouble;

                        tmpdouble = rc->m_xi[1][j+1];
                        rc->m_xi[1][j+1] = rc->m_xi[0][j];
                        rc->m_xi[1][j] = tmpdouble;
                    }
                }
            }
            
            return;

        }        

        void NodalBasisManager::NodalPointExpand3D(BaryNodeContainer * rc)
        {
            ASSERTL0(false, "3D Point Expansion Not Implemented Yet");

            //    npts = (porder+1)*(porder+2)*(porder+3)/6;
            // 
            //    rc->m_xi[0][sum] = 2.0*b - 1.0;
            //	  rc->m_xi[1][sum] = 2.0*c - 1.0;
            //    rc->m_xi[2][sum] = 2.0*d - 1.0;
            //

        }


        bool operator  == (const BaryNodeContainer* x, const BaryNodeContainer &y)
        {
            ASSERTL0(x->m_ntype == y.m_ntype, "Nodal Basis type not the same");

            if(x->m_porder == y.m_porder)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator  == (const BaryNodeContainer& x, const BaryNodeContainer *y)
        {
            ASSERTL0(x.m_ntype == y->m_ntype, "Nodal Basis type not the same");

            if(x.m_porder == y->m_porder)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator  != (const BaryNodeContainer* x, const BaryNodeContainer &y)
        {

            ASSERTL0(x->m_ntype == y.m_ntype, "Nodal Basis type not the same");

            if(x->m_porder == y.m_porder)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator  != (const BaryNodeContainer& x, const BaryNodeContainer *y)
        {

            ASSERTL0(x.m_ntype == y->m_ntype, "Nodal Basis type not the same");

            if(x.m_porder == y->m_porder)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    }//end namespace
}//end namespace

/** 
* $Log: NodalBasisManager.cpp,v $
* Revision 1.6  2006/09/28 20:33:53  kirby
* *** empty log message ***
*
* Revision 1.5  2006/09/10 02:22:32  kirby
* *** empty log message ***
*
* Revision 1.4  2006/09/10 02:20:52  kirby
* *** empty log message ***
*
* Revision 1.3  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.2  2006/06/01 13:43:19  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:30  kirby
* *** empty log message ***
*
* Revision 1.18  2006/04/25 20:23:32  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.17  2006/04/01 21:59:26  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.16  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.15  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.14  2006/03/02 16:28:34  sherwin
*
* Corrected memory error in GenMassMatrix
*
* Revision 1.13  2006/02/19 13:26:12  sherwin
*
* Coding standard revisions so that libraries compile
*
*
**/ 


