///////////////////////////////////////////////////////////////////////////////
//
// File NodalBasisManager.h
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
// Description: Header file of manager for Nodal basis
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALBASISMANAGER_H
#define NODALBASISMANAGER_H

#include <StdRegions/PolyManager.h>
#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
  namespace StdRegions
  {

    static const int perm3A_2d[3][3] = {{0,1,2},{0,2,1},{2,0,1}};
    static const int perm3B_2d[3][3] = {{0,1,2},{1,0,2},{1,2,0}};
    static const int perm6_2d [6][3] = {{0,1,2},{0,2,1},{1,2,0},
					{1,0,2},{2,0,1},{2,1,0}};

    
    class BaryNodeContainer
    {
    public:
      /////////////////////////
      // Container Data
      /////////////////////////
      NodalBasisType m_ntype;     //< Nodal Point Type
      int m_dim;                  //< Dimension of point set
      int m_porder;               //< Polynomial order of point set
      /// Number of barycentric points (with no symmetry duplications)
      int m_bnpts;               
      /// Up to five possible symmetries (in 3D)
      int *m_symm[(int) NekConstants::kMaxSym];  
      /// Up to four barycentric coordinates (in 3D)
      double * m_bpoints[NekConstants::kMaxBary];  

      int m_npts;
      double * m_xi[NekConstants::kMaxDim];

      /////////////////////////////
      // Container Initialization
      // and destruction
      ////////////////////////////

      BaryNodeContainer():
	m_dim(0),
	m_porder(0),
	m_npts(0)
      {
	int i;
	for(i=0;i<NekConstants::kMaxSym;i++)
	{
	  m_symm[i] = (int*) NULL;
	}

	for(i=0;i<NekConstants::kMaxBary;i++)
	{
	  m_bpoints[i] = (double*)NULL;
	}

	for(i=0;i<NekConstants::kMaxDim;i++)
	{
	  m_xi[i] = (double*)NULL;
	}
      }

      ~BaryNodeContainer()
      {
	int i;
	for(i=0;i<NekConstants::kMaxSym;i++)
	{
	  if(m_symm[i])
	  {
	    delete[] m_symm[i];
	  }
	}

	for(i=0;i<NekConstants::kMaxBary;i++)
	{
	  if(m_bpoints[i])
	  {
	    delete[] m_bpoints[i];
	  }
	}

	for(i=0;i<NekConstants::kMaxDim;i++)
	{
	  if(m_xi[i])
	  {
	    delete[] m_xi[i];
	  }
	}
      }

      void set_ntype(NodalBasisType ntype)
      {
	m_ntype = ntype;
      }

      //Overloaded Operators
      friend bool operator  == (const BaryNodeContainer *x,
				const BaryNodeContainer &y);
      friend bool operator  == (const BaryNodeContainer& x,
				const BaryNodeContainer *y);
      friend bool operator  != (const BaryNodeContainer* x,
				const BaryNodeContainer &y);
      friend bool operator  != (const BaryNodeContainer& x,
				const BaryNodeContainer *y);
    };


    class NodalBasisManager: public PolyManager
    {
    public:
      NodalBasisManager();
      ~NodalBasisManager();

      int GetNodePoints(const NodalBasisType ntype, const int order,
			const double* &x, const double* &y, const double* &z);

    protected:

      /// Functions used by constructor to read nodal points From files
      void NodalFileRead2D(std::FILE *fp, BaryNodeContainer * rc);
      void NodalFileRead3D(std::FILE *fp, BaryNodeContainer * rc);
      
      /// Functions used to expand barycenter points to standard elements
      void NodalPointExpand2D(BaryNodeContainer * rc);
      void NodalPointExpand3D(BaryNodeContainer * rc);

    private:
      std::vector<BaryNodeContainer*> *m_nodecontainer;
    };

  } //end of namespace
} //end of namespace

#endif //NODALBASISMANAGER_H


/**
 * $Log: NodalBasisManager.h,v $
 * Revision 1.9  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.8  2006/03/06 12:39:59  sherwin
 *
 * Added NekConstants class for all constants in this library
 *
 * Revision 1.7  2006/03/05 22:11:02  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.6  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.5  2006/02/19 13:26:12  sherwin
 *
 * Coding standard revisions so that libraries compile
 *
 *
 **/
