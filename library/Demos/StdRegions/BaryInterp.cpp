///////////////////////////////////////////////////////////////////////////////
//
// File: NodalDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/Timer.h>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

namespace po = boost::program_options;

template <class myType>
int commoncode(myType*E, 
	       NekMatrix<NekDouble> matT,
	       int n_coeffs)
{
  Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
  t[0] = -0.5;
  t[1] = -0.25;
  t[2] = -0.3;
   
  const Array<OneD, const NekDouble> tt = t;
  int tallyflag = 0;
  for(int uu = 0; uu<n_coeffs; uu++)
  { 	      
     
    // Using existing Nek++ version:
    Timer t1;
    t1.Start();   
    NekDouble val1 = E->PhysEvaluate(t,matT.GetPtr()+uu*(matT.GetRows()));
    t1.Stop();
    NekDouble elapsed  = t1.TimePerTest(1);
    cout<<"\n Original Nek++ PhysEvaluate took = "<<elapsed<<" s. val = "<<val1<<"\n";
      
    // Using (1): separable basis in x- and y- direction
    // case (1A)
    Array<OneD, Array<OneD,  NekDouble> >quadpts(2);
    quadpts[0] = E->GetBasis(0)->GetZ();
    quadpts[1] = E->GetBasis(1)->GetZ();
		
    Array<OneD, Array<OneD,  NekDouble> >indphis(2);
    indphis[0] = E->GetBasis(0)->GetBdata();
    indphis[1] = E->GetBasis(1)->GetBdata();
    t1.Start();   				    
    NekDouble val1A = E->PhysEvaluateBaryInd(tt, quadpts, indphis, 0);
    t1.Stop();
    NekDouble elapsed1A  = t1.TimePerTest(1);
    cout<<"\n 1A case: took "<<elapsed1A<<" s.val = "<<val1A;
      
    //case (1Ba)
    Array<OneD, Array<OneD,  NekDouble> >quadpts2(2);
    Array<OneD, Array<OneD,  NekDouble> >indphis2(2);
    t1.Start();   
    NekDouble val1Ba = E->PhysEvaluateBaryInd(tt, quadpts, indphis, uu);
    t1.Stop();
    NekDouble elapsed1Ba  = t1.TimePerTest(1);
    cout<<"\n 1Ba case: took "<<elapsed1Ba<<" s.val = "<<val1Ba;
    
    //case (1Bb) assert should fail
    t1.Start();   
    NekDouble val1Bb = E->PhysEvaluateBaryInd(tt, quadpts, indphis, n_coeffs+1);
    t1.Stop();
    NekDouble elapsed1Bb  = t1.TimePerTest(1);
    cout<<"\n 1Bb case: took "<<elapsed1Bb<<" s. val = "<<val1Bb<<"\n";

    // Using (2): Barycentriuc interp version of existing physevaluate()
    t1.Start();   
    NekDouble val2 = E->PhysEvaluateBary(tt, matT.GetPtr()+uu*(matT.GetRows()));
    t1.Stop();
    NekDouble elapsed2  = t1.TimePerTest(1);
    cout<<"\n 2 case: took "<<elapsed2<<" s.val = "<<val2;

    //check
    if(abs(val1 - val1A)>1e-9 ||
       abs(val1 - val1Ba)>1e-9 ||
       abs(val1 - val2)>1e-9  )
      {
	tallyflag = 1;
      }
  }
  return tallyflag;
}



// Given \phi_1(x) and \phi_2(x) return tensor prod if tensorp = 1
template <class myType>

void PolyEval2(                
	       myType *expobj,
	       NekMatrix< NekDouble > &Iprod, int m_ncoeffs 
			       )
{
  // if tensorp = 1 return tensor prod
  //          i.e let(z = [x_1, y_1]) 
  //          return \forall i \forall j [\phi(z_i)*\phi(z_j)]
  // else return \phi(x_1)*\phi(y_1)
  int  nquad0  = expobj->GetBasis(0)->GetNumPoints();
  int  nquad1  = expobj->GetBasis(1)->GetNumPoints();
  int  nmodes0 = expobj->GetBasis(0)->GetNumModes();
  //int  nmodes1 = expobj->GetBasis(1)->GetNumModes();
  NekMatrix<NekDouble> matX(nquad0*nquad1,m_ncoeffs );//nmodes0*nmodes1);
  ASSERTL0(nquad0 == nquad1,"failed in polyeval2" );
  int coltemp = 0;
  int ctrtemp = 0, tempi = 0;
  for(int i = 0; i < m_ncoeffs; i++, tempi++)//nmodes0*nmodes1; i++)
    {
	    
      Array<OneD, NekDouble> vals(nquad0*nquad1); 
	    
      expobj->FillMode(i,vals);

      // coltemp = (coltemp+(nmodes0-((tempi)%(nmodes0))+1));	        
      coltemp = (coltemp+(nmodes0-(tempi))+1);	        
      if(coltemp > m_ncoeffs-ctrtemp-1)
	{
	  ctrtemp = ctrtemp + 1;
	  tempi = 0;
	  coltemp = ctrtemp;

	}
      if(i == 0)
	coltemp = 0;



      //	cout<<"\n coltemp="<<coltemp<< " i = "<<i;
      //Vmath::Vcopy(nquad1*nquad0, &vals[0], 1, &(matX.GetPtr()[loc]),1);
      Vmath::Vcopy(nquad1*nquad0, &vals[0], 1, &(matX.GetPtr()[i*nquad0*nquad1]),1);
	      
    }
  Iprod=matX;
	

} 


int main(int argc, char *argv[])
{
  string shape, ntype;
  int baryinterp = 1;
  vector<string> basis(3, "NoBasisType"), pointstype(3, "NoPointsType");
  vector<int> order, points;
  po::options_description desc("Available options");
  desc.add_options()
    ("help,h",
     "Produce this help message and list basis and shape types.")
    ("nodal,n",
     po::value<string>(&ntype),
     "Optional nodal type, autofills shape and basis choices.")
    ("shape,s",
     po::value<string>(&shape),
     "Region shape to project function on.")
    ("basis,b",
     po::value<vector<string>>(&basis)->multitoken(),
     "Basis type, separate by spaces for higher dimensions.")
    ("order,o",
     po::value<vector<int>>(&order)->multitoken()->required(),
     "Order of basis sets, separate by spaces for higher dimensions.")
    ("points,p",
     po::value<vector<int>>(&points)->multitoken()->required(),
     "Number of quadrature points, separate by spaces for "
     "higher dimensions.")
    ("pointstype,P",
     po::value<vector<string>>(&pointstype)->multitoken(),
     "Optional points type, separate by spaces for higher dimensions.")
    ("diff,d",
     "Project derivative.");

  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, desc), vm);
      if (vm.count("help"))
        {
	  cout << desc;
	  cout << endl << "All nodal types, -n [ --nodal ], are:" << endl;
	  for (int i = 22; i < SIZE_PointsType; ++i)
            {
	      cout << kPointsTypeStr[i] << endl;
            };
	  cout << endl << "All shape types, -s [ --shape ], are:" << endl;
	  for (int i = 1; i < SIZE_ShapeType; ++i)
            {
	      cout << ShapeTypeMap[i] << endl;
            };
	  cout << endl << "All basis types, -b [ --basis ], are:" << endl;
	  for (int i = 1; i < SIZE_BasisType; ++i)
            {
	      cout << BasisTypeMap[i] << endl;
            };
	  cout << endl << "All points types, -P [ --pointstype ], are:"
	       << endl;
	  for (int i = 1; i < SIZE_PointsType; ++i)
            {
	      cout << kPointsTypeStr[i] << endl;
            };
	  return 1;
        }
      po::notify(vm);
    }
  catch (const exception &e)
    {
      cerr << "Error: " << e.what() << endl << desc;
      return 0;
    }

  vector<PointsType> ptype;
  if (vm.count("pointstype"))
    {
      for (auto &p : pointstype)
        {
	  PointsType tmp = eNoPointsType;
	  for (int i = 1; i < SIZE_PointsType; ++i) // starts at nodal points
            {
	      if (boost::iequals(kPointsTypeStr[i], p))
                {
		  tmp = static_cast<PointsType>(i);
		  break;
                }
	      ASSERTL0(i != SIZE_PointsType - 1,
		       "The points type '" + p + "' does not exist");
            }

	  ptype.push_back(tmp);
        }
    }

  //Convert string input argument to nodal type
  PointsType nodaltype = eNoPointsType;
  ShapeType stype = eNoShapeType;
  vector<BasisType> btype(3, eNoBasisType);
  if (vm.count("nodal"))
    {
      for (int i = 22; i < SIZE_PointsType; ++i) // starts at nodal points
        {
	  if (boost::iequals(kPointsTypeStr[i], ntype))
            {
	      nodaltype = static_cast<PointsType>(i);
	      break;
            }
	  ASSERTL0(i != SIZE_PointsType - 1, ("The nodal type '" + ntype +
					      "' does not exist"))
	    }
      switch (nodaltype)
        {
	case eNodalTriElec:
	case eNodalTriFekete:
	case eNodalTriSPI:
	case eNodalTriEvenlySpaced:
	  btype[0] = eOrtho_A;
	  btype[1] = eOrtho_B;
	  stype = eTriangle;
	  break;
	case eNodalQuadElec:
	  btype[0] = eOrtho_A;
	  btype[1] = eOrtho_B;
	  stype = eQuadrilateral;
	  break;
	case eNodalTetElec:
	case eNodalTetSPI:
	case eNodalTetEvenlySpaced:
	  btype[0] = eOrtho_A;
	  btype[1] = eOrtho_B;
	  btype[2] = eOrtho_C;
	  stype = eTetrahedron;
	  break;
	case eNodalPrismElec:
	case eNodalPrismSPI:
	case eNodalPrismEvenlySpaced:
	  btype[0] = eOrtho_A;
	  btype[1] = eOrtho_A;
	  btype[2] = eOrtho_B;
	  stype = ePrism;
	  break;
	case eNodalHexElec:
	  btype[0] = eOrtho_A;
	  btype[1] = eOrtho_A;
	  btype[2] = eOrtho_A;
	  stype = eHexahedron;
	  break;
	default:
	  ASSERTL0(!nodaltype, ("The nodal type '" + ntype +
				"' is invalid for StdProject."));
	  break;
        }
    }

  //Convert string input argument to shape type
  if (stype == eNoShapeType)
    {
      for (int i = 1; i < SIZE_ShapeType; ++i)
        {
	  if (boost::iequals(ShapeTypeMap[i], shape))
            {
	      stype = static_cast<ShapeType>(i);
	      break;
            }
	  ASSERTL0(i != SIZE_ShapeType - 1, ("The shape type '" + shape +
					     "' does not exist"))
	    }
    }

    if (stype == ePoint && vm.count("diff"))
    {
        NEKERROR(ErrorUtil::efatal,
                 "It is not possible to run the diff version for shape: point")
    }

    //Check arguments supplied equals dimension
    const int dimension = (stype == ePoint) ? 1 : ShapeTypeDimMap[stype];
    ASSERTL0(order.size() == dimension,
             "Number of orders supplied should match shape dimension");
    ASSERTL0(points.size() == dimension,
             "Number of points supplied should match shape dimension");
    ASSERTL0(ptype.size() == dimension || ptype.size() == 0,
             "Number of points types should match shape dimension if "
             "supplied.");

    if (!vm.count("nodal"))
    {
        ASSERTL0(basis.size() == dimension,
                 "Number of bases supplied should match shape dimension");
        //Convert string input argument to basis types
        for (int i = 0; i < dimension; ++i)
        {
            for (int j = 1; j < SIZE_BasisType; ++j)
            {
                if (boost::iequals(BasisTypeMap[j], basis[i]))
                {
                    btype[i] = static_cast<BasisType>(j);
                    break;
                }
                ASSERTL0(j != SIZE_BasisType - 1, ("The basis type '" + basis[i]
                                                   + "' does not exist"))
            }
        }
    }

  //check basis selection is permitted for chosen shape
  map<ShapeType, vector<vector<BasisType>>> allowableBasis;
  allowableBasis[ePoint] = {
    {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
     eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
     eFourierHalfModeRe, eFourierHalfModeIm}
  };
  allowableBasis[eSegment] = { allowableBasis[ePoint][0] };
  allowableBasis[eTriangle] = {
    {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
  };
  allowableBasis[eQuadrilateral] = {
    allowableBasis[eSegment][0], allowableBasis[eSegment][0]
  };
  allowableBasis[eTetrahedron] = {
    {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_C, eModified_C, eGLL_Lagrange, eGauss_Lagrange}
  };
  allowableBasis[ePyramid] = {
    {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
    {eOrthoPyr_C, eModifiedPyr_C, eGLL_Lagrange, eGauss_Lagrange}
  };
  allowableBasis[ePrism] = {
    {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
    {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
  };
  allowableBasis[eHexahedron] = {
    {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
     eLegendre, eChebyshev, eMonomial},
    {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
     eLegendre, eChebyshev, eMonomial},
    {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
     eLegendre, eChebyshev, eMonomial}
  };
    

  for (int i = 0; i < dimension; ++i)
    {
      const unsigned int basisListLength = allowableBasis[stype][i].size();
      for (int j = 0; j < basisListLength; ++j)
        {
	  if (allowableBasis[stype][i][j] == btype[i])
            {
	      break;
            }
	  ASSERTL0(j != basisListLength - 1,
		   ("The basis type '" +
		    static_cast<string>(BasisTypeMap[btype[i]]) +
		    "' is invalid for basis argument " + to_string(i + 1) +
		    " for shape '" + ShapeTypeMap[stype] + "'."))
	    }
    }

  //Declaration of other variables needed
  StdExpansion *E = nullptr;

  // Assign points type according to basis type selection, if not already
  // assigned.
  if (ptype.size() == 0)
    {
      ptype.resize(dimension);
      for (int i = 0; i < dimension; ++i)
        {
	  if (btype[i] == eFourier)
            {
	      ptype[i] = eFourierEvenlySpaced;
            }
	  else if (btype[i] == eFourierSingleMode ||
		   btype[i] == eFourierHalfModeRe ||
		   btype[i] == eFourierHalfModeIm)
            {
	      ptype[i] = eFourierSingleModeSpaced;
            }
	  else
            {
	      if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
                {
		  ptype[i] = eGaussRadauMAlpha1Beta0;
                }
	      else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
                {
		  ptype[i] = eGaussRadauMAlpha2Beta0;
                }
	      else if (i == 2 && stype == ePrism)
                {
		  ptype[i] = eGaussRadauMAlpha1Beta0;
                }
	      else
                {
		  ptype[i] = eGaussLobattoLegendre;
                }
            }
        }
    }

  vector<PointsKey> pkey;
  vector<BasisKey> bkey;
  for (int i = 0; i < dimension; ++i)
    {
      pkey.emplace_back(PointsKey(points[i], ptype[i]));
      bkey.emplace_back(BasisKey(btype[i], order[i], pkey[i]));
    }

  switch (stype)
    {
      cout<<"stype = "<<stype;
    case ePoint:
      {
	E = new StdPointExp(bkey[0]);
	break;
      }
    case eSegment:
      {
	E = new StdSegExp(bkey[0]);
	break;
      }
    case eTriangle:
      {
	E = nodaltype != eNoPointsType ? new StdNodalTriExp(bkey[0],
							    bkey[1],
							    nodaltype)
	  : new StdTriExp(bkey[0], bkey[1]);
	break;
      }
    case eQuadrilateral:
      {
	E = new StdQuadExp(bkey[0], bkey[1]);
	break;
      }
    case eTetrahedron:
      {
	E = nodaltype != eNoPointsType ? new StdNodalTetExp(bkey[0],
							    bkey[1],
							    bkey[2],
							    nodaltype)
	  : new StdTetExp(bkey[0], bkey[1],
			  bkey[2]);
	break;
      }
    case ePyramid:
      {
	E = new StdPyrExp(bkey[0], bkey[1], bkey[2]);
	break;
      }
    case ePrism:
      {
	E = nodaltype != eNoPointsType ? new StdNodalPrismExp(bkey[0],
							      bkey[1],
							      bkey[2],
							      nodaltype)
	  : new StdPrismExp(bkey[0], bkey[1],
			    bkey[2]);
	break;
      }
    case eHexahedron:
      {
	E = new StdHexExp(bkey[0], bkey[1], bkey[2]);
	break;
      }
    default:
      break;
    }

    const auto totPoints = (unsigned) E->GetTotPoints();
    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(totPoints);

    switch (dimension)
    {
        case 1:
        {
            E->GetCoords(x);
            break;
        }

        case 2:
        {
            E->GetCoords(x, y);
            break;
        }

        case 3:
        {
            E->GetCoords(x, y, z);
            break;
        }
        default:
            break;
    }


  // Array<OneD, NekDouble> custptval(n_coeffs),custptvalbary(n_coeffs);


  NekMatrix<NekDouble> matT;
  int n_coeffs;
  int tallyflag = 0;
  if(( strcmp(ShapeTypeMap[stype],"Triangle") == 0
      || strcmp(ShapeTypeMap[stype],"Quadrilateral") == 0 ) &&( baryinterp == 1))
    {
      if( strcmp(ShapeTypeMap[stype], "Triangle") == 0 )
	{
	  StdTriExp exp1(bkey[0],bkey[1]);
	  int nmodes0 = E->GetBasis(0)->GetNumModes();
	  int nmodes1 = E->GetBasis(1)->GetNumModes();

	  n_coeffs = LibUtilities::StdTriData::getNumberOfCoefficients(nmodes0,nmodes1);
	  PolyEval2(&exp1, matT, n_coeffs);
	  tallyflag = commoncode(&exp1,matT,n_coeffs);
	}
      else  if( strcmp(ShapeTypeMap[stype], "Quadrilateral") == 0 )
	{
	  StdQuadExp exp1(bkey[0],bkey[1]);

	  int nmodes0 = E->GetBasis(0)->GetNumModes();
	  int nmodes1 = E->GetBasis(1)->GetNumModes();

	  n_coeffs = LibUtilities::StdQuadData::getNumberOfCoefficients(nmodes0,nmodes1);
	  PolyEval2(&exp1, matT, n_coeffs);

	  tallyflag = commoncode(&exp1,matT,n_coeffs);		
	}
      else
	{
	  cout<<"\n Error! Please enter Triangle or Quads only";
	  exit(0);
	}
	    
    }
  if(tallyflag == 1)
    {
	cout<<"\n failed!";
	exit(0);
  }

	  
return 0;
}

