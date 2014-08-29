#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sys/dir.h>
#include <cmath>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <SolverUtils/SolverUtilsDeclspec.h>


using namespace Nektar;

NekDouble m_Re;
NekDouble m_Mach;
NekDouble L;
NekDouble m_xo;
NekDouble m_Tinf;
NekDouble m_Suth;
NekDouble m_Tw;
NekDouble m_Twall;
NekDouble m_Gamma;
NekDouble m_Pr;
NekDouble m_longeur;
NekDouble m_uInf;
NekDouble m_rhoInf;
NekDouble m_R;
NekDouble m_vInf;
NekDouble m_To=273.11;
NekDouble m_mu;

int m_xpoints;
int m_ypoints;

int m_option;
int m_option2;

const NekDouble Nvisc = 1;
const NekDouble Omega = 1;

const NekDouble etamax = 10.0;
NekDouble errtol       = 1e-5;

	
/**
 * Calculate the compressible boundary layer using the similarity 
 * solution
 */
void COMPBL(
	Array<OneD, NekDouble> v, 
	Array<OneD, NekDouble> dv)
{	
	NekDouble c, dcdg, cp;

	if (Nvisc == 1) 
    {
		c = sqrt(v[3]) * (1.0 + m_Suth) / (v[3] + m_Suth);
		dcdg = 1.0 / (2. * sqrt(v[3])) - sqrt(v[3]) / (v[3]+m_Suth);
		dcdg = dcdg * (1.0 + m_Suth) / (v[3] + m_Suth);
		cp = dcdg * v[4];
	}
	if (Nvisc == 2)
	{
		c = pow(v[3], (Omega-1.0));
		dcdg = (Omega - 1.0) * pow(v[3], (Omega - 2.0));
		cp = dcdg * v[4];
    }
	if (Nvisc == 3)
	{
		c = sqrt(m_Twall) * (1.0 + m_Suth) / (m_Suth + m_Twall);
		cp = 0.0;
	}

		dv[0] = v[1];
		dv[1] = v[2];
		dv[2] = - v[2] * (cp + v[0]) / c;
		dv[3] = v[4];
		dv[4] = - v[4] * (cp + m_Pr * v[0]) / c - 
				m_Pr * (m_Gamma - 1.0) * pow(m_Mach, 2.0) * 
				pow(v[2], 2);
}

/**
 * Perform the RK4 integration
 */	
void RK4(
	Array<OneD, NekDouble> y, 
	Array<OneD, NekDouble> dydx, 
	int n, 
	NekDouble x, 
	NekDouble h, 
	Array<OneD, NekDouble> yout)
{	
	int nmax = 5;
		
	Array<OneD, NekDouble> yt(nmax, 0.0);
	Array<OneD, NekDouble> dyt(nmax, 0.0);
	Array<OneD, NekDouble> dym(nmax, 0.0);
	NekDouble hh = h * 0.5;
	NekDouble h6 = h / 6;
			
	for (int i = 0; i < n ; i++)
	{
		yt[i] = y[i] + hh * dydx[i];
	}
	
	COMPBL(yt, dyt);
	
	for (int i = 0; i < n; i++)
	{
		yt[i] = y[i] + hh * dyt[i];
	}
			
	COMPBL(yt, dym);
			
	for (int i = 0; i < n; i++)
	{
		yt[i] = y[i] + h * dym[i];
		dym[i] = dyt[i] + dym[i];
	}
		
	COMPBL(yt, dyt);
			
	for (int i = 0; i < n; i++)
	{
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2 * dym[i]);
	}
}

	
/**
 * Calculate ....
 */	
void RKDUMB(
	Array<OneD, NekDouble> vstart, 
	int nvar, 
	NekDouble x1, 
	NekDouble x2, 
	int m_xpoints , 
	Array<OneD, NekDouble> xx,
	Array<OneD, Array<OneD, NekDouble> > y)
{
	int nmax = 5;
	NekDouble x, h;
	Array<OneD, NekDouble> v (nmax, 0.0);
	Array<OneD, NekDouble> dv(nmax, 0.0);

	for (int i = 0; i < nvar; i++)
	{
		v[i] 	= vstart[i];
		y[i][0] = v[i];	
	}
		
	xx[0] = x1;
	x 	  = x1;
	h 	  = (x2-x1) / m_xpoints; 
					
	for (int k = 0; k < m_xpoints; k++)
	{
		COMPBL(v, dv);
		RK4	  (v, dv, nvar, x, h, v);
		
		if (x + h == x)
		{
			cout << "bug" << endl;
		}
		x 		= x + h;
		xx[k+1] = x;
		
		for (int i = 0; i < nvar; i++)
		{
			y[i][k+1] = v[i];
		}
	}
}

/**
 * Creating the output file 
 */
void OUTPUT(
	int n, 
	Array <OneD, NekDouble > xx, 
	Array<OneD, Array<OneD, NekDouble> > ff, 
	int nQuadraturePts,
	Array <OneD, NekDouble > x_QuadraturePts, 
	Array <OneD, NekDouble > y_QuadraturePts,
	Array <OneD, NekDouble > u_QuadraturePts, 
	Array <OneD, NekDouble > v_QuadraturePts,
	Array <OneD, NekDouble > rho_QuadraturePts, 
	Array <OneD, NekDouble > T_QuadraturePts)
{
	cout << "jusque la ca marche1 " << endl;
	
	Array <OneD, NekDouble > z  (m_xpoints, 0.0);
	Array <OneD, NekDouble > v  (m_xpoints, 0.0);
	Array <OneD, NekDouble > dv (m_xpoints, 0.0);
	Array <OneD, NekDouble > y  (m_xpoints, 0.0);
	Array <OneD, NekDouble > u  (m_xpoints, 0.0);
	Array <OneD, NekDouble > t  (m_xpoints, 0.0);
	Array <OneD, NekDouble > rho(m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > fx (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > fe (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > du (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > d2u(m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > dt (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > d2t(m_xpoints, 0.0);
	Array <OneD, NekDouble > mu (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > dmu(m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > y2 (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > u2 (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > t2 (m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > fx2(m_xpoints, 0.0);
	//~ Array <OneD, NekDouble > fe2(m_xpoints, 0.0);
	Array <OneD, NekDouble > vv(m_xpoints, 0.0);

	Array <OneD, Array <OneD, NekDouble > > valeur_u  (m_xpoints);
	Array <OneD, Array <OneD, NekDouble > > valeur_v  (m_xpoints);
	Array <OneD, Array <OneD, NekDouble > > valeur_rho(m_xpoints);
	//Array <OneD, Array <OneD, NekDouble > > valeur_t  (m_xpoints);
	
	cout << "jusque la ca marche2 " << endl;
		
	for (int i = 0; i < m_xpoints ; i++)
	{
		valeur_u[i]   = Array <OneD, NekDouble > (m_ypoints, 0.0);	
	}
	
	cout << "jusque la ca marche3 " << endl;
	
	for (int i = 0; i < m_xpoints ; i++)
	{
	valeur_v[i]   = Array <OneD, NekDouble > (m_ypoints, 0.0);	
	}
	
	cout << "jusque la ca marche4 " << endl;
	
	//~ for (int i = 0; i < m_xpoints ; i++)
	//~ {
	//~ valeur_t[i]   = Array <OneD, NekDouble > (m_ypoints, 0.0);	
	//~ }
		//~ 
	//~ cout << "jusque la ca marche6 " << endl;
	
	for (int i = 0; i < m_xpoints ; i++)
	{
	valeur_rho[i] = Array <OneD, NekDouble > (m_ypoints, 0.0);		
	}
	
	cout << "jusque la ca marche5 " << endl;
	
	
			
	
	
	
		
	NekDouble flg, scale, dd, dm, delta , marois, DELT;
	NekDouble xin, rex, delsx, xcher, ycher, dlta;
	NekDouble inter, inter2,sum, approx;
			
	int indice, indice2, indice3,compteur,compteur2,compteur3;

	z[0] = 0.0;
	NekDouble sumd = 0.0;
	

	
		
	for (int i=1; i < n ; i++)
	{
		z[i] = z[i-1] + 
				0.5 * (xx[i] - xx[i-1]) * (ff[3][i] + ff[3][i-1]);
		dm = ff[3][i-1] - ff[1][i-1];
		dd = ff[3][i] - ff[1][i];
		sumd = sumd + 0.5 * (xx[i] - xx[i-1]) * (dd + dm);
		
		if ((ff[1][i] > 0.999) && (flg < 1.0))
		{
			dlta = z[i];
			flg	 = 2.0;
		}
	}
	
	
		
	scale = sumd;
		
	ofstream file1;
	file1.open("test_u.dat");
			
	ofstream file2;
	file2.open("test_vv.dat");
			
	ofstream file3;
	file3.open("test_rho.dat");
			
	ofstream file4;
	file4.open("test_t.dat");	
		
	DELT = 4.91 * (sqrt(m_longeur)) / (sqrt(m_Re / (m_longeur)));
	cout << DELT << "   " << "delta" << endl;
		
	if (m_option == 1)
	{
		L = x_QuadraturePts[nQuadraturePts - 1] - m_xo;
	}
	else if (m_option == 2)
	{
		m_xo = x_QuadraturePts[0];
		L 	 = x_QuadraturePts[nQuadraturePts-1] - m_xo;
	}
		
	// X initialisation
	for (int kk = 0; kk < m_xpoints; kk++)
	{
		scale = sumd;
		xin   = m_xo + (L / (m_xpoints - 1)) * (kk + 1);
			
		cout << xin << endl;
		rex   = 0.5 * pow((m_Re / scale), 2) + m_Re * xin;
		delsx = sqrt(2.0 / rex) * scale * (xin)* m_Pr;
		scale = scale / delsx;
		delta = 4.91 * (sqrt(xin)) / (sqrt(m_Re / (m_longeur)));
		
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < 5; k++)
			{
				v[k] = ff[k][i];
			}		
			COMPBL(v, dv);
			y[i]  = (z[i] / scale) / delta;
			u[i]  = ff[1][i];
			t[i]  = ff[3][i];
			rho[i] = (1.0 / ff[3][i]);
			vv[i]=  -ff[0][i]/sqrt(m_uInf);
			
			mu[i]  = pow(t[i], 1.5) * 
						(1 + m_Suth) / (t[i] + m_Suth) / m_Re;
		}		

		for (int j = 0; j < m_ypoints; j++)
		{
			marois = 10 * DELT / ((m_ypoints - 1)) * j;
							
			for (int k = 0; k<n-2; k++)
			{			
				if ((y[k] <= marois)&&(y[k+1] > marois))
				{
					indice = k;
				}
			}				
			valeur_u[kk][j]   = (u[indice] + u[indice+1]) / 2;
			valeur_rho[kk][j] = (rho[indice] + rho[indice+1]) / 2;
			//valeur_t[kk][j]   = (t[indice] + t[indice+1]) / 2;
			
		}						
	}
	
	if (m_option2==0)
	{
		for (int i = 0; i < m_xpoints; i++)
		{
			
			for (int j = 0; j < m_ypoints; j++)
			{	
				valeur_v[i][j]=0;
			}
		}
	}
	
	if (m_option2==1)
	{
		
		for (int i = 0; i < m_xpoints-1; i++)
		{
			valeur_v[i][0]=0;
			for (int j = 0; j < m_ypoints-1; j++)
			{	
				valeur_v[i][j+1] = (1.0 / (valeur_rho[i][j+1])) *
				(valeur_v[i][j] * valeur_rho[i][j]) - 
			((10*DELT / m_ypoints) / (L / m_xpoints)) * (valeur_u[i+1][j] * 
			valeur_rho[i+1][j] - valeur_u[i][j] * valeur_rho[i][j]) ;
			}
		}
	}
	else
	{		
		
	
		for (int i = 0; i < m_xpoints-1; i+=int(m_xpoints/10))
		{
			valeur_v[i][0]=0;
			for (int j = 0; j < m_ypoints-1; j++)
			{	
			valeur_v[i][j+1] = (1.0 / (valeur_rho[i][j+1])) *
			(valeur_v[i][j] * valeur_rho[i][j]) - 
			((10*DELT / m_ypoints) / (L / m_xpoints)) * (valeur_u[i+1][j] * 
			valeur_rho[i+1][j] - valeur_u[i][j] * valeur_rho[i][j]) ;
			}
		}
					
			for (int j = 0; j < m_ypoints; j++)
			{	
			approx= 0.1 * (
	valeur_v[0][j] + valeur_v[int(m_xpoints/10)][j] +
	valeur_v[int(m_xpoints/10)*2][j] + valeur_v[int(m_xpoints/10)*3][j] +
	valeur_v[int(m_xpoints/10)*4][j] + valeur_v[int(m_xpoints/10)*5][j] +
	valeur_v[int(m_xpoints/10)*6][j] + valeur_v[int(m_xpoints/10)*7][j] +
	valeur_v[int(m_xpoints/10)*8][j] + valeur_v[int(m_xpoints/10)*9][j] );
	
				for (int i = 0; i < m_xpoints; i++)
				{
					valeur_v[i][j]= approx;
				}
		}
				

}




	
	for (int i = 0; i < nQuadraturePts; i++)
	{
		if (i%100000 == 0)
		{
			cout << "i" << "  " << i << "/" << nQuadraturePts << endl;
		}
	
		xcher = x_QuadraturePts[i];
		ycher = y_QuadraturePts[i];
					
		for (int j = 0; j < m_xpoints-1; j++)
		{							
			if ((m_xo + (L / (m_xpoints-1)) * j<=xcher) && (m_xo + 
				(L / (m_xpoints-1)) * (j+1)>xcher))
			{
			indice2 = j;
			break;					
			}
		}	
				
		//~ file4 << i  << "    "  << indice2 <<  "     "   
			  //~ << x_QuadraturePts[i] << "    \t" ;
		//~ file4 <<  m_xo + (L / (m_xpoints-1)) * indice2 << "       "  ;
					
		if (ycher >= 10 * DELT)
		{				
			u_QuadraturePts[i] = valeur_u[indice2][m_ypoints-10];
			v_QuadraturePts[i] = valeur_v[indice2][m_ypoints-10];						
			rho_QuadraturePts[i] = valeur_rho[indice2][m_ypoints-10];
			//T_QuadraturePts[i] = valeur_t[indice2][m_ypoints-10];	
			T_QuadraturePts[i] = 1.0/rho_QuadraturePts[i];	

			//~ file4 << "    "  << endl;				
		}
		else
		{			
			for (int j = 0; j < m_ypoints; j++)
			{				
				if (((10 * DELT / (m_ypoints-1)) * j <= ycher)&&
					((10 * DELT / (m_ypoints-1)) * (j+1) > ycher))
				{
					indice3 = j;
				break;					
				}
			}
														
			u_QuadraturePts[i] = valeur_u[indice2][indice3];
			v_QuadraturePts[i] = valeur_v[indice2][indice3];
			rho_QuadraturePts[i] = valeur_rho[indice2][indice3];
			//T_QuadraturePts[i]   = valeur_t[indice2][indice3];	
			T_QuadraturePts[i] = 1.0/rho_QuadraturePts[i];
			
			//~ file4 << y_QuadraturePts[i] << "  \t  "  
					//~ << (10 * DELT / (m_ypoints-1))  * indice3  << endl;	
			//~ file4 << "    "  << endl;			
		}
	}
		
	for (int j = 0; j < nQuadraturePts; j++)
	{		
		if ((x_QuadraturePts[j] <= 0.05) &&
			(x_QuadraturePts[j+1] > 0.05))
		{
			file2 << y_QuadraturePts[j] << "     "  
				  << v_QuadraturePts[j] << endl;	
		}				
	} 
	
	
	cout << "compteur v" << "   "  << compteur << endl;
	cout << "compteur u" << "   "  << compteur2 << endl;
	cout << "compteur rho" << "   "  << compteur3 << endl;
}














int main(int argc, char *argv[])
{
	// Variable initialisation
	int nmax  = 5;
	int maxit = 10;

	int  i, j, numModes;
	
	Array<OneD, NekDouble> 				 xx(m_xpoints, 0.0);
	Array<OneD, Array<OneD, NekDouble> > ff(nmax);
	Array<OneD, NekDouble>               parameter(9, 0.0);
	
	for (int i = 0; i < nmax; i++)
	{
		ff[i] = Array<OneD, NekDouble> (m_xpoints);	
	}
	
	Array<OneD, NekDouble > vstart(nmax, 0.0);
	Array<OneD, NekDouble > v(2);
	Array<OneD, NekDouble > dv(2);
	Array<OneD, NekDouble > f(2);
	Array<OneD, NekDouble > f1(2);
	Array<OneD, NekDouble > f2(2);
		
	NekDouble al11, al21, al12, al22, det;
	string opt;
	
	// Reading the session file
    LibUtilities::SessionReaderSharedPtr vSession = 
		LibUtilities::SessionReader::CreateInstance(argc, argv);	
	// Read in mesh from input file and create an object 
	// of class MeshGraph2D
    SpatialDomains::MeshGraphSharedPtr graphShPt; 
    graphShPt = MemoryManager<SpatialDomains::MeshGraph2D>
		::AllocateSharedPtr(vSession);
	//  Feed our spatial discretisation object
    MultiRegions::ContField2DSharedPtr Domain;
    Domain = MemoryManager<MultiRegions::ContField2D>
		::AllocateSharedPtr(vSession, graphShPt, 
							vSession->GetVariable(0));
							
    // Get the total number of elements
	int nElements;
	nElements = Domain->GetExpSize();
    std::cout << "Number of elements                 = " 
			  << nElements << std::endl;

    // Get the total number of quadrature points (depends on n. modes)
	int nQuadraturePts;
	int nQuadraturePts2;
    nQuadraturePts = Domain->GetTotPoints();
    nQuadraturePts2 = Domain->GetTotPoints();
    std::cout << "Number of quadrature points        = " 
			  << nQuadraturePts << std::endl;	
			  
	// Coordinates of the quadrature points
	Array<OneD,NekDouble> x_QuadraturePts;
	Array<OneD,NekDouble> y_QuadraturePts;
    Array<OneD,NekDouble> x_QuadraturePts2;
    Array<OneD,NekDouble> y_QuadraturePts2;
    Array<OneD,NekDouble> z_QuadraturePts;
    x_QuadraturePts  = Array<OneD,NekDouble>(nQuadraturePts);
    y_QuadraturePts  = Array<OneD,NekDouble>(nQuadraturePts);
    x_QuadraturePts2 = Array<OneD,NekDouble>(nQuadraturePts2);
    y_QuadraturePts2 = Array<OneD,NekDouble>(nQuadraturePts2);
    z_QuadraturePts  = Array<OneD,NekDouble>(nQuadraturePts);
    Domain->GetCoords(x_QuadraturePts, 
					  y_QuadraturePts, 
					  z_QuadraturePts);
    Domain->GetCoords(x_QuadraturePts2, 
					  y_QuadraturePts2, 
                      z_QuadraturePts);
    
	vSession->LoadParameter("Re",			m_Re,		1.0);
	vSession->LoadParameter("Mach",			m_Mach,		1.0);
	vSession->LoadParameter("xo",			m_xo,		1.0);
	vSession->LoadParameter("TInf", 		m_Tinf,		1.0);
	vSession->LoadParameter("Twall",		m_Twall,	1.0);
	vSession->LoadParameter("Gamma",		m_Gamma,	1.0);
	vSession->LoadParameter("Pr",			m_Pr,	 	1.0);
	vSession->LoadParameter("L",			m_longeur,  1.0);
	vSession->LoadParameter("rhoInf",		m_rhoInf,   1.0);
	vSession->LoadParameter("uInf",			m_uInf,		1.0);
	vSession->LoadParameter("GasConstant",	m_R,		1.0);
	vSession->LoadParameter("vInf",			m_vInf,		1.0);
	vSession->LoadParameter("MESH_USE",		m_option,	1);
	vSession->LoadParameter("mu",			m_mu,		1.0);
	vSession->LoadParameter("v_calcul",		m_option2,	1);
	vSession->LoadParameter("xpoints",		m_xpoints,	1);
	vSession->LoadParameter("ypoints",		m_ypoints,	1);
	
	m_Re 	= m_Re / m_longeur;
	m_Suth 	= 110.4 / m_Tinf;
	m_Tw 	= m_Twall / m_Tinf;
	
	cout << "option" << "   " << m_option << endl;

	Array<OneD,NekDouble> u_QuadraturePts;
	u_QuadraturePts = Array<OneD,NekDouble>  (nQuadraturePts, 0.0);
	Array<OneD,NekDouble> v_QuadraturePts;
	v_QuadraturePts = Array<OneD,NekDouble>  (nQuadraturePts, 0.0);
	Array<OneD,NekDouble> rho_QuadraturePts;
	rho_QuadraturePts = Array<OneD,NekDouble>(nQuadraturePts, 0.0);
	Array<OneD,NekDouble> T_QuadraturePts;
	T_QuadraturePts = Array<OneD,NekDouble>	 (nQuadraturePts, 0.0);
	
	Array<OneD,NekDouble> u_QuadraturePts2;
	u_QuadraturePts2 = Array<OneD,NekDouble>  (nQuadraturePts2, 0.0);
	Array<OneD,NekDouble> v_QuadraturePts2;
	v_QuadraturePts2 = Array<OneD,NekDouble>  (nQuadraturePts2, 0.0);
	Array<OneD,NekDouble> rho_QuadraturePts2;
	rho_QuadraturePts2 = Array<OneD,NekDouble>(nQuadraturePts2, 0.0);
	Array<OneD,NekDouble> T_QuadraturePts2;
	T_QuadraturePts2 = Array<OneD,NekDouble>  (nQuadraturePts2, 0.0);

	Array<OneD, MultiRegions::ExpListSharedPtr> Exp(4);

	if(m_Tw > 0) 
	{
		vstart[3] = m_Tw;
	}
    if (m_Tw < 0.0)
    {
		v[1] = 1.0 + 0.5 * 0.84 * (m_Gamma- 1) * (m_Mach * m_Mach);    
		v[0] = 0.47 * pow(v[1], 0.21);
    }
    else
    {
		v[1] = 0.062 * pow(m_Mach, 2) - 0.1 * (m_Tw - 1.0) * 
				(10 + m_Mach) / (0.2 + m_Mach);
		v[0] = 0.45 - 0.01 * m_Mach + (m_Tw - 1.0) * 0.06;
		m_Twall = m_Tw;
	}

	dv[0] = v[0] * 0.01;
	
	if (m_Tw < 0.0)
	{
		dv[1] = v[1] * 0.01;
	}
	else
	{
		dv[1] = 0.1;
	}      
		
	vstart[2] = v[0];
		
	if (m_Tw < 0)
	{
		vstart[3] = v[1];
		m_Twall = vstart[3];
	}
	else
	{
		vstart[4] = v[1];
	}
	
	RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);
	
	for (int k = 0; k < maxit; k++)
	{
		vstart[2] = v[0];
	
		if (m_Tw < 0)
		{
			vstart[3] = v[1];
			m_Twall   = vstart[3];
		}
		else
		{
			vstart[4] = v[1];
		}

		RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);
	
		NekDouble err = fabs(ff[1][m_xpoints] - 1) + 
						fabs(ff[3][m_xpoints] - 1);
						
		cout << "err" << scientific 
			 << setprecision(9) 
			 << "   " << err << endl;
	
		if (err < errtol)
		{
			cout << "fini" << endl;
			OUTPUT(m_xpoints, xx, ff, 
				   nQuadraturePts, 
				   x_QuadraturePts, 
				   y_QuadraturePts, 
				   u_QuadraturePts, 
				   v_QuadraturePts, 
				   rho_QuadraturePts, 
				   T_QuadraturePts);
			break;
		}
		else
		{
			f[0] = ff[1][m_xpoints] - 1;
			f[1] = ff[3][m_xpoints] - 1;
			vstart[2] = v[0] + dv[0];
			
			if (m_Tw < 0)
			{
				vstart[3] = v[1];
				m_Twall = vstart[3];
			}
			else
			{
				vstart[4] = v[1];
			}
			
			RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);
			
			f1[0] = ff[1][m_xpoints] - 1;
			f1[1] = ff[3][m_xpoints] - 1;
			
			vstart[2] = v[0];
			
			if (m_Tw < 0)
			{
				vstart[3] = v[1] + dv[1];
				m_Twall = vstart[3];
			}
			else
			{
				vstart[4] = v[1] + dv[1];
			}
		
			RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);
			
			f2[0] = ff[1][m_xpoints] - 1;
			f2[1] = ff[3][m_xpoints] - 1;
			
			al11 = (f1[0] - f[0]) / dv[0];
			al21 = (f1[1] - f[1]) / dv[0];
			al12 = (f2[0] - f[0]) / dv[1];
			al22 = (f2[1] - f[1]) / dv[1];
			det = al11 * al22 - al21 * al12;
			
			dv[0] = ( - al22 * f[0] + al12 * f[1]) / det;
			dv[1] = (al21 * f[0] - al11 * f[1]) / det;
			v[0]  = v[0] + dv[0];
			v[1]  = v[1] + dv[1];
		}
	}
	
	ofstream file15;
	file15.open("fin_verif2.dat");
			
	for (int i=0; i< nQuadraturePts; i++)
	{
		file15 << scientific << setprecision(9) << x_QuadraturePts[i] 
				<< "  \t  " << y_QuadraturePts[i] << "  \t " ;
		file15 << scientific << setprecision(9) << u_QuadraturePts[i] 
				<< "  \t  " << v_QuadraturePts[i] << "  \t " ;
		file15 << scientific << setprecision(9) << rho_QuadraturePts[i] 
				<< "  \t  " << T_QuadraturePts[i] << endl;
	}	
			
	
			
	for (int i = 0; i < nQuadraturePts; i++)
	{
		rho_QuadraturePts[i] = rho_QuadraturePts[i] * m_rhoInf;
		u_QuadraturePts[i] = u_QuadraturePts[i] * m_uInf;
		v_QuadraturePts[i] = v_QuadraturePts[i] * m_uInf;		
		T_QuadraturePts[i] = T_QuadraturePts[i] * m_Tinf;
								
		T_QuadraturePts[i] = T_QuadraturePts[i] * rho_QuadraturePts[i] 
							* m_R;
		T_QuadraturePts[i] = T_QuadraturePts[i] / (m_Gamma-1);
		T_QuadraturePts[i] = T_QuadraturePts[i] + 0.5 * 
							rho_QuadraturePts[i] * 
							(pow(u_QuadraturePts[i], 2.0) + 
							pow(v_QuadraturePts[i], 2.0));
					
		u_QuadraturePts[i] = u_QuadraturePts[i] * rho_QuadraturePts[i];
		v_QuadraturePts[i] = v_QuadraturePts[i] * rho_QuadraturePts[i];
	}

	MultiRegions::ExpList2DSharedPtr Exp2D_uk;
    Exp2D_uk = MemoryManager<MultiRegions::ExpList2D>::
		AllocateSharedPtr(vSession,graphShPt);
 
    MultiRegions::ExpList2DSharedPtr Exp2D_vk;
    Exp2D_vk = MemoryManager<MultiRegions::ExpList2D>::
		AllocateSharedPtr(vSession,graphShPt);
    
    MultiRegions::ExpList2DSharedPtr Exp2D_rhok;
    Exp2D_rhok = MemoryManager<MultiRegions::ExpList2D>::
		AllocateSharedPtr(vSession,graphShPt);
    
    MultiRegions::ExpList2DSharedPtr Exp2D_Tk;
    Exp2D_Tk = MemoryManager<MultiRegions::ExpList2D>::
		AllocateSharedPtr(vSession,graphShPt);
    
    // Filling the 2D expansion using a recursive   
    // algorithm based on the mesh ordering
    LibUtilities::BasisSharedPtr Basis;
    Basis    = Domain->GetExp(0)->GetBasis(0);
    numModes = Basis->GetNumModes();

    std::cout<< "Number of modes = " << numModes << std::endl;
    
    // Copying the ukGlobal vector in m_phys 
    // (with the same pattern of m_phys) 
    Vmath::Vcopy(nQuadraturePts, u_QuadraturePts ,  1, 
		Exp2D_uk->UpdatePhys(), 1);
    Vmath::Vcopy(nQuadraturePts, v_QuadraturePts,   1, 
		Exp2D_vk->UpdatePhys(), 1);
    Vmath::Vcopy(nQuadraturePts, rho_QuadraturePts, 1, 
		Exp2D_rhok->UpdatePhys(), 1);
    Vmath::Vcopy(nQuadraturePts, T_QuadraturePts ,  1, 
		Exp2D_Tk->UpdatePhys(), 1);
    
	// Initialisation of the ExpList Exp
    Exp[0] = Exp2D_rhok;
    Exp[1] = Exp2D_uk;
    Exp[2] = Exp2D_vk;
    Exp[3] = Exp2D_Tk;
    
    // Expansion coefficient extraction 
    // (necessary to write the .fld file)    
    Exp[0]->FwdTrans(Exp2D_rhok->GetPhys(),Exp[0]->UpdateCoeffs());
    Exp[1]->FwdTrans(Exp2D_uk->GetPhys(),Exp[1]->UpdateCoeffs());
    Exp[2]->FwdTrans(Exp2D_vk->GetPhys(),Exp[2]->UpdateCoeffs());
    Exp[3]->FwdTrans(Exp2D_Tk->GetPhys(),Exp[3]->UpdateCoeffs());
    // Generation .FLD file with one field only (at the moment)
    // Definition of the name of the .fld file
    string FalknerSkan = "resultat_fin.fld";
    // Definition of the Field
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef = 
		Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
       
    for (j = 0; j < 4; j++)
    {
        for (i = 0; i < FieldDef.size(); i++)
        {
            if (j == 0)
            {
                FieldDef[i]->m_fields.push_back("rho");
            }
            else if (j == 1)
            {
                FieldDef[i]->m_fields.push_back("rhou");
            }
            else if (j == 2 )
            {
                FieldDef[i]->m_fields.push_back("rhov");
            }
            else if (j == 3 )
            {
                FieldDef[i]->m_fields.push_back("E");
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }    
    
    LibUtilities::Write(FalknerSkan, FieldDef, FieldData);

    std::cout<<"----------------------------------------------------\n";
	std::cout <<"\n=================================================\n";
    std::cout <<"Similarity solution \n"; 
    std::cout <<"===================================================\n";
    std::cout <<"***************************************************\n";
    std::cout <<"DATA FROM THE SESSION FILE:\n";
    std::cout << "Reynolds number                  = " << m_Re      
			  << "\t[-]"   << std::endl;
    std::cout << "Mach number                      = " << m_Mach               
			  << "\t[-]"   << std::endl;
    std::cout << "Characteristic length            = " << m_longeur                
			  << "\t\t[m]" << std::endl;
    std::cout << "U_infinity                       = " << m_uInf            
			  << "\t[m/s]" << std::endl;
    std::cout << "Position x_0 to start the BL [m] = " << m_xo              
			  << "\t\t[m]" << std::endl;
    std::cout <<"***************************************************\n";
    std::cout <<"---------------------------------------------------\n";
    std::cout <<"MESH and EXPANSION DATA:\n";
    std::cout << "Done." << std::endl;
						
    return 0;
}




