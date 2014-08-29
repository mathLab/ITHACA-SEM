    #include <cstdio>
    #include <cstdlib>
    #include <string>
    #include <iostream>
    #include <iomanip>
    #include <sys/dir.h>

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

    // #include <SolverUtils/Filters/Filter.h>

    using namespace Nektar;

    int main(int argc, char *argv[])
    {
		
	string fname1 = std::string(argv[1]);
	string fname2 = std::string(argv[2]);
	string fname3 = std::string(argv[3]);
	NekDouble x_in=atof(argv[4]);
	
	cout << fname1 << "  " << fname2 << "  " << fname3 << endl;
	
	vector <string> monTableau;
	int Nprof,Npoint,Iprof,indice,k;
	monTableau.push_back("");
	NekDouble rubbish;
	
	NekDouble delta, delta2, delta3, Re,c1;
	


// Premier fichier
// ---------------------------
// ---------------------------
// ---------------------------

			ofstream vitesse_u;
			vitesse_u.open("u.dat");
	
			ofstream vitesse_v;
			vitesse_v.open("v.dat");
			
			ofstream temperature;
			temperature.open("temperature.dat");
			
			ofstream test;
			test.open("test.dat");
			
			ofstream energie;
			energie.open("energy.dat");
			
			ofstream rho_1;
			rho_1.open("rho.dat");

			ofstream deltax;
			deltax.open("deltax.dat");

			std::ifstream fichier(argv[1], ios::in);
			
			getline(fichier,monTableau.back());
			fichier >> rubbish >> rubbish >> Re >> rubbish >> c1 ; 
			getline(fichier,monTableau.back());
			fichier >> Nprof;
			fichier >> Iprof >> Npoint ;
			
			cout << "Nprof" << "   " << Nprof << "   " << "Npoint" << "   " << Npoint << endl;	
			cout << "Re" << "   " << Re << endl;

			Array< OneD, Array<OneD, NekDouble> > u(Nprof);
			Array< OneD, Array<OneD, NekDouble> > v(Nprof);
			Array< OneD, Array<OneD, NekDouble> > t(Nprof);
			Array< OneD, Array<OneD, NekDouble> > y(Nprof);
			Array< OneD, Array<OneD, NekDouble> > rho(Nprof);
			Array< OneD, Array<OneD, NekDouble> > energy(Nprof);
			Array<OneD, NekDouble>  x(Nprof);
			
			cout << "c" << "   " << c1 << endl;
	
			for (int i=0; i< Nprof; i++)
			{
			u[i]= Array<OneD, NekDouble> (Npoint);	
			v[i]= Array<OneD, NekDouble> (Npoint);	
			t[i]= Array<OneD, NekDouble> (Npoint);	
			y[i]= Array<OneD, NekDouble> (Npoint);	
			rho[i]= Array<OneD, NekDouble> (Npoint);	
			energy[i]= Array<OneD, NekDouble> (Npoint);	
			}
			
			for (int i=0; i< Nprof; i++)
			{
				fichier >>  x[i] ;
				getline(fichier,monTableau.back());
				getline(fichier,monTableau.back());
				delta=(4.91/sqrt(Re))*sqrt(c1*x[i]);
				deltax << delta << endl;
			
				
				for (int j=0; j< Npoint ; j++)
				{
					fichier >> y[i][j] >> u[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> v[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> t[i][j] >> rho[i][j] >> energy[i][j];
					getline(fichier,monTableau.back());
					test << scientific << setprecision(9) << y[i][j] << "   "  ;
					y[i][j]=y[i][j]*delta;
				
					
					test << y[i][j] << endl;
				}
				
				fichier >> Iprof >> Iprof;		
			}
			
			for (int i=0; i< Nprof-1; i++)
			{
					if ((x[i]<=x_in)&&(x[i+1]>x_in))
					{
					indice=i;
					break;
					}
			}
			
			for (int j=0; j< Npoint; j++)
			{
				
			vitesse_u << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (u[indice][j]+u[indice+1][j])/2 << endl ;
			vitesse_v << scientific << setprecision(9) <<(y[indice][j]+y[indice+1][j])/2 << "   " << (v[indice][j]+v[indice+1][j])/2 << endl ;
			temperature << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (t[indice][j]+t[indice+1][j])/2 << endl ;
			rho_1 << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (rho[indice][j]+rho[indice+1][j])/2 << endl ;
			energie << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (energy[indice][j]+energy[indice+1][j])/2 << endl ;
			}
			
			
			fichier.close();		
// ---------------------------
// ---------------------------
// ---------------------------			

cout << "ca marche " << endl;
cout << y[Nprof-1][Npoint-1] << endl;


// Deuxième fichier
// ---------------------------
// ---------------------------
// ---------------------------

			int Nprof2,Npoint2;

			ofstream vitesse_u2;
			vitesse_u2.open("u2.dat");
	
			ofstream vitesse_v2;
			vitesse_v2.open("v2.dat");
			
			ofstream temperature2;
			temperature2.open("temperature2.dat");
		
			ofstream energie2;
			energie2.open("energy2.dat");
			
			ofstream rho_12;
			rho_12.open("rho2.dat");


			std::ifstream fichier2(argv[2], ios::in);
			
			getline(fichier2,monTableau.back());
			
			getline(fichier2,monTableau.back());
			fichier2 >> Nprof2;
			fichier2 >> Iprof >> Npoint2 ;


					
			cout << "Nprof2" << "   " << Nprof2 << "   " << "Npoint2" << "   " << Npoint2 << endl;	
			
			
			
			Array< OneD, Array<OneD, NekDouble> > u2(Nprof2);
			Array< OneD, Array<OneD, NekDouble> > v2(Nprof2);
			Array< OneD, Array<OneD, NekDouble> > t2(Nprof2);
			Array< OneD, Array<OneD, NekDouble> > y2(Nprof2);
			Array< OneD, Array<OneD, NekDouble> > rho2(Nprof2);
			Array< OneD, Array<OneD, NekDouble> > energy2(Nprof2);
			Array<OneD, NekDouble>  x2(Nprof2);
	
			for (int i=0; i< Nprof2; i++)
			{
			u2[i]= Array<OneD, NekDouble> (Npoint2);	
			v2[i]= Array<OneD, NekDouble> (Npoint2);	
			t2[i]= Array<OneD, NekDouble> (Npoint2);	
			y2[i]= Array<OneD, NekDouble> (Npoint2);	
			rho2[i]= Array<OneD, NekDouble> (Npoint2);	
			energy2[i]= Array<OneD, NekDouble> (Npoint2);	
			}
			
			for (int i=0; i< Nprof2; i++)
			{
				fichier2 >>  x2[i] ;
				getline(fichier2,monTableau.back());
				getline(fichier2,monTableau.back());
				delta2=(4.91/sqrt(Re))*sqrt(c1*x2[i]);
				
				for (int j=0; j< Npoint2 ; j++)
				{
					fichier2 >> y2[i][j] >> u2[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> v2[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> t2[i][j] >> rho2[i][j] >> energy2[i][j];
					getline(fichier2,monTableau.back());
					y2[i][j]=y2[i][j]*delta2;
	
				}
				
				fichier2 >> Iprof >> Iprof;		
			}
			
			for (int i=0; i< Nprof2-1; i++)
			{
					if ((x2[i]<=x_in)&&(x2[i+1]>x_in))
					{
					indice=i;
					break;
					}
			}
			
			for (int j=0; j< Npoint2; j++)
			{
				
			vitesse_u2 << scientific << setprecision(9) <<(y2[indice][j]+y2[indice+1][j])/2 << "   " << (u2[indice][j]+u2[indice+1][j])/2 << endl ;
			vitesse_v2 << scientific << setprecision(9) <<(y2[indice][j]+y2[indice+1][j])/2 << "   " << (v2[indice][j]+v2[indice+1][j])/2 << endl ;
			temperature2 << scientific << setprecision(9) <<(y2[indice][j]+y2[indice+1][j])/2 << "   " << (t2[indice][j]+t2[indice+1][j])/2 << endl ;
			rho_12 << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (rho2[indice][j]+rho2[indice+1][j])/2 << endl ;
			energie2 << scientific << setprecision(9) << (y[indice][j]+y[indice+1][j])/2 << "   " << (energy2[indice][j]+energy2[indice+1][j])/2 << endl ;
			}
			
			fichier2.close();		
// ---------------------------
// ---------------------------
// ---------------------------			

cout << "ca marche " << endl;

// Troisième fichier
// ---------------------------
// ---------------------------
// ---------------------------

			int Nprof3,Npoint3;

			ofstream vitesse_u3;
			vitesse_u3.open("u3.dat");
	
			ofstream vitesse_v3;
			vitesse_v3.open("v3.dat");
			
			ofstream temperature3;
			temperature3.open("temperature3.dat");
			
			ofstream test3;
			test3.open("test3.dat");


			std::ifstream fichier3(argv[3], ios::in);
			
			getline(fichier3,monTableau.back());
			
			getline(fichier3,monTableau.back());
			fichier3 >> Nprof3;
			fichier3 >> Iprof >> Npoint3 ;
			
					
			cout << "Nprof3" << "   " << Nprof3 << "   " << "Npoint3" << "   " << Npoint3 << endl;	
			
			
			cout << "c" << "   " << c1 << endl;
			
			Array< OneD, Array<OneD, NekDouble> > u3(Nprof3);
			Array< OneD, Array<OneD, NekDouble> > v3(Nprof3);
			Array< OneD, Array<OneD, NekDouble> > t3(Nprof3);
			Array< OneD, Array<OneD, NekDouble> > y3(Nprof3);
			Array<OneD, NekDouble>  x3(Nprof3);
	
			for (int i=0; i< Nprof3; i++)
			{
			u3[i]= Array<OneD, NekDouble> (Npoint3);	
			v3[i]= Array<OneD, NekDouble> (Npoint3);	
			t3[i]= Array<OneD, NekDouble> (Npoint3);	
			y3[i]= Array<OneD, NekDouble> (Npoint3);	
			}
			
			for (int i=0; i< Nprof3; i++)
			{
				fichier3 >>  x3[i] ;
				getline(fichier3,monTableau.back());
				getline(fichier3,monTableau.back());
				delta3=(1.00/sqrt(Re))*sqrt(c1*x3[i]);
				
				
				
				for (int j=0; j< Npoint3 ; j++)
				{
					fichier3 >> y3[i][j] >> u3[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> v3[i][j] >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> rubbish >> t3[i][j];
					getline(fichier3,monTableau.back());
					test3 << scientific << setprecision(9) << y3[i][j] << "   " ;
					y3[i][j]=y3[i][j]*delta3;
					v3[i][j]=(v3[i][j]);
					test3 << y3[i][j] << endl ;
				}
				
				fichier3 >> Iprof >> Iprof;		
			}
			
			for (int i=0; i< Nprof3-1; i++)
			{
					if ((x3[i]<=x_in)&&(x3[i+1]>x_in))
					{
					indice=i;
					break;
					}
			}
			
			for (int j=0; j< Npoint3; j++)
			{
				
			vitesse_u3 << scientific << setprecision(9) <<(y3[indice][j]+y3[indice+1][j])/2 << "   " << (u3[indice][j]+u3[indice+1][j])/2 << endl ;
			vitesse_v3 << scientific << setprecision(9) <<(y3[indice][j]+y3[indice+1][j])/2 << "   " << (v3[indice][j]+v3[indice+1][j])/2 << endl ;
			temperature3 << scientific << setprecision(9) <<(y3[indice][j]+y3[indice+1][j])/2 << "   " << (t3[indice][j]+t3[indice+1][j])/2 << endl ;
			}
			
			fichier3.close();	
			
cout << "ca marche " << endl;
cout << y3[Nprof3-1][Npoint3-1] << endl;
	
// ---------------------------
// ---------------------------
// ---------------------------	
	
			
							
    return 0;
}
    
    
