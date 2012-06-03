///////////////////////////////////////////////////////////////////////////////
//
// File ExtractCriticalLayer.cpp
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
// Description: Extract location of critical layer from streak file 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2D.h>


void Computestreakpositions(MultiRegions::ExpListSharedPtr &streak,
                            Array<OneD, NekDouble> &xc,  
                            Array<OneD, NekDouble> &yc, 
                            NekDouble cr); 

int main(int argc, char *argv[])
{
    int i,j;
    NekDouble cr = 0;
    
    if(argc !=3)
    {
        fprintf(stderr,"Usage: ./ExtractCriticalLayer  meshfile fieldfile  \n");
        exit(1);
    }
    
    //------------------------------------------------------------
    // Create Session file. 
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    //-----------------------------------------------------------
    
    //-------------------------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    //------------------------------------------------------------

    //-------------------------------------------------------------
    // Define Streak Expansion   
    MultiRegions::ExpListSharedPtr streak;   

    streak = MemoryManager<MultiRegions::ExpList2D>
        ::AllocateSharedPtr(vSession,graphShPt);

    //---------------------------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    string  streak_field("w");
    for(unsigned int i = 0; i < fielddata.size(); ++i)
    {
        streak->ExtractDataToCoeffs(fielddef [i],
                                    fielddata[i],
                                    streak_field);
    }
    //----------------------------------------------
    
    int npts;
    vSession->LoadParameter("NoCriticalLayerPts",npts,30);
    Array<OneD, NekDouble> x_c(npts);
    Array<OneD, NekDouble> y_c(npts);       
    NekDouble xc, yc;

    Computestreakpositions(streak,x_c, y_c,cr);    

    cout << "# x_c y_c" << endl;
    for(i = 0; i < npts; ++i)
    {
        cout << x_c[i] << " " << y_c[i] << endl;
    }
}

void Computestreakpositions(MultiRegions::ExpListSharedPtr &streak,
                            Array<OneD, NekDouble> &xc,  
                            Array<OneD, NekDouble> &yc,
                            NekDouble cr)
{
    int i;
    int npts = xc.num_elements();

    int nq = streak->GetTotPoints();	
    Array<OneD, NekDouble> derstreak(nq);
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);    
    streak->GetCoords(x,y);         

    streak->BwdTrans(streak->GetCoeffs(),streak->UpdatePhys());
    streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), derstreak);

    // set intiial xc to be equispaced over mesh and yc to be zero 
    NekDouble x_max = Vmath::Vmax(nq,x,1);
    NekDouble x_min = Vmath::Vmin(nq,x,1);

    for(i = 0; i < npts; ++i)
    {
        xc[i]  = x_min + (x_max - x_min)*i/((NekDouble)(npts-1));
        yc[i] = 0.0;
    }
    
    
    int elmtid, offset,cnt;    
    NekDouble U,dU;
    NekDouble F;
    NekDouble ConvTol  = 1e-9;
    NekDouble CoordTol = 1e-5;
    int maxiter = 100;
    Array<OneD, NekDouble> coord(2);

    // Do Newton iteration on y direction
    cerr << "[";
    for(int e=0; e<npts; e++)
    {
        coord[0] = xc[e];
        coord[1] = yc[e];

        if(!(e%10))
        {
            cerr << ".";
        }

        F = 1000;
        cnt = 0;
        while((abs(F)> ConvTol)&&(cnt < maxiter))
        {
            elmtid = streak->GetExpIndex(coord,CoordTol);
            offset = streak->GetPhys_Offset(elmtid);

            U  = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
            dU = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);

            coord[1] = coord[1] - (U-cr)/dU;   
            
            F = U-cr;   
            cnt++;
        }
        ASSERTL0(cnt < maxiter, "Failed to converge Newton iteration");

        yc[e] = coord[1];
    }
    cerr << "]" << endl;
}
