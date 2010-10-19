#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 4 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "--- EXERCISE 4 ----" <<endl;
    cout << "-------------------" <<endl;

    // We will not solve this exercise without using the Nektar++ library

    cout << "METHOD 1: low-level Nektar++" << endl;
    { 
        for(int p = 4; p < 11; p++)
        {
            // The mesh is contained in the input file Chapter4Exercise4_Quad_Pp.xml
            // The first step is to read in this file
            stringstream fileName;
            fileName << "Chapter4Exercise4_Quad_P" << p << ".xml";

            string fileNameString = fileName.str();

            SpatialDomains::MeshGraph2D mesh; 
            // Both the geometry and the expansion information should be read
            mesh.ReadGeometry(fileNameString);
            mesh.ReadExpansions(fileNameString);

            // Construct an object from the class ContField2D.
            // This is the class which represents a multi-elemental
            // continuous spectral/hp expansion.
            // This object can be constructed based on the input mesh
            MultiRegions::ContField2DSharedPtr multiElementExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(mesh);

            // During the construction of this object, the mapping arrays from
            // local to global numbering system have been set up. They are all
            // contained in the data member m_locToGloMap. This variable is an object of the
            // class LocalToGlobalMap.

            // CONSTRUCT THE GLOBAL MASS MATRIX
            // Step 1: construct the elemental mass matrix
            DNekScalMatSharedPtr elementalMassMatrix;

            // As the elemental mass matrix will be identical for every element,
            // we will only generate this mass matrix for the first element
            LocalRegions::MatrixKey masskey(StdRegions::eMass,
                                            (multiElementExp->GetExp(0))->DetExpansionType(),
                                            *(multiElementExp->GetExp(0)));

            elementalMassMatrix = (multiElementExp->GetExp(0))->GetLocMatrix(masskey);

            // Step 2: allocate space for the global mass matrix
            int nGlobalDof = (multiElementExp->GetLocalToGlobalMap())->GetNumGlobalCoeffs();

            NekDouble zero = 0.0;
            DNekMatSharedPtr globalMassMatrix = MemoryManager<DNekMat>::
                AllocateSharedPtr(nGlobalDof,nGlobalDof,zero);

            // Step 3: Assemble the global mass matrix from the elemental mass matrix
            int n,i,j;
            int cnt;
            int globalId1;
            int globalId2;
            int sign1;
            int sign2;
            for(n = cnt = 0; n < multiElementExp->GetExpSize(); ++n)
            {		    
                for(i = 0; i < elementalMassMatrix->GetColumns(); ++i)
                {
                    globalId1 = (multiElementExp->GetLocalToGlobalMap())->GetLocalToGlobalMap(cnt + i);
                    sign1     = (multiElementExp->GetLocalToGlobalMap())->GetLocalToGlobalSign(cnt + i);
                
                    for(j = 0; j < elementalMassMatrix->GetColumns(); ++j)
                    {
                        globalId2 = (multiElementExp->GetLocalToGlobalMap())->GetLocalToGlobalMap(cnt + j);
                        sign2     = (multiElementExp->GetLocalToGlobalMap())->GetLocalToGlobalSign(cnt + j);

                        (*globalMassMatrix)(globalId1,globalId2) 
                            += sign1*sign2*(*elementalMassMatrix)(i,j);                  
                    }		                
                }
                cnt += (multiElementExp->GetExp(n))->GetNcoeffs();
            }    

        
            // CONSTRUCT THE RIGHT HAND SIDE VECTOR
            // Evaluate the forcing function at the quadrature points
            int nTotQuadPoints = multiElementExp->GetTotPoints();
            Array<OneD,NekDouble> x1(nTotQuadPoints);
            Array<OneD,NekDouble> x2(nTotQuadPoints);
            multiElementExp->GetCoords(x1,x2);

            Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
        
            for(i = 0; i < nTotQuadPoints; ++i)
            {
                forcingFunction[i] = cos(x1[i])*cos(x2[i]);
            }

            // Take the inner product with all local elemental modes
            Array<OneD,NekDouble> elementalCoeffs;
            Array<OneD,NekDouble> elementalPhys;
        
            int cnt_phys  = 0;
            int cnt_coeff = 0;
            for(i = 0; i < multiElementExp->GetExpSize(); ++i)
            {
                elementalPhys   = forcingFunction + cnt_phys;
                elementalCoeffs = multiElementExp->UpdateCoeffs() + cnt_coeff;
                (multiElementExp->GetExp(i))->IProductWRTBase(elementalPhys,
                                                              elementalCoeffs);

                cnt_phys  += (multiElementExp->GetExp(i))->GetTotPoints();
                cnt_coeff += (multiElementExp->GetExp(i))->GetNcoeffs();
            }

            // Assemble these coefficients to get the inner product with
            // all global modes
            multiElementExp->Assemble(multiElementExp->GetCoeffs(),
                                      multiElementExp->UpdateContCoeffs());

            // SOLVE THE LINEAR SYSTEM
            DNekLinSysSharedPtr globalMassMatrixLinSys = 
                MemoryManager<DNekLinSys>::AllocateSharedPtr(globalMassMatrix);

            NekVector<const NekDouble> in(multiElementExp->GetContNcoeffs(), 
                                          multiElementExp->UpdateContCoeffs(), 
                                          eCopy);
            NekVector<NekDouble>       out(multiElementExp->GetContNcoeffs(), 
                                           multiElementExp->UpdateContCoeffs(), 
                                           eWrapper);
            globalMassMatrixLinSys->Solve(in,out);

            // DO A BACKWARD TRANSFORMATION TO CALCULATE THE EXPANSION
            // EVALUATED AT THE QUADRATURE POINTS
            multiElementExp->GlobalToLocal();

            cnt_phys  = 0;
            cnt_coeff = 0;
            // The backward transformation is evaluated elementally
            for(i= 0; i < multiElementExp->GetExpSize(); ++i)
            {
                elementalPhys   = multiElementExp->UpdatePhys() + cnt_phys;
                elementalCoeffs = multiElementExp->UpdateCoeffs() + cnt_coeff;
                (multiElementExp->GetExp(i))->BwdTrans(elementalCoeffs, 
                                                       elementalPhys);
            
                cnt_coeff  += (multiElementExp->GetExp(i))->GetNcoeffs();
                cnt_phys   += (multiElementExp->GetExp(i))->GetTotPoints();
            }  

            // CALCULATE THE L2 NORM OF THE ERROR
            NekDouble error = 0.0;
            cnt_phys  = 0;
            // the error is calculated elementally
            for(i= 0; i < multiElementExp->GetExpSize(); ++i)
            {
                (multiElementExp->GetExp(i))->
                    SetPhys(multiElementExp->UpdatePhys()+cnt_phys);
            
                error += pow( (multiElementExp->GetExp(i))->
                              L2(forcingFunction+cnt_phys) , 2);

                cnt_phys   += (multiElementExp->GetExp(i))->GetTotPoints();
            }            
            error = sqrt(error);

            // Display the output              
            cout << "P = " << p << " => Error = " << error << endl;         
        }
        cout << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        for(int p = 4; p < 11; p++)
        {
            // The mesh is contained in the input file Chapter4Exercise4_Quad_Pp.xml
            // The first step is to read in this file
            stringstream fileName;
            fileName << "Chapter4Exercise4_Quad_P" << p << ".xml";

            string fileNameString = fileName.str();

            SpatialDomains::MeshGraph2D mesh; 
            // Both the geometry and th expansion information should be read
            mesh.ReadGeometry(fileNameString);
            mesh.ReadExpansions(fileNameString);

            // Construct an object from the class ContExpList2D.
            // This is the class which represents a multi-elemental
            // continuous spectral/hp expansion.
            // This object can be constructed based on the input mesh
            MultiRegions::ExpList2DSharedPtr multiElementExp =
                MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(mesh);

            // Evaluate the forcing function at the quadrature points
            int nTotQuadPoints = multiElementExp->GetTotPoints();
            Array<OneD,NekDouble> x1(nTotQuadPoints);
            Array<OneD,NekDouble> x2(nTotQuadPoints);
            multiElementExp->GetCoords(x1,x2);

            Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
        
            for(int i = 0; i < nTotQuadPoints; ++i)
            {
                forcingFunction[i] = cos(x1[i])*cos(x2[i]);
            }

            // Store the forcing function as the physical values of an
            // object of the class ContExpList2D
            MultiRegions::ExpList2DSharedPtr forcingExp =
                MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*multiElementExp);
            forcingExp->SetPhys(forcingFunction);

            // Do the projection to obtain the coefficients of the expansion
            // The result is stored in the data member m_contCoeffs of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->FwdTrans(forcingExp->GetPhys(),multiElementExp->UpdateCoeffs());

            // Perform a backward transformation to obtain the solution at the quadrature points
            // The result is stored in the data member m_phys of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->BwdTrans(multiElementExp->GetCoeffs(),multiElementExp->UpdatePhys());

            // Calculate the error
            NekDouble error = multiElementExp->L2(forcingExp->GetPhys());

            // Display the output              
            cout << "P = " << p << " => Error = " << error << endl;  
            // In addition, we will use one the Nektar++ output formats to visualise the output.
            // The solution is written to the files Chapter4Exercise4_Quad_Pi.pos
            // which can be opened using the program Gmsh. Given the values of the coefficients of 
            // the expansion, this program then plots the expansion in a high-order fashion.  
            stringstream outfileName;
            outfileName << "Chapter4Exercise4_Quad_P" << p << ".pos";
            ofstream outfile((outfileName.str()).data());      
            multiElementExp->WriteToFile(outfile,eGmsh);     
        }
        cout << "To see a plot of the solution, open the files Chapter4Exercise4_Quad_Pi.pos" <<endl;
        cout << "using the program Gmsh."<<endl;
        cout << endl;       
    }


}

