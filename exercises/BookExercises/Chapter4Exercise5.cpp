#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

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
            // The mesh and boundary condition information are
            // contained in the input file Chapter4Exercise4_Quad_Pp.xml
            // The first step is to read in this file
            stringstream fileName;
            fileName << "Chapter4Exercise4_Quad_P" << p << ".xml";

            string fileNameString = fileName.str();

            SpatialDomains::MeshGraph2D mesh; 
            // Both the geometry and the expansion information should be read
            mesh.ReadGeometry(fileNameString);
            mesh.ReadExpansions(fileNameString);
            // Also read the boundary conditions
            SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
            boundaryConds.Read(fileNameString);

            // Construct an object from the class ContField2D.
            // This is the class which represents a multi-elemental
            // continuous spectral/hp expansion with boundary conditions.
            // This object can be constructed based on the input mesh
            // and the boundary conditions.
            MultiRegions::ContField2DSharedPtr multiElementExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(mesh,boundaryConds);

            // During the construction of this object, the mapping arrays from
            // local to global numbering system have been set up. They are all
            // contained in the data member m_locToGloMap. This variable is an object of the
            // class LocalToGlobalMap.
            // Also the expansions coefficients corresponding to the Dirichlet DOF's
            // have been calculated during the construction of ContField2D object.
            // The boundary conditions are stored in the ContField2D class, using an object
            // of the class ExpList1D, which is a one-dimensional spectral/hp expansion 
            // defined on the boundary of the domain.

            // CONSTRUCT THE GLOBAL MASS MATRIX
            // Step 1: construct the elemental mass matrix
            DNekScalBlkMatSharedPtr elementalMassMatrix;

            // As the elemental mass matrix will be identical for every element,
            // we will only generate this mass matrix for the first element
            // Note that we are generating the elemntal matrix in static condensed format
            LocalRegions::MatrixKey masskey(StdRegions::eMass,
                                            (multiElementExp->GetExp(0))->DetExpansionType(),
                                            *(multiElementExp->GetExp(0)));

            elementalMassMatrix = (multiElementExp->GetExp(0))->GetLocStaticCondMatrix(masskey);

            // Step 2: allocate space for the global mass matrix
            // We will solve the system using the static condensation technique. Therefore,
            // we will assemble the four different blocks required for this technique rather
            // than just one global matrix.

            // First, define some dimensions;
            int nElements         = multiElementExp->GetExpSize();
            int nDirDofs          = (multiElementExp->GetLocalToGlobalMap())->GetNumDirichletDofs();
            int nGlobalBndDofs    = (multiElementExp->GetLocalToGlobalMap())->GetTotGloBndDofs();
            int nGlobalIntDofs    = multiElementExp->getContNcoeffs() - nGlobalBndDofs;
            // As all elements are identical, we can calculate the dimensions below just once
            int nElementalBndDofs = (multiElementExp->GetExp(0))->NumBndryCoeffs();
            int nElementalDofs    = (multiElementExp->GetExp(0))->GetNcoeffs();

            // Allocate space for the Schur complement
            int dimSchurCompl = nGlobalBndDofs - nDirDofs;

            NekDouble zero = 0.0;
            DNekMatSharedPtr schurComplement = MemoryManager<DNekMat>::
                AllocateSharedPtr(dimSchurCompl,dimSchurCompl,zero);

            // Allocate space for the other blocks
            Array<OneD,unsigned int> bndBlockSize(nElements,nElementalBndDofs);
            Array<OneD,unsigned int> intBlockSize(nElements,nElementalDofs-nElementalBndDofs);
            DNekScalBlkMatSharedPtr BinvD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(bndBlockSize,intBlockSize);
            DNekScalBlkMatSharedPtr invD  = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(intBlockSize,intBlockSize);
            DNekScalBlkMatSharedPtr C     = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(intBlockSize,bndBlockSize);

            // Step 3: Assemble the global matrices from the elemental mass matrix
            int n,i,j;
            int cnt;
            int globalId1;
            int globalId2;
            int sign1;
            int sign2;
            DNekScalMatSharedPtr block;
            for(n = cnt = 0; n < nElements; ++n)
            {	  
                BinvD->SetBlock(n,n,block = elementalMassMatrix->GetBlock(0,1));
                invD->SetBlock(n,n,block = elementalMassMatrix->GetBlock(1,1));
                C->SetBlock(n,n,block = elementalMassMatrix->GetBlock(1,0));
                
                for(i = 0; i < nElementalBndDofs; ++i)
                {
                    globalId1 = (multiElementExp->GetLocalToGlobalMap())->GetBndMap(cnt + i);
                    sign1     = (multiElementExp->GetLocalToGlobalMap())->GetBndSign(cnt + i);
                    
                    if(globalId1 >= nDirDofs)
                    {                        
                        for(j = 0; j < nElementalBndDofs; ++j)
                        {
                            globalId2 = (multiElementExp->GetLocalToGlobalMap())->GetBndMap(cnt + j);
                            sign2     = (multiElementExp->GetLocalToGlobalMap())->GetBndSign(cnt + j);
                            
                            if(globalId2 >= nDirDofs)
                            {  
                                (*schurComplement)(globalId1-nDirDofs,globalId2-nDirDofs) 
                                    += sign1*sign2*(*elementalMassMatrix)(i,j);       
                            }           
                        }	
                    }	                
                }
                cnt += nElementalBndDofs;
            }    
        
            // CONSTRUCT THE RIGHT HAND SIDE VECTOR
            // Evaluate the forcing function at the quadrature points
            int nTotQuadPoints = multiElementExp->GetPointsTot();
            Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x3(nTotQuadPoints,0.0);
            multiElementExp->GetCoords(x1,x2);

            Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
            // Read the forcing function equation as defined in the input file
            SpatialDomains::ConstForcingFunctionShPtr forcingFunctionEquation 
                = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));
        
            for(i = 0; i < nTotQuadPoints; ++i)
            {
                forcingFunction[i] = forcingFunctionEquation->Evaluate(x1[i],x2[i],x3[i]);
            }

            // Take the inner product with all local elemental modes
            Array<OneD,NekDouble> elementalCoeffs;
            Array<OneD,NekDouble> elementalPhys;
        
            int cnt_phys  = 0;
            int cnt_coeff = 0;
            for(i = 0; i < nElements; ++i)
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

            // The right hand side, which is now contained in the member m_contCoeffs of the object
            // multiElementExp should now be modified in order to incorporate the forcing due to the 
            // Dirichlet Dofs.
            Array<OneD,NekDouble> dirDofsValues(multiElementExp->getContNcoeffs(),0.0);
            Array<OneD,NekDouble> dirForcing(multiElementExp->getContNcoeffs(),0.0);

            cnt = 0;
            for(i = 0; i < (multiElementExp->GetBndCondExp()).num_elements(); ++i)
            {
                for(j = 0; j < ((multiElementExp->GetBndCondExp())[i])->GetNcoeffs(); j++)
                {
                    dirDofsValues[(multiElementExp->GetLocalToGlobalMap())->GetBndCondMap(cnt++)] = 
                        (((multiElementExp->GetBndCondExp())[i])->GetCoeffs())[j];
                }
            }

            // Multiply the Dirichlet DOFs with the mass matrix to get the
            // rhs due to the dirichlet DOFs
            MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
            multiElementExp->GeneralMatrixOp(key, dirDofsValues, dirForcing);

            // Substract this from the original rhs
            Vmath::Vsub(multiElementExp->getContNcoeffs(), multiElementExp->UpdateContCoeffs(), 1,
                        dirForcing, 1, multiElementExp->UpdateContCoeffs(), 1);

            // Extract the rhs that is going to be used in the solution of the linear system.
            // This corresponds to the member m_contCoeffs, without the Dirichlet DOFs ( = the first
            // entries in this array). As rhs is a shared array, it will still point to the data
            // in m_contCoeffs.
            Array<OneD, NekDouble> rhs = multiElementExp->UpdateContCoeffs() + nDirDofs;

            // SOLVE THE LINEAR SYSTEM
            // This is done using the static condensation technique

            Array<OneD,NekDouble>  offset;  
            NekVector<NekDouble> rhsBndVector(dimSchurCompl,rhs,eCopy);
            NekVector<NekDouble> rhsIntVector(nGlobalIntDofs,rhs + dimSchurCompl,eCopy);
            NekVector<NekDouble> localBndDofsSolutionVector(dimSchurCompl,rhs,eWrapper);
            NekVector<NekDouble> localIntDofsSolutionVector(nGlobalIntDofs,offset = rhs + dimSchurCompl,eWrapper);
            NekVector<NekDouble> tmpVector;
                    
            // Construct boundary forcing 
            tmpVector = (*BinvD)*rhsIntVector;
            (multiElementExp->GetLocalToGlobalMap())->AssembleBnd(tmpVector,localBndDofsSolutionVector,nDirDofs);
            rhsBndVector = rhsBndVector - localBndDofsSolutionVector;
                        
            // Solve Schur Complement system  
            DNekLinSysSharedPtr schurComplementLinSys = 
                MemoryManager<DNekLinSys>::AllocateSharedPtr(schurComplement);
            schurComplementLinSys->Solve(rhsBndVector,localBndDofsSolutionVector);

            // Solve interior system 
            Vmath::Zero(tmpVector.GetDimension(),tmpVector.GetRawPtr(),1);
            (multiElementExp->GetLocalToGlobalMap())->ContToLocalBnd(localBndDofsSolutionVector,tmpVector,nDirDofs);
            rhsIntVector = rhsIntVector - (*C)*tmpVector;

            localIntDofsSolutionVector = (*invD)*rhsIntVector;

            // Recover the entire solution by addinig intial conditons
            Vmath::Zero(nDirDofs, multiElementExp->UpdateContCoeffs(), 1);
            Vmath::Vadd(multiElementExp->getContNcoeffs(), dirDofsValues, 1, multiElementExp->UpdateContCoeffs(), 1,
                multiElementExp->UpdateContCoeffs(), 1);

            // DO A BACKWARD TRANSFORMATION TO CALCULATE THE EXPANSION
            // EVALUATED AT THE QUADRATURE POINTS
            multiElementExp->ContToLocal();

            cnt_phys  = 0;
            cnt_coeff = 0;
            // The backward transformation is evaluated elementally
            for(i= 0; i < nElements; ++i)
            {
                elementalPhys   = multiElementExp->UpdatePhys() + cnt_phys;
                elementalCoeffs = multiElementExp->UpdateCoeffs() + cnt_coeff;
                (multiElementExp->GetExp(i))->BwdTrans(elementalCoeffs, 
                                                       elementalPhys);
            
                cnt_coeff  += (multiElementExp->GetExp(i))->GetNcoeffs();
                cnt_phys   += (multiElementExp->GetExp(i))->GetTotPoints();
            }  
            ofstream outfile("testie.dat");
                multiElementExp->WriteToFile(outfile,eTecplot);
            outfile.close();

            // CALCULATE THE L2 NORM OF THE ERROR
            NekDouble error = 0.0;
            cnt_phys  = 0;
            // the error is calculated elementally
            for(i= 0; i < nElements; ++i)
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
            // The mesh and boundary condition information are
            // contained in the input file Chapter4Exercise4_Quad_Pp.xml
            // The first step is to read in this file
            stringstream fileName;
            fileName << "Chapter4Exercise4_Quad_P" << p << ".xml";

            string fileNameString = fileName.str();

            SpatialDomains::MeshGraph2D mesh; 
            // Both the geometry and th expansion information should be read
            mesh.ReadGeometry(fileNameString);
            mesh.ReadExpansions(fileNameString);
            // Also read the boundary conditions
            SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
            boundaryConds.Read(fileNameString);

            // Construct an object from the class ContField2D.
            // This is the class which represents a multi-elemental
            // continuous spectral/hp expansion with boundary conditions.
            // This object can be constructed based on the input mesh
            // and the boundary conditions.
            MultiRegions::ContField2DSharedPtr multiElementExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(mesh,boundaryConds);

            // Evaluate the forcing function at the quadrature points
            int nTotQuadPoints = multiElementExp->GetPointsTot();
            Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x3(nTotQuadPoints,0.0);
            multiElementExp->GetCoords(x1,x2);

            Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
            // Read the forcing function equation as defined in the input file
            SpatialDomains::ConstForcingFunctionShPtr forcingFunctionEquation 
                = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));
        
            for(int i = 0; i < nTotQuadPoints; ++i)
            {
                forcingFunction[i] = forcingFunctionEquation->Evaluate(x1[i],x2[i],x3[i]);
            }

            // Store the forcing function as the physical values of an
            // object of the class ContExpList2D
            MultiRegions::ContField2DSharedPtr forcingExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*multiElementExp);
            forcingExp->SetPhys(forcingFunction);

            // Do the projection to obtain the coefficients of the expansion
            // The result is stored in the data member m_contCoeffs of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->FwdTrans(*forcingExp);

            // Perform a backward transformation to obtain the solution at the quadrature points
            // The result is stored in the data member m_phys of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->BwdTrans(*multiElementExp);

            // Calculate the error
            NekDouble error = multiElementExp->L2(*forcingExp);

            // Display the output              
            cout << "P = " << p << " => Error = " << error << endl;
            // In addition, we will use one the Nektar++ output formats to visualise the output.
            // The solution is written to the files Chapter4Exercise5_Quad_Pi.pos
            // which can be opened using the program Gmsh. Given the values of the coefficients of 
            // the expansion, this program then plots the expansion in a high-order fashion.  
            stringstream outfileName;
            outfileName << "Chapter4Exercise5_Quad_P" << p << ".pos";
            ofstream outfile((outfileName.str()).data());      
            multiElementExp->WriteToFile(outfile,eGmsh);
        }
        cout << "To see a plot of the solution, open the files Chapter4Exercise5_Quad_Pi.pos" <<endl;
        cout << "using the program Gmsh."<<endl;
        cout << endl;       
    }


}

