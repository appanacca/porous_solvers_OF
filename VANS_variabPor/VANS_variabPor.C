/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "kdtree.H"



#include "fvCFD.H"
#include "pisoControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // this are the numbers of realization in the database metamodel
    int nSamples = 1778498;
    // number of parameters in the metamodel
	int dimx = 4;
    // number of outputs, the 3 elements of the permeability tensor
	int dimy = 3;

    // matrix for the input parameters
	Eigen::Matrix<double, Dynamic, Dynamic>  mat(nSamples, dimx);
    // matrix for the output results
	Eigen::Matrix<double, Dynamic, Dynamic>  maty(nSamples, dimy);

    // name of the file in which the database is saved
	char* res = "result.dat";
	read_file(mat, maty, res);

	KDTreeEigenMatrixAdaptor< Eigen::Matrix<double,Dynamic,Dynamic>,nanoflann::metric_L2> mat_index(mat, 20 );

	kdtree_build<double>(mat_index, nSamples /* samples */, dimx /* dim */, dimy);

	// Query point:  points in wich we will ask for NN in the kdtree
	std::vector<double> query_pt(4);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;


        #include "CourantNo.H"

// this is a flag that "activate", or not, the metamodelling procedure
scalar flag_k_update = readScalar(transportProperties.lookup("k_update"));



if (flag_k_update == 1){ 
    forAll(K, j){
        // search through only the porous cells
        if( K[j].xx() != 0)
        {
            // this check is done for symmetry reason
            // beacuse in this case the 11 and 22 components has to be switched
            if(std::fabs(angle_zx[j]) > 45){
        		query_pt[1] = std::fabs(angle_zx[j]);
        		query_pt[2] = std::fabs(angle_yx[j]); 
        		query_pt[0] = reynolds[j]; 
        		query_pt[3] = por[j]; 

                // search in the database the best fit for the 4 parameters above
                double a = kdtree_search<double>(mat_index, query_pt);
    	        K[j].xx() = maty(a,1);
                K[j].yy() = maty(a,2);
                K[j].zz() = maty(a,0);
            }
            else{
        	    query_pt[1] = std::fabs(angle_zx[j]);  
            	query_pt[2] = std::fabs(angle_yx[j]);  
            	query_pt[0] = reynolds[j];     
            	query_pt[3] = por[j];  
    
                double a = kdtree_search<double>(mat_index, query_pt);      
        	    K[j].xx() = maty(a,0); 
                K[j].yy() = maty(a,2);  
                K[j].zz() = maty(a,1);  
            }   
        }
    }
}
else{Info << "Not updating K" << endl;}

// this two operations need to be done before the matrix assembly
// becouse the operators fvc::Sc
volVectorField KU = K & U;


        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + (1/por)*fvm::div(phi*phiPor, U)
            - fvm::laplacian(nu, U)
            - (nu/por)* (fvc::grad(U) & fvc::grad(por)) 
            // the grad is explicit because the implicit version is not present in OF
            - (nu/por)*fvc::laplacian(por) * U
	        + fvc::Sp(por*nu*(1/(LenRef*LenRef)),KU)
            - gravity 
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();

	/////////////////////////////////

	dimensionedScalar tolerance
	(
	"tolerance",
	dimensionSet(0,1,-1,0,0,0,0),
	SMALL
	); 


	angle_yx = atan(U.component(1) / (U.component(0) + tolerance)) *(180.0 /constant::mathematical::pi);
	angle_zx = atan(U.component(2) / (U.component(0) + tolerance)) *(180.0 / constant::mathematical::pi);
	reynolds = (mag(U)*LenRef)/nu;


	////////////////////////////////
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
