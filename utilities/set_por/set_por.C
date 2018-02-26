/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    scale_U

Description
    Scale the velocity field U by a given scalar.

   
\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include <math.h>       /* fabs */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*
    this function modfify the cell_por_value 
	the idea is to set a scalar value named z_itf that define the interface
	between the porous media and the free flow
*/
template<typename T1, typename T2> 
void itf_function(const double& itf_value, T1& cell_por_value, const T2& cell_height, const double& por_max, const double& len_ref){
	double len_cyl(0.), corr(0.);


	if( ((cell_height) <= ((itf_value +len_ref/2)) ) && ( (cell_height) >= ((itf_value -len_ref/2))) ){
            
    	len_cyl = std::fabs(((-cell_height +itf_value +len_ref*0.5)));
		corr = 1 -((1 -por_max)/ len_ref )*len_cyl;
		cell_por_value = corr;
    }
    else if( ((cell_height) < (itf_value -len_ref/2)) ){ 
        cell_por_value = por_max;
    }
    else if( ((cell_height) > (itf_value +len_ref/2)) ){
        cell_por_value = 1;
    }
}



void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");
       

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    dimensionedScalar LenRef
    (
        transportProperties.lookup("LenRef")
    );


    dimensionedScalar Itf
    (
        transportProperties.lookup("Itf")
    );

    dimensionedScalar porosity
    (
        transportProperties.lookup("porosity")
    );


    IOobject porheader
    (
        "por",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (porheader.headerOk())
    {
        Info<< "Reading por" << endl;
        volScalarField por(porheader, mesh);
        
        double por_max = porosity.value();

        forAll(por, j){

                /*
                    if you want another shape of the interface you need to implement another function like this one
                */
                itf_function(por[j], mesh.C()[j][1], Itf.value(), por_max, LenRef.value());

        }

        if (writeResults)
        {
            por.write();
        }
    }
    else
    {
        Info<< "No por" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
