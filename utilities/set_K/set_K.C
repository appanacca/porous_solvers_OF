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
template<typename T1, typename T2, typename T3> 
void itf_function(const double& itf_value, T1& cell_K_value, const T2& cell_height, T3& K_input_value, const double& len_ref){
    double len_cyl(0.), corr1(0.), corr2(0.), corr3(0.);

    if( ((cell_height) <= ((itf_value +len_ref/2)) ) && ( (cell_height) >= ((itf_value -len_ref/2))) ){
            
        len_cyl = std::fabs(((-cell_height +itf_value +len_ref*0.5)));

        corr1 = (K_input_value[0].value()/ len_ref )*len_cyl;
        corr2 = (K_input_value[1].value()/ len_ref )*len_cyl;
        corr3 = (K_input_value[2].value()/ len_ref )*len_cyl;
        cell_K_value.xx() = corr1;
        cell_K_value.yy() = corr2;
        cell_K_value.zz() = corr3;
    }
    else if( ((cell_height) < (itf_value -len_ref/2)) ){ 
        cell_K_value.xx() = K_input_value[0].value();
        cell_K_value.yy() = K_input_value[1].value();
        cell_K_value.zz() = K_input_value[2].value();
    }
    else if( ((cell_height) > (itf_value +len_ref/2) )){
        cell_K_value.xx() = 0;
        cell_K_value.yy() = 0;
        cell_K_value.zz() = 0;
    }
};


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

    dimensionedVector K_input
    (
        transportProperties.lookup("K_input")
    );

    dimensionedScalar Itf
    (
        transportProperties.lookup("Itf")
    );


    dimensionedScalar porosity
    (
        transportProperties.lookup("porosity")
    );



    IOobject Kheader
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Kheader.headerOk())
    {
        Info<< "Reading K" << endl;
        volTensorField K(Kheader, mesh);

        forAll(K, j){
            itf_function(Itf.value(), K[j], mesh.C()[j][1], K_input, LenRef.value());
        }
        


        if (writeResults)
        {
            K.write();
        } 

    }
      else
      {
        Info<< "No K" << endl;
     }

    Info<< "\nEnd\n" << endl;


}
// ************************************************************************* //
