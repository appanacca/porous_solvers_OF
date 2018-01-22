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

                double z_itf =0;

                if(  ((mesh.C()[j][0]) >= 0) &&  ((mesh.C()[j][0]) <= 0.3214)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*std::min(28., (28. +0.0067751*x*x -0.0021245*x*x*x) );
                }
                else if(  ((mesh.C()[j][0]) >= 0.3214) &&  ((mesh.C()[j][0]) <= 0.5)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(25.073 +0.97548*x -0.10161*x*x +0.0018898*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) >= 0.5) &&  ((mesh.C()[j][0]) <= 0.7143)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(25.796 +0.82067*x -0.090554*x*x +0.0016265*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) >= 0.7143) &&  ((mesh.C()[j][0]) <= 1.0714)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(40.464 -1.3796*x +0.019459*x*x -0.00020703*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) >= 1.0714) &&  ((mesh.C()[j][0]) <= 1.4286)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(17.925 +0.87439*x -0.055674*x*x +0.00062777*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) >= 1.4286) &&  ((mesh.C()[j][0]) <= 1.9286)   ){
                    double x = mesh.C()[j][0]*28;
                    z_itf = (1./28.)*std::max(0., (56.390 -2.0105*x +0.016449*x*x +0.000026749*x*x*x) );
                }
                else if(  ((mesh.C()[j][0]) <= 9) &&  ((mesh.C()[j][0]) >= (9-0.3214))   ){
                    double x = 9-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*std::min(28., (28. +0.0067751*x*x -0.0021245*x*x*x) );
                }
                else if(  ((mesh.C()[j][0]) <= (9-0.3214) ) &&  ((mesh.C()[j][0]) >= (9 -0.5))   ){
                    double x = 9*28-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(25.073 +0.97548*x -0.10161*x*x +0.0018898*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) <= (9 -0.5)) &&  ((mesh.C()[j][0]) >= (9-0.7143))   ){
                    double x = 9*28-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(25.796 +0.82067*x -0.090554*x*x +0.0016265*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) <= (9-0.7143)) &&  ((mesh.C()[j][0]) >= (9-1.0714))   ){
                    double x = 9*28-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(40.464 -1.3796*x +0.019459*x*x -0.00020703*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) <= (9-1.0714)) &&  ((mesh.C()[j][0]) >= (9-1.4286))   ){
                    double x = 9*28-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*(17.925 +0.87439*x -0.055674*x*x +0.00062777*x*x*x);
                }
                else if(  ((mesh.C()[j][0]) <= (9-1.4286)) &&  ((mesh.C()[j][0]) >= (9-1.9286))   ){
                    double x = 9*28-mesh.C()[j][0]*28;
                    z_itf = (1./28.)*std::max(0., (56.390 -2.0105*x +0.016449*x*x +0.000026749*x*x*x) );
                }
               
                if( ((mesh.C()[j][1]) <= ((z_itf +LenRef.value()/2)) ) && ( (mesh.C()[j][1]) >= ((z_itf -LenRef.value()/2))) ){
            
                    double len_cyl = std::fabs(((-mesh.C()[j][1] +z_itf +LenRef.value()*0.5)));

                    double corr1 = (K_input[0].value()/ LenRef.value() )*len_cyl;
                    double corr2 = (K_input[1].value()/ LenRef.value() )*len_cyl;
                    double corr3 = (K_input[2].value()/ LenRef.value() )*len_cyl;
                        K[j].xx() = corr1;
                        K[j].yy() = corr2;
                        K[j].zz() = corr3;
                    }
                else if( ((mesh.C()[j][1]) < (z_itf -LenRef.value()/2)) ){ 
                     K[j].xx() = K_input[0].value();
                     K[j].yy() = K_input[1].value();
                     K[j].zz() = K_input[2].value();
                }
                else if( ((mesh.C()[j][1]) > (z_itf +LenRef.value()/2) )){
                      K[j].xx() = 0;
                      K[j].yy() = 0;
                      K[j].zz() = 0;
                }
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
