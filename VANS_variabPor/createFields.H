    Info<< "Reading transportProperties\n" << endl;

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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar LenRef
    (
        transportProperties.lookup("LenRef")
    );




Info<< "Reading field gravity\n" << endl;
    dimensionedVector gravity
    (
        transportProperties.lookup("gravity")
    );


    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field porosity\n" << endl;
    volScalarField por
    (
        IOobject
        (
            "por",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field K\n" << endl;
    volTensorField K
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );




    Info<< "Reading field angle_yx\n" << endl;
    volScalarField angle_yx
    (
        IOobject
        (
            "angle_yx",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field angl_zx\n" << endl;
    volScalarField angle_zx
    (
        IOobject
        (
            "angle_zx",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field reynolds\n" << endl;
    volScalarField reynolds
    (
        IOobject
        (
            "reynolds",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );





#   include "createPhi.H"

 Info<< "Reading/calculating face flux field phi= U*por\n" << endl;
 
 surfaceScalarField phiPor
 (
     IOobject
     (
         "phiPor",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     linearInterpolate(por) //* mesh.Sf()
 );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
