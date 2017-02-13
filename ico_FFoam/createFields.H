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
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

Info<< "Reading nu_fluid\n" << endl;
dimensionedScalar nu_fluid
(
    "nu_fluid",
    dimViscosity,
    transportProperties.lookup("nu_fluid")
);

Info<< "Reading delta_ij\n" << endl;
    dimensionedVector delta_ij
    (
        transportProperties.lookup("delta_ij")
    );



Info<< "Reading field d\n" << endl;
volScalarField d
(
    IOobject
    (
        "d",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
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


//#include "createPhi.H"
// the following do the same thing in "createPhi.H", but here I can change the name of the vector field
Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(D) & mesh.Sf()
);

surfaceScalarField phiU
(
    IOobject
    (
        "phiU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);




label dRefCell = 0;
scalar dRefValue = 0.0;
setRefCell(d, mesh.solutionDict().subDict("PISO"), dRefCell, dRefValue);
mesh.setFluxRequired(d.name());
