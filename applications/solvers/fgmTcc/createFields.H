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


/*    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();
*/
    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    dimensionedScalar fanzy
    (
       	transportProperties.lookup("fanzy")
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

   volScalarField Dt
   (
    	IOobject
       (
           "Dt",
           runTime.timeName(),
           mesh,
           IOobject::READ_IF_PRESENT,
           IOobject::AUTO_WRITE
       ),
       thermo.Dt()
   );

    volScalarField kappa
    (
     	IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	thermo.kappa()
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

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
