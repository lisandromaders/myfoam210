if (fgmTable.debugOutput())
{
    Info<< "### Write thermo property fields ###" << endl;
    volScalarField CpDebug
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.Cp()
    );
    CpDebug.write();

    volScalarField psiDebug
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.psi()
    );
    psiDebug.write();

    volScalarField muDebug
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.mu()
    );
    muDebug.write();

    volScalarField DtDebug
    (
        IOobject
        (
            "Dt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.Dt()
    );
    DtDebug.write();

    volScalarField kappaDebug
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.kappa()
    );
    kappaDebug.write();

    volScalarField rhoDebug
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );
    rhoDebug.write();
}
