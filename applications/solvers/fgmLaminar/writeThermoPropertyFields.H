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

    volScalarField alphaDebug
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.alpha()
    );
    alphaDebug.write();
}
